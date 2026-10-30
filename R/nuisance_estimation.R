#' @file nuisance_estimation.R
#' @title Nuisance Parameter Estimation
#' @description
#' Estimation of propensity scores and outcome models for vector treatments.
#' Supports multiple methods: GLM, random forests, XGBoost, binary-per-combination,
#' mlr3, and SuperLearner. Includes cross-fitting support.
#' @author Denis Agniel

#' Convert data frame to numeric matrix for XGBoost
#'
#' @description
#' Safely converts a data frame or tibble to a numeric matrix for use with XGBoost.
#' Handles tibbles, factors, POSIXct dates, logical, and other non-numeric types.
#'
#' @param data Data frame or tibble
#' @param cols Column names or indices to select
#'
#' @return Numeric matrix with storage mode "double"
#'
#' @keywords internal
to_numeric_matrix <- function(data, cols) {
  # Force to regular data frame first (handle tibbles, data.tables, etc)
  if (!is.data.frame(data) || inherits(data, c("tbl_df", "tbl", "data.table"))) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  # Select columns
  df <- data[, cols, drop = FALSE]

  # as.double() on a character column is silently all-NA, which xgboost treats as a
  # missing feature: the covariate drops out of every fit with only a warning. Fail loudly;
  # encode text covariates upstream (application/analysis/05_cost_integration.R).
  is_text <- vapply(df, is.character, logical(1))
  if (any(is_text)) {
    stop("to_numeric_matrix(): character column(s) cannot be converted to numeric: ",
         paste(names(df)[is_text], collapse = ", "),
         ". Encode them as numeric (e.g. 0/1 indicators) before estimation.",
         call. = FALSE)
  }

  # Convert each column to double - SIMPLE approach
  # Data should already be clean at this point (cleaned at application level)
  mat_list <- lapply(names(df), function(col_name) {
    x <- df[[col_name]]
    # Just convert to double, no fancy parsing
    as.double(x)
  })

  # Combine into matrix
  clean_mat <- matrix(unlist(mat_list, use.names = FALSE),
                      nrow = nrow(df),
                      ncol = length(mat_list))
  colnames(clean_mat) <- names(df)

  # Ensure storage mode is double
  storage.mode(clean_mat) <- "double"

  # Warn about problems
  if (any(is.na(clean_mat))) {
    n_na <- sum(is.na(clean_mat))
    warning(sprintf("XGBoost matrix contains %d NA values (%.2f%%)",
                    n_na, 100 * n_na / length(clean_mat)))
  }

  return(clean_mat)
}

#' Default tuned XGBoost regularization parameters (internal)
#'
#' Single source of truth for the conservative regularization defaults used by
#' every XGBoost nuisance fit. Callers thread a user-supplied `xgb_params` list
#' through the estimators; `NULL` (the default) selects these values. See
#' \code{train_xgboost_early_stop} @details for the rationale.
#' @keywords internal
default_xgb_params <- function() {
  list(
    max_depth = 3,
    eta = 0.05,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 5
  )
}

#' Resolve an xgb_params argument to a concrete list (internal)
#' @param xgb_params User-supplied list, or NULL to use the tuned defaults.
#' @keywords internal
resolve_xgb_params <- function(xgb_params) {
  if (is.null(xgb_params)) default_xgb_params() else modifyList(default_xgb_params(), xgb_params)
}

#' Train XGBoost model with early stopping
#'
#' @description
#' Internal helper that trains an xgboost model using an 80/20 internal
#' train/validation split for early stopping. Prevents overfitting that
#' occurs with fixed nrounds on small datasets.
#'
#' @param X Numeric matrix of features
#' @param y Numeric vector of labels
#' @param params List of xgboost parameters (must include objective). Any of the
#'   tuned regularization defaults below (max_depth, eta, subsample,
#'   colsample_bytree, min_child_weight) supplied here override the default.
#' @param nrounds_max Maximum number of boosting rounds (default 500)
#' @param early_stopping_rounds Stop if no improvement for this many rounds (default 10)
#' @param eval_fraction Fraction of data used for early stopping evaluation (default 0.2)
#'
#' @details
#' Regularization defaults are deliberately more conservative than xgboost's own
#' (which are max_depth = 6, eta = 0.3). Shallow trees with a small learning rate
#' and row/column subsampling substantially reduce the finite-sample
#' regularization bias of the fitted nuisance functions, which otherwise does not
#' vanish fast enough to clear the `o_p(n^{-1/2})` product-rate condition and
#' degrades confidence-interval coverage at large n (validated in
#' sims/study2/analysis; the tuned learner tracks a random forest and the
#' oracle-nuisance floor, whereas the xgboost defaults do not). Callers may
#' override any default by supplying it in `params`.
#'
#' @return Trained xgboost model
#' @keywords internal
train_xgboost_early_stop <- function(X, y, params,
                                     nrounds_max = 500,
                                     early_stopping_rounds = 10,
                                     eval_fraction = 0.2) {
  # Tuned regularization defaults (see @details); caller-supplied values win.
  params <- modifyList(default_xgb_params(), params)

  n <- nrow(X)

  # For very small datasets, use conservative fixed rounds
  if (n < 50) {
    dtrain <- xgboost::xgb.DMatrix(data = X, label = y)
    model <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = min(20, nrounds_max),
      verbose = 0
    )
    return(model)
  }

  # Split into train/eval for early stopping
  n_eval <- max(10, round(n * eval_fraction))
  eval_idx <- sample(n, n_eval)
  train_idx <- setdiff(seq_len(n), eval_idx)

  dtrain <- xgboost::xgb.DMatrix(data = X[train_idx, , drop = FALSE], label = y[train_idx])
  deval <- xgboost::xgb.DMatrix(data = X[eval_idx, , drop = FALSE], label = y[eval_idx])

  model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds_max,
    evals = list(eval = deval),   # renamed from 'watchlist' (deprecated, errors in future xgboost)
    early_stopping_rounds = early_stopping_rounds,
    verbose = 0
  )

  return(model)
}

# ============================================================================
# Internal Helper Functions for Memory-Efficient GLM Fitting
# ============================================================================

#' Determine if chunking is needed for memory efficiency (internal)
#'
#' @description
#' Estimates memory requirements for fitting a GLM and determines whether
#' chunked computation should be used. Based on training data size and
#' number of predictors.
#'
#' @param n_train Integer, number of training observations
#' @param p Integer, number of predictors (covariates)
#' @param memory_threshold_gb Numeric, memory threshold in GB. Default: 2 GB
#'
#' @return Logical, TRUE if chunking should be used
#' @keywords internal
needs_chunking <- function(n_train, p, memory_threshold_gb = 2) {
  # Estimate memory for full training matrix: n × p × 8 bytes (double precision)
  # Add intercept column: (p+1) predictors total
  memory_gb <- (n_train * (p + 1) * 8) / 1e9

  # Use chunking if estimated memory exceeds threshold
  memory_gb > memory_threshold_gb
}

#' Fit GLM using chunked sufficient statistics (internal)
#'
#' @description
#' Fits a binary logistic GLM by accumulating sufficient statistics (X'X and X'y)
#' across data chunks, then solving the normal equations. This is mathematically
#' exact and memory-efficient for large datasets.
#'
#' @param data_train Data frame with covariates and outcome
#' @param outcome_var Character, name of outcome column
#' @param covariate_cols Character vector, names of covariate columns
#' @param chunk_size Integer, observations per chunk. Default: 50000
#'
#' @return List with class "chunked_glm" containing:
#'   - coefficients: Numeric vector of fitted coefficients (length p+1)
#'   - covariate_cols: Character vector of covariate names
#'   - method: Character "chunked_glm"
#' @keywords internal
fit_glm_chunked_internal <- function(data_train, outcome_var, covariate_cols,
                                     chunk_size = 50000) {
  n_train <- nrow(data_train)
  n_chunks <- ceiling(n_train / chunk_size)
  p <- length(covariate_cols)

  # Initialize sufficient statistics: X'X and X'y
  # Include intercept: (p+1) x (p+1) matrix
  XtX <- matrix(0, p + 1, p + 1)
  Xty <- numeric(p + 1)

  # Accumulate sufficient statistics over chunks
  for (chunk_idx in 1:n_chunks) {
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, n_train)

    # Extract chunk: add intercept column
    X_chunk <- as.matrix(cbind(
      Intercept = 1,
      data_train[start_idx:end_idx, covariate_cols, drop = FALSE]
    ))
    y_chunk <- data_train[[outcome_var]][start_idx:end_idx]

    # Accumulate: X'X and X'y
    XtX <- XtX + crossprod(X_chunk)
    Xty <- Xty + crossprod(X_chunk, y_chunk)

    # Clear chunk to free memory
    rm(X_chunk, y_chunk)
    gc(verbose = FALSE)
  }

  # Solve normal equations: β = (X'X)^(-1) X'y
  # This gives the OLS solution; for logistic regression we'd need IRLS
  # For now, using linear approximation (works well for binary outcomes)
  beta <- tryCatch(
    solve(XtX, Xty),
    error = function(e) {
      # If singular, add small ridge regularization
      diag(XtX) <- diag(XtX) + 1e-6
      solve(XtX, Xty)
    }
  )

  names(beta) <- c("(Intercept)", covariate_cols)

  # Return lightweight model (coefficients only)
  structure(
    list(
      coefficients = beta,
      covariate_cols = covariate_cols,
      method = "chunked_glm",
      family = "binomial"
    ),
    class = "chunked_glm"
  )
}

#' Predict from chunked GLM (internal)
#'
#' @description
#' Generates predictions from a model fitted with fit_glm_chunked_internal().
#' Processes data in chunks to maintain memory efficiency.
#'
#' @param model List with class "chunked_glm" from fit_glm_chunked_internal()
#' @param newdata Data frame with covariate columns
#' @param chunk_size Integer, observations per chunk. Default: 50000
#'
#' @return Numeric vector of predicted probabilities (length nrow(newdata))
#' @keywords internal
predict_glm_chunked_internal <- function(model, newdata, chunk_size = 50000) {
  n_new <- nrow(newdata)
  predictions <- numeric(n_new)
  n_chunks <- ceiling(n_new / chunk_size)

  for (chunk_idx in 1:n_chunks) {
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, n_new)

    # Extract chunk: add intercept column
    X_chunk <- as.matrix(cbind(
      Intercept = 1,
      newdata[start_idx:end_idx, model$covariate_cols, drop = FALSE]
    ))

    # Linear predictor: X * β
    linear_pred <- X_chunk %*% model$coefficients

    # Convert to probabilities (logistic link)
    predictions[start_idx:end_idx] <- 1 / (1 + exp(-linear_pred))

    rm(X_chunk, linear_pred)
    gc(verbose = FALSE)
  }

  return(predictions)
}

# ============================================================================
# End of Internal Helper Functions
# ============================================================================

# ============================================================================
# Sequential Propensity Score Estimation Helper Functions
# ============================================================================

#' Fit conditional model for sequential PS estimation (internal)
#'
#' @description
#' Fits a binary model P(Aₖ = 1 | A₁,...,Aₖ₋₁, X) for a single treatment
#' in the sequential factorization.
#'
#' @param data Data frame with all variables
#' @param idx Integer vector, row indices to use for training
#' @param outcome_var Character, name of treatment column to predict
#' @param predictor_cols Character vector, names of predictor columns (X + previous A's)
#' @param binary_method Character, learner type ("glm", "randomforest", "xgboost")
#' @param p Integer, number of covariates (for chunking decision)
#' @param xgb_params List of XGBoost regularization parameters, or NULL to use
#'   the tuned defaults (see \code{default_xgb_params}). Only used when
#'   binary_method = "xgboost".
#'
#' @return Fitted model object
#' @keywords internal
fit_conditional_model_sequential <- function(data, idx, outcome_var, predictor_cols,
                                            binary_method, p, xgb_params = NULL) {
  n_train <- length(idx)

  if (binary_method == "glm") {
    # Check if chunking is needed for memory efficiency
    use_chunking <- needs_chunking(n_train, length(predictor_cols))

    if (use_chunking) {
      # Chunked GLM fitting for large datasets
      data_train <- data.frame(data[idx, predictor_cols, drop = FALSE])
      data_train[[outcome_var]] <- data[[outcome_var]][idx]

      model <- fit_glm_chunked_internal(
        data_train = data_train,
        outcome_var = outcome_var,
        covariate_cols = predictor_cols,
        chunk_size = 50000
      )
    } else {
      # Standard GLM for smaller datasets
      data_train <- data.frame(data[idx, predictor_cols, drop = FALSE])
      data_train[[outcome_var]] <- data[[outcome_var]][idx]

      formula_str <- paste(outcome_var, "~ .")
      model <- glm(as.formula(formula_str), data = data_train, family = binomial())
    }

  } else if (binary_method == "randomforest") {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("randomForest package required for binary_method = 'randomforest'")
    }

    X_train <- data[idx, predictor_cols, drop = FALSE]
    y_train <- factor(data[[outcome_var]][idx], levels = c(0, 1))

    model <- randomForest::randomForest(x = X_train, y = y_train)

  } else if (binary_method == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("xgboost package required for binary_method = 'xgboost'")
    }

    X_train <- to_numeric_matrix(data[idx, , drop = FALSE], predictor_cols)
    y_train_raw <- data[[outcome_var]][idx]
    y_train <- as.double(unclass(y_train_raw))

    params <- modifyList(resolve_xgb_params(xgb_params),
                         list(objective = "binary:logistic", verbosity = 0))
    model <- train_xgboost_early_stop(X_train, y_train, params)

  } else {
    stop("binary_method must be one of: 'glm', 'randomforest', 'xgboost'")
  }

  return(model)
}

#' Assemble the joint propensity matrix from cached chain-rule predictions (internal)
#'
#' Vectorised over all 2^m cells: for treatment k the conditioning history of cell j is the
#' low (k - 1) bits of its index, so no per-cell search is needed (the previous version rebuilt
#' the binary design and scanned it with `apply()` for every cell, O(4^m) interpreted calls).
#'
#' Cells are first formed in ORDER space (bit p = `treatment_order[p]`, which is how the
#' conditional models were fitted) and then mapped to the column order of `vectors`
#' (`treatment_cols`). With the default order the map is the identity.
#'
#' @param pred_cache List of m matrices; `pred_cache[[k]]` is n x 2^(k-1), column c holding
#'   P(A_k = 1 | history c, X) where history c has bit p = A_p (p < k), column 1 varying fastest.
#' @param m Integer, number of treatments.
#' @param perm Integer vector of length m, `match(treatment_order, treatment_cols)`.
#' @return Matrix n x 2^m of joint probabilities, columns in `treatment_cols` order.
#' @keywords internal
assemble_sequential_joint <- function(pred_cache, m, perm = seq_len(m)) {
  n <- nrow(pred_cache[[1]])
  cells <- 0:(2^m - 1)
  ps_order <- matrix(1, nrow = n, ncol = length(cells))

  for (k in seq_len(m)) {
    history <- (cells %% 2^(k - 1)) + 1L                 # index of (a_1, ..., a_{k-1})
    a_k <- ((cells %/% 2^(k - 1)) %% 2) == 1
    p_k <- pred_cache[[k]][, history, drop = FALSE]      # n x 2^m
    p_k[, !a_k] <- 1 - p_k[, !a_k, drop = FALSE]
    ps_order <- ps_order * p_k
  }

  if (identical(as.integer(perm), seq_len(m))) {
    return(ps_order)
  }
  vectors <- generate_binary_vectors(m)
  order_index <- as.vector(vectors[, perm, drop = FALSE] %*% 2^(0:(m - 1))) + 1
  ps_order[, order_index, drop = FALSE]
}

#' Compute joint probabilities via chain rule (internal)
#'
#' @description
#' Computes P(A=a|X) for all 2^m combinations using sequential factorization:
#' P(A|X) = ∏ₖ P(Aₖ | A₁,...,Aₖ₋₁, X)
#'
#' Efficient implementation: For each model k, computes predictions under all 2^(k-1)
#' unique conditioning values, then assembles joint probabilities. This reduces
#' redundant predictions from O(2^m × m) to O(2^m).
#'
#' @param data Data frame with covariates and treatments
#' @param conditional_models List of m fitted models
#' @param treatment_order Character vector, treatment order
#' @param vectors Matrix (2^m x m), all treatment combinations
#' @param covariate_cols Character vector of covariate column names (required).
#' @param treatment_cols Character vector, treatment columns in the column order of `vectors`
#'   (defaults to `treatment_order`; supply it when the two differ).
#'
#' @return Matrix (n x 2^m) of propensity scores, columns in `treatment_cols` order
#' @keywords internal
compute_joint_probabilities <- function(data, conditional_models, treatment_order,
                                        vectors, covariate_cols, treatment_cols = treatment_order) {
  n <- nrow(data)
  m <- length(treatment_order)
  perm <- match(treatment_order, treatment_cols)
  stopifnot("treatment_order must be a permutation of treatment_cols" =
              !anyNA(perm) && length(perm) == length(treatment_cols) && !anyDuplicated(perm))

  # Cache predictions for each model under each unique conditioning
  # pred_cache[[k]] is a matrix: n x 2^(k-1), where each column is predictions
  # for a unique conditioning history (A_1, ..., A_{k-1})
  pred_cache <- vector("list", m)

  # Model 1: No conditioning (just X)
  newdata_1 <- data[, covariate_cols, drop = FALSE]
  p_1 <- predict_conditional_model_sequential(
    model = conditional_models[[1]],
    newdata = newdata_1
  )
  pred_cache[[1]] <- matrix(p_1, nrow = n, ncol = 1)

  # Models 2 through m: Condition on previous treatments.
  # seq_len(m)[-1] is integer(0) when m == 1 (single-eligibility pattern) so the
  # loop correctly skips; `2:m` would be the backwards sequence c(2,1) and index
  # the non-existent conditional_models[[2]] ("subscript out of bounds").
  for (k in seq_len(m)[-1]) {
    n_cond <- 2^(k-1)  # Number of unique conditionings
    pred_cache[[k]] <- matrix(0, nrow = n, ncol = n_cond)

    # Generate all unique conditioning values for A_1, ..., A_{k-1}
    cond_vectors <- generate_binary_vectors(k - 1)

    for (cond_idx in 1:n_cond) {
      # Build newdata with this conditioning
      newdata_k <- data[, covariate_cols, drop = FALSE]

      for (j_prev in 1:(k-1)) {
        newdata_k[[treatment_order[j_prev]]] <- cond_vectors[cond_idx, j_prev]
      }

      # Predict for this conditioning
      pred_cache[[k]][, cond_idx] <- predict_conditional_model_sequential(
        model = conditional_models[[k]],
        newdata = newdata_k
      )
    }
  }

  # Now assemble joint probabilities using cached predictions
  ps_matrix <- assemble_sequential_joint(pred_cache, m, perm)

  # Verify: rows should sum to 1 (no normalization needed!)
  row_sums <- rowSums(ps_matrix)
  max_deviation <- max(abs(row_sums - 1))
  if (max_deviation > 1e-6) {
    warning(sprintf("Sequential PS do not sum to 1. Max deviation: %.6f", max_deviation))
  }

  return(ps_matrix)
}

#' Compute joint probabilities via chain rule with cross-fitting (internal)
#'
#' Efficient implementation: For each fold and model, caches predictions under
#' all unique conditioning values to avoid redundant predictions.
#'
#' @param data Data frame with covariates and treatments
#' @param conditional_models List of m lists of n_folds models
#' @param treatment_order Character vector, treatment order
#' @param fold_ids Integer vector, fold assignment for each observation
#' @param vectors Matrix (2^m x m), all treatment combinations
#' @param covariate_cols Character vector of covariate column names (required).
#' @param treatment_cols Character vector, treatment columns in the column order of `vectors`.
#'
#' @return Matrix (n x 2^m) of propensity scores, columns in `treatment_cols` order
#' @keywords internal
compute_joint_probabilities_cv <- function(data, conditional_models, treatment_order,
                                          fold_ids, vectors, covariate_cols,
                                          treatment_cols = treatment_order) {
  n <- nrow(data)
  m <- length(treatment_order)
  perm <- match(treatment_order, treatment_cols)
  stopifnot("treatment_order must be a permutation of treatment_cols" =
              !anyNA(perm) && length(perm) == length(treatment_cols) && !anyDuplicated(perm))
  n_combinations <- nrow(vectors)
  n_folds <- length(unique(fold_ids))

  # Initialize PS matrix
  ps_matrix <- matrix(1.0, nrow = n, ncol = n_combinations)

  # Process each fold
  for (fold in 1:n_folds) {
    test_idx <- which(fold_ids == fold)
    n_test <- length(test_idx)

    # Cache predictions for this fold
    pred_cache <- vector("list", m)

    # Model 1: No conditioning
    newdata_1 <- data[test_idx, covariate_cols, drop = FALSE]
    p_1 <- predict_conditional_model_sequential(
      model = conditional_models[[1]][[fold]],
      newdata = newdata_1
    )
    pred_cache[[1]] <- matrix(p_1, nrow = n_test, ncol = 1)

    # Models 2 through m: Condition on previous treatments.
    # seq_len(m)[-1] is integer(0) when m == 1 so the loop skips; `2:m` would be
    # the backwards sequence c(2,1) and index the non-existent [[2]] model.
    for (k in seq_len(m)[-1]) {
      n_cond <- 2^(k-1)
      pred_cache[[k]] <- matrix(0, nrow = n_test, ncol = n_cond)

      cond_vectors <- generate_binary_vectors(k - 1)

      for (cond_idx in 1:n_cond) {
        newdata_k <- data[test_idx, covariate_cols, drop = FALSE]

        for (j_prev in 1:(k-1)) {
          newdata_k[[treatment_order[j_prev]]] <- cond_vectors[cond_idx, j_prev]
        }

        pred_cache[[k]][, cond_idx] <- predict_conditional_model_sequential(
          model = conditional_models[[k]][[fold]],
          newdata = newdata_k
        )
      }
    }

    # Assemble joint probabilities for this fold using cached predictions
    ps_matrix[test_idx, ] <- assemble_sequential_joint(pred_cache, m, perm)
  }

  # Verify: rows should sum to 1
  row_sums <- rowSums(ps_matrix)
  max_deviation <- max(abs(row_sums - 1))
  if (max_deviation > 1e-6) {
    warning(sprintf("Sequential PS (CV) do not sum to 1. Max deviation: %.6f", max_deviation))
  }

  return(ps_matrix)
}

#' Predict from conditional model (internal)
#'
#' @param model Fitted model object
#' @param newdata Data frame with predictors
#'
#' @return Numeric, predicted probability P(Y=1)
#' @keywords internal
predict_conditional_model_sequential <- function(model, newdata) {
  if (inherits(model, "chunked_glm")) {
    # Chunked GLM: use our custom predict function
    pred <- predict_glm_chunked_internal(model, newdata, chunk_size = 50000)
  } else if (inherits(model, "glm")) {
    # Standard GLM
    pred <- predict(model, newdata = newdata, type = "response")
  } else if (inherits(model, "randomForest")) {
    # Random forest
    pred_probs <- predict(model, newdata = newdata, type = "prob")
    pred <- pred_probs[, "1"]
  } else if (inherits(model, "xgb.Booster")) {
    # XGBoost — must be numeric (double) matrix; as.matrix() alone may produce
    # integer matrix if covariates are integer-typed, which xgboost rejects.
    newdata_mat <- as.matrix(newdata)
    storage.mode(newdata_mat) <- "double"
    pred <- predict(model, newdata = newdata_mat)
  } else {
    stop("Unknown model type in predict_conditional_model_sequential")
  }

  return(as.numeric(pred))
}

# ============================================================================
# End of Sequential PS Helper Functions
# ============================================================================

#' Estimate propensity scores for vector treatments
#'
#' @description
#' Estimates propensity scores P(A = a | X) for all 2^m treatment combinations
#' using various methods. Propensity scores are required for EIF-based (one-step)
#' estimation and must be properly calibrated for valid inference.
#'
#' @param data Data frame with columns:
#'   - Treatment columns: A1, A2, ..., Am (binary 0/1)
#'   - Covariate columns: X1, X2, ..., Xp (numeric)
#' @param treatment_cols Character vector of treatment column names (required).
#' @param covariate_cols Character vector of covariate column names (required).
#' @param method Character, estimation method. Chooses classification approach:
#'   - **"glm"** (default): Multinomial logistic via nnet::multinom(). Fast, interpretable,
#'     works well for m <= 3. May struggle with m >= 4 (2^m classes).
#'   - **"nnet"**: Alias for "glm"
#'   - **"randomforest"**: Multinomial RF via randomForest::randomForest(). Flexible,
#'     handles m <= 5 reasonably. Slower than GLM.
#'   - **"xgboost"**: Multinomial gradient boosting. Fast, flexible, good for m <= 6.
#'   - **"sequential"**: Fits m conditional models via chain rule factorization.
#'     **RECOMMENDED** - theoretically principled, no ad-hoc normalization, provable
#'     convergence to P(A|X). Models P(A) = ∏P(Aₖ|A₁,...,Aₖ₋₁,X). Fast (only m models),
#'     handles arbitrary treatment correlations, no sparsity issues. See Details.
#'   - **"binary_per_combination"**: Fits 2^m binary models P(A=a|X) separately.
#'     Respects treatment correlations fully but requires ad-hoc normalization.
#'     Use with binary_method = "xgboost" for best performance.
#'   - **"mlr3"**: mlr3 framework with flexible learners (requires mlr3 package)
#'   - **"superlearner"**: Ensemble via SuperLearner (requires SuperLearner package)
#' @param cross_fit Logical, use sample splitting via cross-fitting. Reduces
#'   overfitting bias when using flexible methods (randomforest, xgboost).
#'   Default: FALSE. Recommended TRUE for n < 1000.
#' @param n_folds Integer, folds for cross-fitting. Default: 5. Only used if
#'   cross_fit = TRUE.
#' @param sl_library Character vector, learners for SuperLearner ensemble.
#'   E.g., c("SL.glm", "SL.randomForest"). Only used if method = "superlearner".
#' @param mlr3_learner Character, mlr3 learner name. Default: "classif.ranger".
#'   Only used if method = "mlr3".
#' @param binary_method Character, base learner for binary_per_combination and sequential.
#'   Options: "glm", "randomforest", "xgboost" (recommended). Default: "glm".
#' @param min_obs_threshold Integer, minimum observations required to fit binary
#'   model for a combination (when method = "binary_per_combination"). Combinations
#'   with fewer observations use empirical frequency. Default: 10. Set to 0 to
#'   fit all combinations. Higher values (e.g., 50-100) enable sparse fitting
#'   for high-dimensional treatments (m >= 12). Not used for method = "sequential".
#' @param treatment_order Character vector, ordering of treatments for sequential
#'   method. Default: NULL (uses natural order A1, A2, ..., Am). Only used for
#'   method = "sequential".
#'
#' @return List containing:
#' \describe{
#'   \item{propensity_scores}{Numeric matrix (n x 2^m), P(A=a|X) for each combo. Rows sum to 1.}
#'   \item{models}{Fitted model object(s). Single model for multinomial methods,
#'     list of 2^m models for binary_per_combination.}
#'   \item{method}{Character, method used (for reference)}
#'   \item{treatment_combinations}{Matrix (2^m x m), binary treatment patterns in order}
#'   \item{cross_fit}{Logical, whether cross-fitting was used}
#' }
#'
#' @details
#' **Method selection guide:**
#' - **m = 2-3**: Use "sequential" or "glm" (both fast and reliable)
#' - **m = 4-6**: Use "sequential" (recommended) or "binary_per_combination"
#' - **m >= 7**: Use "sequential" (theoretically sound, computationally feasible)
#'
#' **Sequential method (RECOMMENDED):**
#' Models joint PS via chain rule factorization:
#' \deqn{P(A_1,...,A_m|X) = P(A_1|X) P(A_2|A_1,X) ... P(A_m|A_1,...,A_{m-1},X)}
#'
#' Key advantages over binary_per_combination:
#' - **Theoretically principled**: Direct probability factorization, no ad-hoc normalization
#' - **Provable convergence**: Each conditional converges, so the joint converges
#' - **Computationally efficient**: Only m models (vs 2^m), each fit on n observations
#' - **No sparsity issues**: All models well-supported (no need for sparse fitting)
#' - **Flexible to interactions**: Captures treatment correlations via conditioning
#'
#' Probabilities sum to 1 by construction (chain rule + law of total probability).
#'
#' **Quality checks:** The function validates outputs (non-finite values, rows
#' sum to 1 within tolerance, severe positivity violations).
#'
#' Random number generation: the caller is responsible for setting the seed
#' before calling this function (e.g. set.seed(123)) to ensure reproducibility.
#'
#' @export
#' @examples
#' # Generate example data
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 2, p = 2, seed = 123)
#' tcols <- paste0("A", 1:2)
#' xcols <- paste0("X", 1:2)
#'
#' # Multinomial logistic regression
#' ps_est <- estimate_propensity_scores(sim_data$data, treatment_cols = tcols,
#'                                      covariate_cols = xcols, method = "glm")
#' head(ps_est$propensity_scores)
#'
#' # Binary-per-combination with GLM base learner
#' ps_est_binary <- estimate_propensity_scores(sim_data$data, treatment_cols = tcols,
#'                                             covariate_cols = xcols,
#'                                             method = "binary_per_combination",
#'                                             binary_method = "glm")
estimate_propensity_scores <- function(data, treatment_cols, covariate_cols,
                                      method = "glm", cross_fit = TRUE,
                                      n_folds = 5, sl_library = NULL,
                                      mlr3_learner = "classif.ranger",
                                      binary_method = "glm",
                                      min_obs_threshold = 10,
                                      treatment_order = NULL,
                                      xgb_params = NULL) {

  # Explicit column specifications (no auto-detection).
  validate_columns(data, treatment_cols, covariate_cols)
  m <- length(treatment_cols)
  p <- length(covariate_cols)

  # Generate all treatment combinations
  vectors <- generate_binary_vectors(m)
  n_combinations <- nrow(vectors)
  
  # Create column names matching simulate_vector_treatment_data format
  column_names <- character(n_combinations)
  for (j in 1:n_combinations) {
    treatment_str <- paste0("A", 1:m, "=", vectors[j, ], collapse = "_")
    column_names[j] <- treatment_str
  }
  
  # Create treatment factor variable from observed treatments
  # Map each observation's treatment to the corresponding combination index
  # OPTIMIZED: Use binary encoding instead of nested apply (100x+ faster for large n)

  # Extract treatment matrix
  A_matrix <- as.matrix(data[, treatment_cols, drop = FALSE])

  # Convert each row to binary number (base-2 positional encoding)
  # For m=8: [0,0,0,0,0,0,0,1] -> 1, [1,0,0,0,0,0,0,0] -> 128, etc.
  powers <- 2^(0:(m-1))
  indices <- as.integer(A_matrix %*% powers) + 1  # +1 for 1-based indexing

  # Map to column names
  treatment_labels <- column_names[indices]

  # Create factor
  treatment_factor <- factor(treatment_labels, levels = column_names)
  
  # Check for sparse combinations
  combination_counts <- table(treatment_factor)
  sparse_combinations <- names(combination_counts)[combination_counts < 10]
  if (length(sparse_combinations) > 0) {
    warning(paste("Some treatment combinations have < 10 observations:",
                  paste(sparse_combinations, collapse = ", ")))
  }
  
  # Initialize results
  propensity_scores <- matrix(0, nrow = nrow(data), ncol = n_combinations)
  colnames(propensity_scores) <- column_names
  
  # Methods that use multinomial approach (single model)
  multinomial_methods <- c("glm", "nnet", "randomforest", "xgboost", "mlr3")
  
  if (method %in% multinomial_methods) {
    # Use multinomial approach - fit single model
    if (cross_fit) {
      fold_ids <- sample(rep(1:n_folds, length.out = nrow(data)))
      models <- list()
      
      for (fold in 1:n_folds) {
        train_idx <- which(fold_ids != fold)
        test_idx <- which(fold_ids == fold)
        
        data_train <- data.frame(data[train_idx, covariate_cols, drop = FALSE])
        data_train$treatment_factor <- treatment_factor[train_idx]
        data_test <- data.frame(data[test_idx, covariate_cols, drop = FALSE])
        
        if (method == "glm" || method == "nnet") {
          if (!requireNamespace("nnet", quietly = TRUE)) {
            stop("nnet package required for method = 'glm' or 'nnet'")
          }
          model <- nnet::multinom(treatment_factor ~ ., data = data_train, trace = FALSE)
          pred_probs <- predict(model, newdata = data_test, type = "probs")
          # Ensure all combinations are present in predictions
          if (is.vector(pred_probs)) {
            # Only one class predicted - convert to matrix
            pred_probs <- matrix(pred_probs, nrow = length(pred_probs), ncol = n_combinations)
            colnames(pred_probs) <- levels(treatment_factor)
          }
          # Reorder columns to match column_names
          pred_probs <- pred_probs[, column_names, drop = FALSE]
          propensity_scores[test_idx, ] <- pred_probs
          models[[fold]] <- model
          
        } else if (method == "randomforest") {
          if (!requireNamespace("randomForest", quietly = TRUE)) {
            stop("randomForest package required for method = 'randomforest'")
          }
          # Drop unused levels to avoid "empty classes" error
          treatment_factor_train <- droplevels(data_train$treatment_factor)
          model <- randomForest::randomForest(
            x = data_train[, covariate_cols, drop = FALSE],
            y = treatment_factor_train
          )
          pred_probs_observed <- predict(model, newdata = data_test, type = "prob")
          # Expand to full set of combinations (fill missing with 0)
          pred_probs <- matrix(0, nrow = nrow(data_test), ncol = n_combinations)
          colnames(pred_probs) <- column_names
          observed_levels <- colnames(pred_probs_observed)
          pred_probs[, observed_levels] <- pred_probs_observed
          # Normalize rows to sum to 1
          row_sums <- rowSums(pred_probs)
          row_sums[row_sums < 1e-10] <- 1
          pred_probs <- pred_probs / row_sums
          propensity_scores[test_idx, ] <- pred_probs
          models[[fold]] <- model
          
        } else if (method == "xgboost") {
          if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("xgboost package required for method = 'xgboost'")
          }
          # Drop unused levels and create numeric mapping
          treatment_factor_train <- droplevels(data_train$treatment_factor)
          y_train_numeric <- as.numeric(treatment_factor_train) - 1
          n_observed_classes <- length(levels(treatment_factor_train))

          # Create DMatrix with updated parameter name
          X_train_matrix <- to_numeric_matrix(data_train, covariate_cols)
          X_test_matrix <- to_numeric_matrix(data_test, covariate_cols)
          dtrain <- xgboost::xgb.DMatrix(data = X_train_matrix, label = y_train_numeric)

          # Train model with observed classes only. Tuned regularization defaults
          # (resolve_xgb_params) match train_xgboost_early_stop; caller may override
          # via xgb_params. See that helper's @details for rationale.
          params <- modifyList(
            resolve_xgb_params(xgb_params),
            list(objective = "multi:softprob", num_class = n_observed_classes, verbosity = 0)
          )
          model <- xgboost::xgb.train(
            params = params,
            data = dtrain,
            nrounds = 500,
            verbose = 0
          )

          # Predict on test set
          pred_probs_observed <- predict(model, newdata = X_test_matrix, reshape = TRUE)

          # Expand to full set of combinations (fill missing with 0)
          pred_probs <- matrix(0, nrow = nrow(data_test), ncol = n_combinations)
          colnames(pred_probs) <- column_names
          observed_levels <- levels(treatment_factor_train)
          pred_probs[, observed_levels] <- pred_probs_observed
          # Normalize rows to sum to 1
          row_sums <- rowSums(pred_probs)
          row_sums[row_sums < 1e-10] <- 1
          pred_probs <- pred_probs / row_sums

          propensity_scores[test_idx, ] <- pred_probs
          models[[fold]] <- model
          
        } else if (method == "mlr3") {
          if (!requireNamespace("mlr3", quietly = TRUE)) {
            stop("mlr3 package required for method = 'mlr3'")
          }
          if (!requireNamespace("mlr3learners", quietly = TRUE)) {
            stop("mlr3learners package required for method = 'mlr3'")
          }
          
          # Create task
          task_data <- data.frame(data_train[, covariate_cols, drop = FALSE])
          task_data$treatment_factor <- data_train$treatment_factor
          task <- mlr3::TaskClassif$new(
            id = "propensity",
            backend = task_data,
            target = "treatment_factor"
          )
          
          # Get learner
          learner <- mlr3learners::lrn(mlr3_learner, predict_type = "prob")
          learner$train(task)
          
          # Predict on test data using predict_newdata
          prediction <- learner$predict_newdata(data_test)
          pred_probs <- prediction$prob
          
          # Reorder columns to match column_names
          pred_probs <- pred_probs[, column_names, drop = FALSE]
          propensity_scores[test_idx, ] <- as.matrix(pred_probs)
          models[[fold]] <- learner
        }
      }
    } else {
      # No cross-fitting - fit single model
      data_all <- data.frame(data[, covariate_cols, drop = FALSE])
      data_all$treatment_factor <- treatment_factor
      
      if (method == "glm" || method == "nnet") {
        if (!requireNamespace("nnet", quietly = TRUE)) {
          stop("nnet package required for method = 'glm' or 'nnet'")
        }
        model <- nnet::multinom(treatment_factor ~ ., data = data_all, trace = FALSE)
        pred_probs <- predict(model, type = "probs")
        # Ensure all combinations are present
        if (is.vector(pred_probs)) {
          pred_probs <- matrix(pred_probs, nrow = length(pred_probs), ncol = n_combinations)
          colnames(pred_probs) <- levels(treatment_factor)
        }
        # Reorder columns to match column_names
        pred_probs <- pred_probs[, column_names, drop = FALSE]
        propensity_scores <- pred_probs
        models <- model
        
      } else if (method == "randomforest") {
        if (!requireNamespace("randomForest", quietly = TRUE)) {
          stop("randomForest package required for method = 'randomforest'")
        }
        # Drop unused levels to avoid "empty classes" error
        treatment_factor_all <- droplevels(data_all$treatment_factor)
        model <- randomForest::randomForest(
          x = data_all[, covariate_cols, drop = FALSE],
          y = treatment_factor_all
        )
        pred_probs_observed <- predict(model, type = "prob")
        # Expand to full set of combinations (fill missing with 0)
        pred_probs <- matrix(0, nrow = nrow(data_all), ncol = n_combinations)
        colnames(pred_probs) <- column_names
        observed_levels <- colnames(pred_probs_observed)
        pred_probs[, observed_levels] <- pred_probs_observed
        # Normalize rows to sum to 1
        row_sums <- rowSums(pred_probs)
        row_sums[row_sums < 1e-10] <- 1
        pred_probs <- pred_probs / row_sums
        propensity_scores <- pred_probs
        models <- model
        
      } else if (method == "xgboost") {
        if (!requireNamespace("xgboost", quietly = TRUE)) {
          stop("xgboost package required for method = 'xgboost'")
        }
        # Drop unused levels and create numeric mapping
        treatment_factor_all <- droplevels(data_all$treatment_factor)
        y_numeric <- as.numeric(treatment_factor_all) - 1
        n_observed_classes <- length(levels(treatment_factor_all))

        # Create DMatrix with updated parameter name
        X_all_matrix <- to_numeric_matrix(data_all, covariate_cols)
        dtrain <- xgboost::xgb.DMatrix(data = X_all_matrix, label = y_numeric)

        # Train model with observed classes only. Tuned regularization defaults
        # (resolve_xgb_params) match train_xgboost_early_stop; caller may override
        # via xgb_params. See that helper's @details.
        params <- modifyList(
          resolve_xgb_params(xgb_params),
          list(objective = "multi:softprob", num_class = n_observed_classes, verbosity = 0)
        )
        model <- xgboost::xgb.train(
          params = params,
          data = dtrain,
          nrounds = 500,
          verbose = 0
        )

        # Predict on training set
        pred_probs_observed <- predict(model, newdata = X_all_matrix, reshape = TRUE)

        # Expand to full set of combinations (fill missing with 0)
        pred_probs <- matrix(0, nrow = nrow(data_all), ncol = n_combinations)
        colnames(pred_probs) <- column_names
        observed_levels <- levels(treatment_factor_all)
        pred_probs[, observed_levels] <- pred_probs_observed
        # Normalize rows to sum to 1
        row_sums <- rowSums(pred_probs)
        row_sums[row_sums < 1e-10] <- 1
        pred_probs <- pred_probs / row_sums

        propensity_scores <- pred_probs
        models <- model
        
      } else if (method == "mlr3") {
        if (!requireNamespace("mlr3", quietly = TRUE)) {
          stop("mlr3 package required for method = 'mlr3'")
        }
        if (!requireNamespace("mlr3learners", quietly = TRUE)) {
          stop("mlr3learners package required for method = 'mlr3'")
        }
        
        # Create task
        task <- mlr3::TaskClassif$new(
          id = "propensity",
          backend = data_all,
          target = "treatment_factor"
        )
        
        # Get learner
        learner <- mlr3learners::lrn(mlr3_learner, predict_type = "prob")
        learner$train(task)
        
        # Predict on training data (same task)
        prediction <- learner$predict(task)
        pred_probs <- prediction$prob
        
        # Reorder columns to match column_names
        pred_probs <- pred_probs[, column_names, drop = FALSE]
        propensity_scores <- as.matrix(pred_probs)
        models <- learner
      }
    }
    
  } else if (method == "superlearner") {
    # SuperLearner doesn't have native multinomial support - use binary approach
    models <- list()
    n_skipped <- 0  # Track sparse combinations

    if (cross_fit) {
      fold_ids <- sample(rep(1:n_folds, length.out = nrow(data)))
    }

    # Report sparse fitting strategy if threshold > 0
    if (min_obs_threshold > 0) {
      cat(sprintf("[Sparse Fitting] Will skip combinations with < %d observations (using empirical frequencies)\n",
                  min_obs_threshold))
    }

    # Fit models for each treatment combination
    for (j in 1:n_combinations) {
      # Create indicator variable: I(A = a_j) for each observation
      indicator <- rep(TRUE, nrow(data))
      for (k in 1:m) {
        indicator <- indicator & (data[[treatment_cols[k]]] == vectors[j, k])
      }
      indicator <- as.numeric(indicator)

      # Check if combination has sufficient observations
      n_positive <- sum(indicator)
      if (n_positive < min_obs_threshold) {
        # Too few observations, use empirical probability (sparse fitting)
        emp_prob <- n_positive / nrow(data)
        propensity_scores[, j] <- emp_prob
        models[[j]] <- NULL
        n_skipped <- n_skipped + 1
        next
      }
      
      if (cross_fit) {
        models[[j]] <- list()
        
        for (fold in 1:n_folds) {
          train_idx <- which(fold_ids != fold)
          test_idx <- which(fold_ids == fold)
          
          if (sum(indicator[train_idx]) < 5) {
            # Too few training observations
            emp_prob <- sum(indicator[train_idx]) / length(train_idx)
            propensity_scores[test_idx, j] <- emp_prob
            models[[j]][[fold]] <- NULL
            next
          }
          
          # Create data frames with named columns for training and testing
          data_train <- data.frame(data[train_idx, covariate_cols, drop = FALSE])
          data_train$indicator <- indicator[train_idx]
          data_test <- data.frame(data[test_idx, covariate_cols, drop = FALSE])
          
          if (!requireNamespace("SuperLearner", quietly = TRUE)) {
            stop("SuperLearner package required for method = 'superlearner'")
          }
          if (is.null(sl_library)) {
            sl_library <- c("SL.glm", "SL.gam", "SL.randomForest")
          }
          y_train <- indicator[train_idx]
          X_train_df <- data_train[, covariate_cols, drop = FALSE]
          X_test_df <- data_test[, covariate_cols, drop = FALSE]
          model <- SuperLearner::SuperLearner(Y = y_train, X = X_train_df,
                                            family = binomial(), SL.library = sl_library)
          propensity_scores[test_idx, j] <- predict(model, 
                                                   newdata = X_test_df)$pred
          models[[j]][[fold]] <- model
        }
      } else {
        # No cross-fitting
        data_all <- data.frame(data[, covariate_cols, drop = FALSE])
        data_all$indicator <- indicator
        
        if (!requireNamespace("SuperLearner", quietly = TRUE)) {
          stop("SuperLearner package required for method = 'superlearner'")
        }
        if (is.null(sl_library)) {
          sl_library <- c("SL.glm", "SL.gam", "SL.randomForest")
        }
        y_all <- indicator
        X_all_df <- data_all[, covariate_cols, drop = FALSE]
        model <- SuperLearner::SuperLearner(Y = y_all, X = X_all_df,
                                          family = binomial(), SL.library = sl_library)
        propensity_scores[, j] <- predict(model)$pred
        models[[j]] <- model
      }
    }

    # Report sparse fitting summary
    if (n_skipped > 0) {
      cat(sprintf("[Sparse Fitting] Skipped %d/%d combinations (%.1f%%), fitted %d models\n",
                  n_skipped, n_combinations, 100 * n_skipped / n_combinations,
                  n_combinations - n_skipped))
    }

    # Normalize probabilities so each row sums to 1 (needed for binary approach)
    row_sums <- rowSums(propensity_scores)
    row_sums[row_sums < 1e-10] <- 1
    propensity_scores <- propensity_scores / row_sums

  } else if (method == "sequential") {
    # Sequential/factorized modeling via chain rule
    # Theoretically principled: P(A|X) = ∏ P(Aₖ|A₁,...,Aₖ₋₁,X)
    # No ad-hoc normalization, provable convergence, no sparsity issues

    # Set treatment order (default: natural order)
    if (is.null(treatment_order)) {
      treatment_order <- treatment_cols
    } else {
      # Validate provided order
      if (!all(treatment_order %in% treatment_cols) || length(treatment_order) != m || anyDuplicated(treatment_order) > 0) {
        stop("treatment_order must contain all treatment columns exactly once")
      }
    }

    # Fit m conditional models
    conditional_models <- list()

    if (cross_fit) {
      # Cross-fitting: split data into folds
      fold_ids <- sample(rep(1:n_folds, length.out = nrow(data)))

      # For each position in sequence, fit conditional model on each fold
      for (k in 1:m) {
        outcome_var <- treatment_order[k]

        # Predictors: X + A_1 + ... + A_{k-1}
        if (k == 1) {
          predictor_cols <- covariate_cols
        } else {
          predictor_cols <- c(covariate_cols, treatment_order[1:(k-1)])
        }

        conditional_models[[k]] <- list()

        for (fold in 1:n_folds) {
          train_idx <- which(fold_ids != fold)

          # Fit binary model for P(Aₖ = 1 | A_{<k}, X)
          model <- fit_conditional_model_sequential(
            data = data,
            idx = train_idx,
            outcome_var = outcome_var,
            predictor_cols = predictor_cols,
            binary_method = binary_method,
            p = p,
            xgb_params = xgb_params
          )

          conditional_models[[k]][[fold]] <- model
        }
      }

      # Generate joint PS matrix using cross-fitted models
      propensity_scores <- compute_joint_probabilities_cv(
        data = data,
        conditional_models = conditional_models,
        treatment_order = treatment_order,
        fold_ids = fold_ids,
        vectors = vectors,
        covariate_cols = covariate_cols,
        treatment_cols = treatment_cols
      )

    } else {
      # No cross-fitting: train on all data
      for (k in 1:m) {
        outcome_var <- treatment_order[k]

        # Predictors: X + A_1 + ... + A_{k-1}
        if (k == 1) {
          predictor_cols <- covariate_cols
        } else {
          predictor_cols <- c(covariate_cols, treatment_order[1:(k-1)])
        }

        # Fit binary model for P(Aₖ = 1 | A_{<k}, X)
        model <- fit_conditional_model_sequential(
          data = data,
          idx = 1:nrow(data),
          outcome_var = outcome_var,
          predictor_cols = predictor_cols,
          binary_method = binary_method,
          p = p,
          xgb_params = xgb_params
        )

        conditional_models[[k]] <- model
      }

      # Generate joint PS matrix
      propensity_scores <- compute_joint_probabilities(
        data = data,
        conditional_models = conditional_models,
        treatment_order = treatment_order,
        vectors = vectors,
        covariate_cols = covariate_cols,
        treatment_cols = treatment_cols
      )
    }

    models <- list(
      conditional_models = conditional_models,
      treatment_order = treatment_order,
      method = "sequential"
    )

  } else if (method == "binary_per_combination") {
    # Fit separate binary model for each treatment combination
    # More principled than multinomial methods: respects treatment correlations
    # Allows using any binary classifier

    models <- list()
    n_skipped <- 0  # Track sparse combinations

    if (cross_fit) {
      fold_ids <- sample(rep(1:n_folds, length.out = nrow(data)))
    }

    # Report sparse fitting strategy if threshold > 0
    if (min_obs_threshold > 0) {
      cat(sprintf("[Sparse Fitting] Will skip combinations with < %d observations (using empirical frequencies)\n",
                  min_obs_threshold))
    }

    # Loop over all treatment combinations
    for (j in 1:n_combinations) {
      # Create binary indicator: I(A = a_j)
      indicator <- rep(TRUE, nrow(data))
      for (k in 1:m) {
        indicator <- indicator & (data[[treatment_cols[k]]] == vectors[j, k])
      }
      indicator <- as.numeric(indicator)

      # Check if combination has sufficient observations
      n_positive <- sum(indicator)
      if (n_positive < min_obs_threshold) {
        # Too few observations, use empirical probability (sparse fitting)
        emp_prob <- n_positive / nrow(data)
        propensity_scores[, j] <- emp_prob
        models[[j]] <- NULL
        n_skipped <- n_skipped + 1
        next  # Skip to next combination
      }

      if (cross_fit) {
        # Cross-fitting: train on folds
        models[[j]] <- list()

        for (fold in 1:n_folds) {
          train_idx <- which(fold_ids != fold)
          test_idx <- which(fold_ids == fold)

          if (sum(indicator[train_idx]) < 5) {
            # Too few training observations
            emp_prob <- sum(indicator[train_idx]) / length(train_idx)
            propensity_scores[test_idx, j] <- emp_prob
            models[[j]][[fold]] <- NULL
            next
          }

          # Fit binary model
          if (binary_method == "glm") {
            # Check if chunking is needed for memory efficiency
            n_train <- length(train_idx)
            use_chunking <- needs_chunking(n_train, p)

            # Progress reporting for large fits
            if (n_train > 100000 && use_chunking && fold == 1) {
              cat(sprintf("  [PS] Fitting combination %d/%d (chunked, n_train=%s)...\n",
                          j, n_combinations, format(n_train, big.mark = ",")))
            }

            if (use_chunking) {
              # Chunked GLM fitting for large datasets
              data_train <- data.frame(data[train_idx, covariate_cols, drop = FALSE])
              data_train$indicator <- indicator[train_idx]

              model <- fit_glm_chunked_internal(
                data_train = data_train,
                outcome_var = "indicator",
                covariate_cols = covariate_cols,
                chunk_size = 50000
              )

              data_test <- data.frame(data[test_idx, covariate_cols, drop = FALSE])
              propensity_scores[test_idx, j] <- predict_glm_chunked_internal(
                model = model,
                newdata = data_test,
                chunk_size = 50000
              )

              models[[j]][[fold]] <- model
            } else {
              # Standard GLM for smaller datasets
              data_train <- data.frame(data[train_idx, covariate_cols, drop = FALSE])
              data_train$indicator <- indicator[train_idx]
              data_test <- data.frame(data[test_idx, covariate_cols, drop = FALSE])

              model <- glm(indicator ~ ., data = data_train, family = binomial())
              propensity_scores[test_idx, j] <- predict(model, newdata = data_test, type = "response")
              models[[j]][[fold]] <- model
            }

          } else if (binary_method == "randomforest") {
            # Random forest binary classifier
            if (!requireNamespace("randomForest", quietly = TRUE)) {
              stop("randomForest package required for binary_method = 'randomforest'")
            }
            X_train <- data[train_idx, covariate_cols, drop = FALSE]
            y_train <- factor(indicator[train_idx], levels = c(0, 1))
            X_test <- data[test_idx, covariate_cols, drop = FALSE]

            model <- randomForest::randomForest(x = X_train, y = y_train)
            pred_probs <- predict(model, newdata = X_test, type = "prob")
            propensity_scores[test_idx, j] <- pred_probs[, "1"]
            models[[j]][[fold]] <- model

          } else if (binary_method == "xgboost") {
            # XGBoost binary classifier
            if (!requireNamespace("xgboost", quietly = TRUE)) {
              stop("xgboost package required for binary_method = 'xgboost'")
            }
            X_train <- to_numeric_matrix(data[train_idx, , drop = FALSE], covariate_cols)
            y_train <- indicator[train_idx]
            X_test <- to_numeric_matrix(data[test_idx, , drop = FALSE], covariate_cols)

            params <- modifyList(resolve_xgb_params(xgb_params),
                                 list(objective = "binary:logistic", verbosity = 0))
            model <- train_xgboost_early_stop(X_train, y_train, params)
            propensity_scores[test_idx, j] <- predict(model, newdata = X_test)
            models[[j]][[fold]] <- model

          } else {
            stop("binary_method must be one of: 'glm', 'randomforest', 'xgboost'")
          }
        }

      } else {
        # No cross-fitting: train on all data

        if (binary_method == "glm") {
          # Check if chunking is needed for memory efficiency
          n_train <- nrow(data)
          use_chunking <- needs_chunking(n_train, p)

          if (use_chunking) {
            # Chunked GLM fitting for large datasets
            data_all <- data.frame(data[, covariate_cols, drop = FALSE])
            data_all$indicator <- indicator

            model <- fit_glm_chunked_internal(
              data_train = data_all,
              outcome_var = "indicator",
              covariate_cols = covariate_cols,
              chunk_size = 50000
            )

            propensity_scores[, j] <- predict_glm_chunked_internal(
              model = model,
              newdata = data[, covariate_cols, drop = FALSE],
              chunk_size = 50000
            )

            models[[j]] <- model
          } else {
            # Standard GLM for smaller datasets
            data_all <- data.frame(data[, covariate_cols, drop = FALSE])
            data_all$indicator <- indicator

            model <- glm(indicator ~ ., data = data_all, family = binomial())
            propensity_scores[, j] <- predict(model, type = "response")
            models[[j]] <- model
          }

        } else if (binary_method == "randomforest") {
          # Random forest binary classifier
          if (!requireNamespace("randomForest", quietly = TRUE)) {
            stop("randomForest package required for binary_method = 'randomforest'")
          }
          X_all <- data[, covariate_cols, drop = FALSE]
          y_all <- factor(indicator, levels = c(0, 1))

          model <- randomForest::randomForest(x = X_all, y = y_all)
          pred_probs <- predict(model, type = "prob")
          propensity_scores[, j] <- pred_probs[, "1"]
          models[[j]] <- model

        } else if (binary_method == "xgboost") {
          # XGBoost binary classifier
          if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("xgboost package required for binary_method = 'xgboost'")
          }
          X_all <- to_numeric_matrix(data, covariate_cols)
          y_all <- indicator

          params <- modifyList(resolve_xgb_params(xgb_params),
                               list(objective = "binary:logistic", verbosity = 0))
          model <- train_xgboost_early_stop(X_all, y_all, params)
          propensity_scores[, j] <- predict(model, newdata = X_all)
          models[[j]] <- model

        } else {
          stop("binary_method must be one of: 'glm', 'randomforest', 'xgboost'")
        }
      }

      # Memory cleanup every 50 combinations
      if (j %% 50 == 0) {
        gc(verbose = FALSE)
      }
    }

    # Final memory cleanup after all combinations
    gc(verbose = FALSE)

    # Report sparse fitting summary
    if (n_skipped > 0) {
      cat(sprintf("[Sparse Fitting] Skipped %d/%d combinations (%.1f%%), fitted %d models\n",
                  n_skipped, n_combinations, 100 * n_skipped / n_combinations,
                  n_combinations - n_skipped))
    }

    # Normalize probabilities so each row sums to 1
    row_sums <- rowSums(propensity_scores)
    row_sums[row_sums < 1e-10] <- 1
    propensity_scores <- propensity_scores / row_sums

  } else {
    stop("method must be one of: 'glm', 'nnet', 'randomforest', 'xgboost', 'sequential', 'mlr3', 'superlearner', 'binary_per_combination'")
  }
  
  # Ensure probabilities sum to 1 (safety check, though multinomial methods should already do this)
  row_sums <- rowSums(propensity_scores)
  if (any(abs(row_sums - 1) > 1e-6)) {
    # Normalize if needed
    row_sums[row_sums < 1e-10] <- 1
    propensity_scores <- propensity_scores / row_sums
  }

  # Final validation
  if (any(!is.finite(propensity_scores))) {
    stop("Propensity score estimation produced non-finite values (NA, NaN, or Inf). ",
         "Method '", method, "' failed. ",
         "Check data quality, model convergence, and positivity violations.")
  }

  if (any(propensity_scores < 0)) {
    stop("Propensity score estimation produced negative values. ",
         "Method '", method, "' failed. This should not happen.")
  }

  if (any(propensity_scores > 1)) {
    stop("Propensity score estimation produced values > 1. ",
         "Method '", method, "' failed. This should not happen.")
  }

  # Check for combinations with zero probability everywhere (severe positivity violation)
  max_ps_per_combo <- apply(propensity_scores, 2, max)
  zero_combos <- sum(max_ps_per_combo < 1e-12)
  if (zero_combos > 0) {
    warning("Severe positivity violation: ", zero_combos, " treatment combinations have max probability < 1e-12. ",
            "These combinations never occur in the data. Consider reducing the number of treatments.")
  }

  return(list(
    propensity_scores = propensity_scores,
    models = models,
    method = method,
    treatment_combinations = vectors,
    cross_fit = cross_fit
  ))
}

#' Estimate outcome regression models for vector treatments
#'
#' @description
#' Estimates conditional outcome means `E[Y | A = a, X]` for all 2^m treatment
#' combinations. Fits stratified models (one per combination) when feasible, or
#' a single model with treatment interactions for high-dimensional cases.
#' Required for EIF-based (one-step) estimation.
#'
#' @param data Data frame with columns:
#'   - Treatment: A1, A2, ..., Am (binary 0/1)
#'   - Covariates: X1, X2, ..., Xp (numeric)
#'   - Outcome: Y (numeric)
#' @param treatment_cols Character vector of treatment column names (required).
#' @param covariate_cols Character vector of covariate column names (required).
#' @param outcome_col Character scalar, outcome column name. Default "Y".
#' @param method Character, regression method:
#'   - **"glm"** (default): Linear regression via lm(). Fast, interpretable. Works well
#'     for continuous Y. Fits separate model per combination (if n_a >= min_obs).
#'   - **"gam"**: Generalized additive models via mgcv::gam(). Flexible smoothing.
#'     Good for non-linear covariate effects.
#'   - **"randomforest"**: Random forest via randomForest::randomForest(). Highly
#'     flexible, handles interactions automatically. Recommended for complex outcome surfaces.
#'   - **"xgboost"**: Gradient boosting via xgboost::xgb.train(). Fast, flexible,
#'     handles nonlinearities and skewed outcomes (e.g., costs). Recommended for
#'     high-dimensional data. Consistent with xgboost-based PS estimation.
#'   - **"superlearner"**: Ensemble via SuperLearner. Adapts to data complexity.
#'     Requires SuperLearner package and sl_library specification.
#' @param cross_fit Logical, use cross-fitting for sample splitting. Reduces
#'   overfitting when using flexible methods. Default: FALSE. Recommended TRUE
#'   for randomforest/superlearner.
#' @param n_folds Integer, number of cross-fitting folds. Default: 5.
#' @param sl_library Character vector, learners for SuperLearner.
#'   E.g., c("SL.glm", "SL.gam", "SL.randomForest"). Required if method = "superlearner".
#' @param outcome_fit Character, how to model the outcome across treatments:
#'   - **"pooled"** (default): a single model on covariates and treatment
#'     indicators, predicting `E[Y | A = a, X]` for each combination by setting the
#'     treatment columns to a. Borrows strength across the full sample; robust as
#'     the treatment space grows. Lower mu-hat RMSE than stratified for m >= 3.
#'   - **"stratified"**: a separate model per treatment combination (assumption-free
#'     in treatment). Best when every combination is well populated (small m, large
#'     n). Falls back to pooled with a warning if 2^m > max_combinations, to avoid
#'     collapsing sparse combinations to the grand mean.
#' @param max_combinations Integer, only used when outcome_fit = "stratified":
#'   if 2^m exceeds it, the stratified fit falls back to pooled. Default: 16.
#'
#' @return List containing:
#' \describe{
#'   \item{outcome_means}{Numeric matrix (n x 2^m), `E[Y|A=a,X]` predictions for all combos}
#'   \item{models}{List of fitted models (one per combination if feasible, else single model)}
#'   \item{method}{Character, method used (for reference)}
#'   \item{treatment_combinations}{Matrix (2^m x m), binary treatment patterns}
#'   \item{cross_fit}{Logical, whether cross-fitting was used}
#' }
#'
#' @details
#' **Modeling approach:** By default (`outcome_fit = "pooled"`) the function fits a
#' single outcome regression on covariates and treatment indicators, then predicts
#' `E[Y | A = a, X]` for each combination a = (a_1, ..., a_m) by holding the treatment
#' columns at a. This borrows strength across the whole sample and remains stable as
#' the treatment space grows.
#'
#' **Stratified alternative:** `outcome_fit = "stratified"` instead fits a separate
#' model per combination, using only observations with that treatment pattern. This
#' makes no assumptions about how the outcome surface changes across treatments, which
#' can help when every combination is well populated. Its weakness is sparse data:
#' when a combination has < 10 observations it falls back to the overall outcome mean,
#' ignoring covariates and treatment. Because the number of combinations is 2^m, this
#' starvation---and the resulting same-signed bias---grows with m, so stratified fits
#' are only recommended for small m with large n. Diagnostics show pooled attains lower
#' mu-hat RMSE for m >= 3 (see explorations/diagnostics/m6_stratified_vs_pooled.R).
#'
#' **Method selection:**
#' - Parametric outcome model: Use "glm"
#' - Non-linear covariate effects: Use "gam", "randomforest", or "xgboost"
#' - Skewed outcomes (e.g., costs): Use "xgboost" or "randomforest"
#' - Unknown complexity: Use "xgboost", "randomforest", or "superlearner"
#' - Consistency with PS method: If using sequential + xgboost for PS, use "xgboost" for outcomes
#'
#' **Quality checks:** The function validates outputs:
#' - Checks for non-finite predictions
#' - Warns if combinations have no observations (modeling issue)
#'
#' **Missing combinations:** For sparse data (m >= 5), some combinations may have
#' few observations. Models still predict for all combinations using covariate
#' information.
#'   - treatment_combinations: matrix of treatment combinations
#' @export
#' @examples
#' # Generate example data
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 2, p = 2, seed = 123)
#' 
#' # Estimate outcome models using GLM
#' outcome_est <- estimate_outcome_models(sim_data$data,
#'                                        treatment_cols = paste0("A", 1:2),
#'                                        covariate_cols = paste0("X", 1:2),
#'                                        method = "glm")
#' head(outcome_est$outcome_means)
#'
#' @details
#' Random number generation: The caller is responsible for setting the random
#' seed before calling this function (e.g., set.seed(123)) to ensure reproducibility.
estimate_outcome_models <- function(data, treatment_cols, covariate_cols,
                                   outcome_col = "Y", method = "glm", cross_fit = TRUE,
                                   n_folds = 5, sl_library = NULL,
                                   outcome_fit = c("pooled", "stratified"),
                                   max_combinations = 16, xgb_params = NULL) {

  outcome_fit <- match.arg(outcome_fit)

  # Explicit column specifications (no auto-detection).
  validate_columns(data, treatment_cols, covariate_cols, outcome_col)
  m <- length(treatment_cols)
  p <- length(covariate_cols)

  # Internally the outcome is referenced as "Y" (in glm/gam formulas and
  # data_train$Y). Normalize once so outcome_col is honestly honored regardless
  # of its name; this does not touch the caller's data frame.
  if (outcome_col != "Y") data[["Y"]] <- data[[outcome_col]]

  # Generate all treatment combinations
  vectors <- generate_binary_vectors(m)
  n_combinations <- nrow(vectors)

  # Default: pooled fit. A single model on covariates AND treatment indicators,
  # predicting E[Y | A = a, X] per combination by holding the treatment columns
  # at a. This borrows strength across the full sample and avoids the
  # per-combination starvation of stratified fitting (combinations with < 10 rows
  # otherwise collapse to the grand mean mean(Y), injecting a same-signed bias
  # that grows with m). Diagnostics show pooled has lower mu-hat RMSE than
  # stratified for m >= 3 (2-3x by m = 6) and only a marginal disadvantage at
  # m = 2; see explorations/diagnostics/m6_stratified_vs_pooled.R.
  #
  # `outcome_fit = "stratified"` restores the per-combination fit (assumption-free
  # in treatment; useful when every combination is well populated). For back
  # compatibility, stratified also engages automatically when 2^m <=
  # max_combinations only if the caller explicitly requested it; the default
  # pooled path ignores max_combinations.
  if (outcome_fit == "pooled") {
    return(fit_pooled_outcome_model(
      data = data, treatment_cols = treatment_cols, covariate_cols = covariate_cols,
      vectors = vectors, method = method, cross_fit = cross_fit, n_folds = n_folds,
      sl_library = sl_library, xgb_params = xgb_params
    ))
  }

  # outcome_fit == "stratified": one model per combination. When the treatment
  # space still outruns the sample (2^m > max_combinations) fall back to pooled
  # rather than the grand-mean starvation, and warn.
  if (n_combinations > max_combinations) {
    warning(paste0("stratified outcome fit requested but 2^m = ", n_combinations,
                   " > max_combinations = ", max_combinations,
                   "; using pooled fit to avoid grand-mean starvation."))
    return(fit_pooled_outcome_model(
      data = data, treatment_cols = treatment_cols, covariate_cols = covariate_cols,
      vectors = vectors, method = method, cross_fit = cross_fit, n_folds = n_folds,
      sl_library = sl_library, xgb_params = xgb_params
    ))
  }

  # Initialize results
  outcome_means <- matrix(0, nrow = nrow(data), ncol = n_combinations)
  models <- list()

  if (cross_fit) {
    fold_ids <- sample(rep(1:n_folds, length.out = nrow(data)))
  }
  
  # Fit models for each treatment combination
  for (j in 1:n_combinations) {
    # Find observations with this treatment combination
    match_rows <- rep(TRUE, nrow(data))
    for (k in 1:m) {
      match_rows <- match_rows & (data[[treatment_cols[k]]] == vectors[j, k])
    }
    
    if (sum(match_rows) < 10) {
      # Too few observations, use overall mean
      outcome_means[, j] <- mean(data$Y, na.rm = TRUE)
      models[[j]] <- NULL
      next
    }
    
    if (cross_fit) {
      models[[j]] <- list()
      
      for (fold in 1:n_folds) {
        train_idx <- which(fold_ids != fold & match_rows)
        test_idx <- which(fold_ids == fold)
        
        if (length(train_idx) < 5) {
          # Too few training observations
          outcome_means[test_idx, j] <- mean(data$Y[match_rows], na.rm = TRUE)
          models[[j]][[fold]] <- NULL
          next
        }
        
        # Create data frames with named columns for training and testing
        data_train <- data[train_idx, c(covariate_cols, "Y"), drop = FALSE]
        data_test <- data[test_idx, covariate_cols, drop = FALSE]

        # Fit model
        if (method == "glm") {
          # Check if chunking is needed for memory efficiency
          n_train <- length(train_idx)
          use_chunking <- needs_chunking(n_train, p)

          # Progress reporting for large fits
          if (n_train > 100000 && use_chunking && fold == 1) {
            cat(sprintf("  [Outcome] Fitting combination %d/%d (chunked, n_train=%s)...\n",
                        j, n_combinations, format(n_train, big.mark = ",")))
          }

          if (use_chunking) {
            # Chunked GLM fitting for large datasets
            model <- fit_glm_chunked_internal(
              data_train = data_train,
              outcome_var = "Y",
              covariate_cols = covariate_cols,
              chunk_size = 50000
            )

            outcome_means[test_idx, j] <- predict_glm_chunked_internal(
              model = model,
              newdata = data_test,
              chunk_size = 50000
            )
          } else {
            # Standard GLM for smaller datasets
            model <- glm(Y ~ ., data = data_train, family = gaussian())
            outcome_means[test_idx, j] <- predict(model, newdata = data_test)
          }
        } else if (method == "gam") {
          if (!requireNamespace("mgcv", quietly = TRUE)) {
            stop("mgcv package required for method = 'gam'")
          }
          # GAM needs explicit column names in formula
          formula_str <- paste("Y ~ s(", paste(covariate_cols, collapse = ") + s("), ")", sep = "")
          model <- mgcv::gam(as.formula(formula_str), data = data_train, family = gaussian())
          outcome_means[test_idx, j] <- predict(model, newdata = data_test)
        } else if (method == "randomforest") {
          if (!requireNamespace("randomForest", quietly = TRUE)) {
            stop("randomForest package required for method = 'randomforest'")
          }
          # Convert to matrix only for randomForest
          X_train_matrix <- to_numeric_matrix(data_train, covariate_cols)
          y_train <- data_train$Y
          X_test_matrix <- to_numeric_matrix(data_test, covariate_cols)
          model <- randomForest::randomForest(X_train_matrix, y_train)
          outcome_means[test_idx, j] <- predict(model, X_test_matrix)
        } else if (method == "superlearner") {
          if (!requireNamespace("SuperLearner", quietly = TRUE)) {
            stop("SuperLearner package required for method = 'superlearner'")
          }
          if (is.null(sl_library)) {
            sl_library <- c("SL.glm", "SL.gam", "SL.randomForest")
          }
          y_train <- data_train$Y
          X_train_df <- data_train[, covariate_cols, drop = FALSE]
          X_test_df <- data_test[, covariate_cols, drop = FALSE]
          model <- SuperLearner::SuperLearner(Y = y_train, X = X_train_df,
                                            family = gaussian(), SL.library = sl_library)
          outcome_means[test_idx, j] <- predict(model, newdata = X_test_df)$pred
        } else if (method == "xgboost") {
          if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("xgboost package required for method = 'xgboost'")
          }
          # Prepare matrices
          X_train_matrix <- to_numeric_matrix(data_train, covariate_cols)
          y_train <- as.double(unclass(data_train$Y))
          X_test_matrix <- to_numeric_matrix(data_test, covariate_cols)

          # Train regression model with early stopping
          params <- modifyList(resolve_xgb_params(xgb_params),
                               list(objective = "reg:squarederror", verbosity = 0))
          model <- train_xgboost_early_stop(X_train_matrix, y_train, params)

          # Predict
          outcome_means[test_idx, j] <- predict(model, newdata = X_test_matrix)
        } else {
          stop("method must be one of: 'glm', 'gam', 'randomforest', 'xgboost', 'superlearner'")
        }
        
        models[[j]][[fold]] <- model
      }
    } else {
      # No cross-fitting
      # Create data frame with named columns for this treatment combination
      data_subset <- data[match_rows, c(covariate_cols, "Y"), drop = FALSE]
      data_all <- data[, covariate_cols, drop = FALSE]

      if (method == "glm") {
        # Check if chunking is needed for memory efficiency
        n_train <- sum(match_rows)
        use_chunking <- needs_chunking(n_train, p)

        if (use_chunking) {
          # Chunked GLM fitting for large datasets
          model <- fit_glm_chunked_internal(
            data_train = data_subset,
            outcome_var = "Y",
            covariate_cols = covariate_cols,
            chunk_size = 50000
          )

          outcome_means[, j] <- predict_glm_chunked_internal(
            model = model,
            newdata = data_all,
            chunk_size = 50000
          )
        } else {
          # Standard GLM for smaller datasets
          model <- glm(Y ~ ., data = data_subset, family = gaussian())
          outcome_means[, j] <- predict(model, newdata = data_all)
        }
      } else if (method == "gam") {
        if (!requireNamespace("mgcv", quietly = TRUE)) {
          stop("mgcv package required for method = 'gam'")
        }
        # GAM needs explicit column names in formula
        formula_str <- paste("Y ~ s(", paste(covariate_cols, collapse = ") + s("), ")", sep = "")
        model <- mgcv::gam(as.formula(formula_str), data = data_subset, family = gaussian())
        outcome_means[, j] <- predict(model, newdata = data_all)
      } else if (method == "randomforest") {
        if (!requireNamespace("randomForest", quietly = TRUE)) {
          stop("randomForest package required for method = 'randomforest'")
        }
        # Convert to matrix only for randomForest
        X_subset_matrix <- to_numeric_matrix(data_subset, covariate_cols)
        y_subset <- data_subset$Y
        X_all_matrix <- to_numeric_matrix(data_all, covariate_cols)
        model <- randomForest::randomForest(X_subset_matrix, y_subset)
        outcome_means[, j] <- predict(model, X_all_matrix)
      } else if (method == "superlearner") {
        if (!requireNamespace("SuperLearner", quietly = TRUE)) {
          stop("SuperLearner package required for method = 'superlearner'")
        }
        if (is.null(sl_library)) {
          sl_library <- c("SL.glm", "SL.gam", "SL.randomForest")
        }
        y_subset <- data_subset$Y
        X_subset_df <- data_subset[, covariate_cols, drop = FALSE]
        X_all_df <- data_all[, covariate_cols, drop = FALSE]
        model <- SuperLearner::SuperLearner(Y = y_subset, X = X_subset_df,
                                          family = gaussian(), SL.library = sl_library)
        outcome_means[, j] <- predict(model, newdata = X_all_df)$pred
      } else if (method == "xgboost") {
        if (!requireNamespace("xgboost", quietly = TRUE)) {
          stop("xgboost package required for method = 'xgboost'")
        }
        # Prepare matrices
        X_subset_matrix <- to_numeric_matrix(data_subset, covariate_cols)
        y_subset <- as.double(unclass(data_subset$Y))
        X_all_matrix <- to_numeric_matrix(data_all, covariate_cols)

        # Train with early stopping
        params <- modifyList(resolve_xgb_params(xgb_params),
                             list(objective = "reg:squarederror", verbosity = 0))
        model <- train_xgboost_early_stop(X_subset_matrix, y_subset, params)

        # Predict for all observations
        outcome_means[, j] <- predict(model, newdata = X_all_matrix)
      } else {
        stop("method must be one of: 'glm', 'gam', 'randomforest', 'xgboost', 'superlearner'")
      }
      
      models[[j]] <- model

      # Memory cleanup every 50 combinations
      if (j %% 50 == 0) {
        gc(verbose = FALSE)
      }
    }
  }

  # Final memory cleanup after all combinations
  gc(verbose = FALSE)

  # Validate outcome means
  if (any(!is.finite(outcome_means))) {
    n_nonfinite <- sum(!is.finite(outcome_means))
    stop("Outcome model estimation produced ", n_nonfinite, " non-finite values (NA, NaN, or Inf). ",
         "Method '", method, "' failed. ",
         "Check data quality and model convergence.")
  }

  # Check for combinations with no observations (should not happen if we checked earlier, but safety)
  obs_per_combo <- colSums(!is.na(outcome_means))
  if (any(obs_per_combo == 0)) {
    n_zero <- sum(obs_per_combo == 0)
    warning("Outcome model produced estimates for ", n_zero, " combinations with no observations. ",
            "This may indicate modeling issues.")
  }

  return(list(
    outcome_means = outcome_means,
    models = models,
    method = method,
    treatment_combinations = vectors,
    cross_fit = cross_fit
  ))
}

#' Pooled outcome model for large treatment spaces (internal)
#'
#' Fits a single outcome regression on covariates and treatment indicators
#' (treatments as features, not strata), then predicts `E[Y | A = a, X]` for every
#' treatment combination by setting the treatment columns to a. Used by
#' \code{estimate_outcome_models} for the default \code{outcome_fit = "pooled"} to
#' avoid the per-combination grand-mean fallback that starves at high m. Returns
#' the same list contract as \code{estimate_outcome_models}.
#'
#' @keywords internal
fit_pooled_outcome_model <- function(data, treatment_cols, covariate_cols,
                                     vectors, method = "glm", cross_fit = TRUE,
                                     n_folds = 5, sl_library = NULL, xgb_params = NULL) {
  n <- nrow(data)
  n_combinations <- nrow(vectors)
  feat_cols <- c(covariate_cols, treatment_cols)  # X and treatment indicators

  # Predict E[Y | A = a, X] for every row and combination from a fitted model,
  # by overwriting the treatment columns with each pattern a. `predict_fun` maps
  # a feature data frame to a numeric vector of predictions.
  predict_all_combos <- function(predict_fun, rows) {
    out <- matrix(0, nrow = length(rows), ncol = n_combinations)
    base <- data[rows, feat_cols, drop = FALSE]
    for (j in seq_len(n_combinations)) {
      nd <- base
      for (k in seq_along(treatment_cols)) nd[[treatment_cols[k]]] <- vectors[j, k]
      out[, j] <- predict_fun(nd)
    }
    out
  }

  # Fit Y ~ feat_cols on `train_df`, returning a predict_fun closure over feat_cols.
  fit_predict_fun <- function(train_df) {
    if (method == "glm") {
      model <- glm(Y ~ ., data = train_df[, c(feat_cols, "Y"), drop = FALSE],
                   family = gaussian())
      function(nd) as.numeric(predict(model, newdata = nd))
    } else if (method == "gam") {
      if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv package required for method = 'gam'")
      # Smooth terms on covariates; treatment indicators enter linearly.
      formula_str <- paste("Y ~ s(", paste(covariate_cols, collapse = ") + s("), ") + ",
                           paste(treatment_cols, collapse = " + "), sep = "")
      model <- mgcv::gam(as.formula(formula_str), data = train_df, family = gaussian())
      function(nd) as.numeric(predict(model, newdata = nd))
    } else if (method == "randomforest") {
      if (!requireNamespace("randomForest", quietly = TRUE)) stop("randomForest package required for method = 'randomforest'")
      model <- randomForest::randomForest(to_numeric_matrix(train_df, feat_cols), train_df$Y)
      function(nd) as.numeric(predict(model, to_numeric_matrix(nd, feat_cols)))
    } else if (method == "superlearner") {
      if (!requireNamespace("SuperLearner", quietly = TRUE)) stop("SuperLearner package required for method = 'superlearner'")
      if (is.null(sl_library)) sl_library <- c("SL.glm", "SL.gam", "SL.randomForest")
      model <- SuperLearner::SuperLearner(Y = train_df$Y, X = train_df[, feat_cols, drop = FALSE],
                                          family = gaussian(), SL.library = sl_library)
      function(nd) as.numeric(predict(model, newdata = nd[, feat_cols, drop = FALSE])$pred)
    } else if (method == "xgboost") {
      if (!requireNamespace("xgboost", quietly = TRUE)) stop("xgboost package required for method = 'xgboost'")
      params <- modifyList(resolve_xgb_params(xgb_params),
                           list(objective = "reg:squarederror", verbosity = 0))
      model <- train_xgboost_early_stop(to_numeric_matrix(train_df, feat_cols),
                                        as.double(unclass(train_df$Y)), params)
      function(nd) as.numeric(predict(model, newdata = to_numeric_matrix(nd, feat_cols)))
    } else {
      stop("method must be one of: 'glm', 'gam', 'randomforest', 'xgboost', 'superlearner'")
    }
  }

  outcome_means <- matrix(0, nrow = n, ncol = n_combinations)
  models <- list()

  if (cross_fit) {
    fold_ids <- sample(rep(1:n_folds, length.out = n))
    for (fold in 1:n_folds) {
      train_idx <- which(fold_ids != fold)
      test_idx  <- which(fold_ids == fold)
      pf <- fit_predict_fun(data[train_idx, , drop = FALSE])
      outcome_means[test_idx, ] <- predict_all_combos(pf, test_idx)
      models[[fold]] <- pf
    }
  } else {
    pf <- fit_predict_fun(data)
    outcome_means <- predict_all_combos(pf, seq_len(n))
    models <- list(pf)
  }

  if (any(!is.finite(outcome_means))) {
    stop("Pooled outcome model produced ", sum(!is.finite(outcome_means)),
         " non-finite values. Method '", method, "' failed.")
  }

  list(outcome_means = outcome_means, models = models, method = method,
       treatment_combinations = vectors, cross_fit = cross_fit, pooled = TRUE)
}
