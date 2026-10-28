#' @file estimation.R
#' @title Estimation Functions for Vector Incremental Effects
#' @description
#' Main estimation functions for computing psi_delta = `E[Y^{Q_delta}]`, the
#' counterfactual mean under incremental interventions. Implements outcome
#' regression, IPW, and one-step (EIF-based) estimators with full inference.
#' Also includes an experimental TMLE-style estimator with no supporting theory.
#' @author Denis Agniel

#' Compute h function value for divergence-based tilt
#'
#' @description
#' Computes the h(a, lambda) function for different divergence measures.
#' This function is used in the tilting mechanism for Q_delta distributions.
#'
#' @param a Numeric vector, binary treatment combination (length m)
#' @param lambda Numeric vector, Lagrange multipliers (length m+1: intercept + m marginal terms)
#' @param divergence Character, divergence type: "kl", "reverse_kl", "chi2", "hellinger"
#'
#' @return Numeric, h(a, lambda) value for the specified divergence
#'
#' @details
#' For a linear term `z = lambda[1] + a' * lambda[2:(m+1)]`:
#' \itemize{
#'   \item KL divergence: h(z) = exp(z)
#'   \item Reverse KL: h(z) = -1/z (with safeguarding for z < 0, see Table 1 in paper)
#'   \item Chi-squared: h(z) = z/2 + 1
#'   \item Hellinger: h(z) = 1/(1-z)^2 (with domain checks)
#' }
#'
#' Domain violations raise an error rather than saturating (unlike the hot
#' Newton-solver/EIF path in \code{divergence_h_pair()}, which this function
#' otherwise delegates to for the actual h-value): this is a user-facing
#' helper, so failing loudly on a genuinely infeasible \code{(a, lambda)}
#' pair is more useful than silently returning a saturated value.
#'
#' @keywords internal
compute_h_function_value <- function(a, lambda, divergence) {
  # Compute linear term: lambda[1] + a' * lambda[2:end]
  m <- length(a)
  if (length(lambda) != m + 1) {
    stop("lambda must have length m+1 (intercept + m marginal terms)")
  }

  linear_term <- lambda[1] + sum(a * lambda[2:(m+1)])

  # Explicit, fail-loud domain checks (kept stricter than the canonical
  # primitive's internal clip buffer, so these always trigger before any
  # clipping would): DOMAIN_VALIDATION_STRICT_TOL = 1e-8 vs. the canonical
  # primitive's 1e-12 buffers.
  if (divergence == "reverse_kl" && linear_term >= -DOMAIN_VALIDATION_STRICT_TOL) {
    stop("Invalid domain for reverse KL: linear_term >= 0 (must be negative)")
  }
  if (divergence == "hellinger" && linear_term >= 1 - DOMAIN_VALIDATION_STRICT_TOL) {
    stop("Invalid domain for Hellinger: linear_term >= 1")
  }

  # Canonical h-function (see divergence_h_functions.R) -- single source of
  # truth, shared with the Newton solver and EIF computation.
  divergence_h_pair(linear_term, divergence)$h
}

#' Compute plugin estimate of `E[Y^{Q_delta}]`
#'
#' @description
#' Calculates the causal estimand `E[Y^{Q_delta}]` using the identification
#' results from the paper. Implements three estimation approaches: outcome
#' regression, inverse probability weighting (IPW), and the efficient-influence-
#' function-based estimator (named "doubly_robust" for backward compatibility;
#' see Details for why that name is a misnomer).
#' This is typically called internally by estimate_qdelta_mean() but can
#' be used directly for custom workflows.
#'
#' @param data Data frame with treatment (A1, ..., Am), covariate (X1, ..., Xp),
#'   and outcome (Y) columns
#' @param delta Numeric vector of length m, target marginal changes (odds ratios)
#' @param treatment_cols Character vector of treatment column names (required).
#' @param covariate_cols Character vector of covariate column names (required).
#' @param outcome_col Character scalar, outcome column name. Default "Y".
#' @param qdelta_results List of length n, per-observation Q_delta solver results.
#'   Each element from compute_qdelta_per_obs() containing q_delta, lambda, converged.
#'   If NULL, will compute automatically (slower).
#' @param outcome_means Numeric matrix (n x 2^m), estimated `E[Y|A=a,X]` for each
#'   observation and treatment combination. Rows are observations, columns are
#'   treatment combinations in binary order.
#' @param propensity_scores Numeric matrix (n x 2^m), estimated P(A=a|X) for each
#'   observation and treatment combination. Must sum to 1 across columns.
#' @param vectors Matrix (2^m x m), treatment combinations in binary order.
#'   Optional if m is provided. Generated automatically if NULL.
#' @param m Integer, number of treatments. Only needed if vectors = NULL.
#'   Defaults to length(treatment_cols).
#' @param method Character, estimation approach:
#'   - "outcome_regression": `E_n[sum_a mu(a,X) * p_0(a|X) * h(a, lambda)]`
#'   - "ipw": `E_n[Y * I(A=a) * h(a, lambda)]`
#'   - "doubly_robust": Uses the full EIF (recommended; the name is kept for
#'     backward compatibility but is a misnomer -- see Details)
#' @param divergence Character, divergence type ("kl", "reverse_kl", "chi2", "hellinger")
#' @param return Character, what to return: "mean" (default, returns scalar) or
#'   "vector" (returns observation-level contributions)
#'
#' @return If return = "mean": Numeric scalar, estimated `E[Y^{Q_delta}]`.
#'   If return = "vector": Numeric vector (length n) of observation-level contributions.
#'
#' @details
#' The three methods have different robustness properties:
#' - **outcome_regression**: Consistent if outcome models correct, ignores propensity misspecification
#' - **ipw**: Consistent if propensity models correct, ignores outcome misspecification
#' - **doubly_robust** (recommended): Uses the full EIF. Requires **both** nuisances
#'   to be consistent, but nuisance error enters only through a second-order
#'   product-rate remainder
#'   \code{R_n = max_a||pihat_a - pi_a||_2 * max_a||muhat_a - mu_a||_2 + max_a||pihat_a - pi_a||_2^2},
#'   so a slowly-converging outcome model is tolerated when the propensity model
#'   converges quickly. The method name is retained for backward compatibility
#'   and is a misnomer: this estimator is **not** doubly robust. It is not
#'   consistent when the propensity model is inconsistent but the outcome model
#'   is correct, because the estimand depends on the propensity score through
#'   the divergence tilt and not only through inverse-probability weighting.
#'   Note also that the squared term requires
#'   \code{||pihat_a - pi_a||_2 = o_p(n^(-1/4))} even when the outcome model is
#'   exactly correct. See Theorem 1 of the accompanying paper.
#'
#' For the "doubly_robust" method, the function uses the full EIF formula with
#' Omega matrix correction, which is the estimator analyzed in the paper.
#'
#' @seealso
#' \code{\link{estimate_qdelta_mean}} for end-to-end estimation with inference,
#' \code{\link{compute_qdelta_per_obs}} for per-observation Q_delta computation
#' @export
#' @examples
#' # Generate example data
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 2, p = 2, seed = 123)
#'
#' # Estimate nuisance parameters
#' tcols <- paste0("A", 1:2)
#' xcols <- paste0("X", 1:2)
#' ps_est <- estimate_propensity_scores(sim_data$data, treatment_cols = tcols,
#'                                      covariate_cols = xcols, method = "glm")
#' outcome_est <- estimate_outcome_models(sim_data$data, treatment_cols = tcols,
#'                                        covariate_cols = xcols, method = "glm")
#'
#' # Compute per-observation Q_delta
#' vectors <- generate_binary_vectors(2)
#' delta <- c(1.5, 2.0)
#' qdelta_results <- vector("list", nrow(sim_data$data))
#' for (i in 1:nrow(sim_data$data)) {
#'   qdelta_results[[i]] <- compute_qdelta_per_obs(
#'     ps_est$propensity_scores[i, ], delta, "kl"
#'   )
#' }
#'
#' # Compute psi_delta
#' psi_est <- compute_psi_delta(sim_data$data, delta,
#'                             treatment_cols = tcols, covariate_cols = xcols,
#'                             qdelta_results = qdelta_results,
#'                             outcome_means = outcome_est$outcome_means,
#'                             propensity_scores = ps_est$propensity_scores,
#'                             vectors = vectors, method = "doubly_robust")
compute_psi_delta <- function(data, delta, treatment_cols, covariate_cols,
                             outcome_col = "Y", qdelta_results = NULL, outcome_means,
                             propensity_scores, vectors = NULL, m = NULL,
                             method = "doubly_robust", divergence = "kl",
                             return = "mean") {

  validate_columns(data, treatment_cols, covariate_cols, outcome_col)

  # Handle vectors parameter (combinatorics; m is the number of treatments)
  if (is.null(vectors)) {
    if (is.null(m)) {
      m <- length(treatment_cols)
    }
    vectors <- generate_binary_vectors(m)
  } else {
    m <- ncol(vectors)
  }
  if (is.null(qdelta_results)) {
    qdelta_results = list()
    for (i in 1:nrow(data)) {
      qdelta_results[[i]] <- compute_qdelta_per_obs(
        propensity_scores[i, ], delta, divergence, tol = 1e-8
      )
    }
  }

  n <- nrow(data)
  n_combinations <- nrow(vectors)

  # Validate dimensions
  if (!is.matrix(outcome_means) || nrow(outcome_means) != n || ncol(outcome_means) != n_combinations) {
    stop("outcome_means must be an n x 2^m matrix (", n, " x ", n_combinations, "). ",
         "Got: ", paste(dim(outcome_means), collapse = " x "), ".", call. = FALSE)
  }

  if (!is.matrix(propensity_scores) || nrow(propensity_scores) != n || ncol(propensity_scores) != n_combinations) {
    stop("propensity_scores must be an n x 2^m matrix (", n, " x ", n_combinations, "). ",
         "Got: ", paste(dim(propensity_scores), collapse = " x "), ".", call. = FALSE)
  }

  # Validate qdelta_results
  if (length(qdelta_results) != n) {
    stop("qdelta_results must have length equal to number of observations (", n, ").", call. = FALSE)
  }

  if (method == "outcome_regression") {
    # Outcome regression path: E[sum_a mu_a(X) * p_0(a|X) * h(a, lambda)]
    # Compute h values for all observations and combinations
    h_values_matrix <- matrix(0, nrow = n, ncol = n_combinations)
    for (i in 1:n) {
      lambda_i <- qdelta_results[[i]]$lambda
      # Vectorized over treatment combinations
      linear_terms <- lambda_i[1] + vectors %*% lambda_i[2:(m+1)]
      if (divergence == "kl") {
        h_values_matrix[i, ] <- exp(linear_terms)
      } else {
        # For other divergences, compute individually (less common)
        for (j in 1:n_combinations) {
          h_values_matrix[i, j] <- compute_h_function_value(vectors[j, ], lambda_i, divergence)
        }
      }
    }
    # Vectorized sum: sum over all i and j
    psi_est <- sum(outcome_means * propensity_scores * h_values_matrix) / n

  } else if (method == "ipw") {
    # IPW path: E[Y * sum_a I(A=a) * h(a, lambda)]
    # Vectorized: extract treatment assignments as matrix
    A_matrix <- as.matrix(data[, treatment_cols, drop = FALSE])

    # Find combination indices for each observation (vectorized)
    combination_indices <- apply(A_matrix, 1, function(a_i) {
      which(apply(vectors, 1, function(x) all(x == a_i)))
    })

    # Compute h values for each observation's treatment
    h_values <- numeric(n)
    for (i in 1:n) {
      lambda_i <- qdelta_results[[i]]$lambda
      a_i <- A_matrix[i, ]
      h_values[i] <- compute_h_function_value(a_i, lambda_i, divergence)
    }

    psi_est <- sum(data[[outcome_col]] * h_values) / n

  } else if (method == "doubly_robust") {
    # Doubly-robust path: uses the full EIF with Omega matrix correction
    # This is the recommended approach and uses the correct EIF formula from main.tex

    # Extract X, A, Y from data (explicit column specifications)
    X_matrix <- if (length(covariate_cols) > 0) {
      as.matrix(data[, covariate_cols, drop = FALSE])
    } else {
      NULL
    }
    A_matrix <- as.matrix(data[, treatment_cols, drop = FALSE])
    Y_vector <- data[[outcome_col]]

    # Compute full EIF with pre-computed qdelta_results
    eif_result <- efficient_influence_function(
      X = X_matrix,
      A = A_matrix,
      Y = Y_vector,
      propensity_scores = propensity_scores,
      delta = delta,
      outcome_means = outcome_means,
      divergence = divergence,
      qdelta_results = qdelta_results
    )
    gc(verbose = FALSE)  # Clear memory after EIF computation

    # The EIF mean equals the true causal effect when nuisance parameters are correct
    if (return == "mean") {
      psi_est <- mean(eif_result$eif_values)
    } else if (return == "ifvals") {
      psi_est <- eif_result
    } else {
      stop("return must be 'mean' or 'ifvals'", call. = FALSE)
    }

  } else {
    stop("method must be one of: 'outcome_regression', 'ipw', 'doubly_robust'")
  }
  return(psi_est)
}

#' Target outcome model for TMLE
#'
#' Performs the TMLE targeting step by updating the outcome model using a
#' parametric submodel. The clever covariate is computed from the EIF components.
#'
#' @param data tibble, observed data
#' @param outcome_means Matrix, initial outcome means (n × 2^m)
#' @param propensity_scores Matrix, propensity scores (n × 2^m)
#' @param qdelta_results List, per-observation Q_delta results
#' @param vectors Matrix, treatment combinations (2^m × m)
#' @param treatment_cols Character vector of treatment column names.
#' @param outcome_col Character scalar, outcome column name. Default "Y".
#' @param divergence Character, divergence type
#' @return Matrix, targeted outcome means (n × 2^m)
#' @keywords internal
target_outcome_model_tmle <- function(data, outcome_means, propensity_scores,
                                     qdelta_results, vectors, treatment_cols,
                                     outcome_col = "Y", divergence = "kl") {

  n <- nrow(data)
  n_combinations <- nrow(vectors)
  m <- ncol(vectors)

  # Extract treatment assignments (explicit column specifications)
  A_matrix <- as.matrix(data[, treatment_cols, drop = FALSE])
  Y_vector <- data[[outcome_col]]

  # Initialize targeted outcome means
  targeted_outcome_means <- outcome_means

  # For each treatment combination, perform targeting
  for (a_idx in 1:n_combinations) {
    a <- vectors[a_idx, ]

    # Find observations with this treatment combination (vectorized)
    match_rows <- rowSums(A_matrix == matrix(rep(a, each = n), nrow = n, ncol = m)) == m

    if (sum(match_rows) < 2) {
      # Too few observations, skip targeting for this combination
      next
    }

    # Compute clever covariate for this treatment combination
    # The clever covariate H(a, X) is h(a, λ) for each observation (vectorized)
    clever_covariate <- numeric(n)
    for (i in 1:n) {
      lambda_i <- qdelta_results[[i]]$lambda
      clever_covariate[i] <- compute_h_function_value(a, lambda_i, divergence)
    }

    # Parametric submodel: μ*_a(X) = μ_a(X) + ε · H(a, X)
    # Fit GLM with offset = μ_a(X) and covariate = H(a, X)
    # Only use observations with this treatment combination
    mu_a_initial <- outcome_means[match_rows, a_idx]
    H_a <- clever_covariate[match_rows]
    Y_a <- Y_vector[match_rows]

    # Fit targeting model: Y = offset(mu_a) + epsilon * H_a
    # Use glm with offset
    targeting_data <- data.frame(
      Y = Y_a,
      H = H_a,
      offset = mu_a_initial
    )

    # Fit model: Y ~ offset(mu_a) + epsilon * H
    # This is equivalent to: Y - mu_a = epsilon * H
    tryCatch({
      targeting_model <- glm(Y ~ 0 + H + offset(offset),
                            data = targeting_data,
                            family = gaussian())

      epsilon <- coef(targeting_model)["H"]
      if (is.na(epsilon)) {
        epsilon <- 0
      }

      # Update outcome means for all observations (not just those with A=a)
      # μ*_a(X) = μ_a(X) + ε · H(a, X)
      targeted_outcome_means[, a_idx] <- outcome_means[, a_idx] + epsilon * clever_covariate

    }, error = function(e) {
      # If targeting fails, keep original outcome means
      warning("TMLE targeting failed for treatment combination ", a_idx,
              ": ", e$message, ". Using initial outcome model.")
    })
  }

  return(targeted_outcome_means)
}

#' Estimate counterfactual mean under incremental intervention distribution
#'
#' @description
#' This is the primary user-facing function that estimates `E[Y^{Q_delta}]`, the
#' counterfactual mean under the incremental intervention distribution Q_delta.
#' Implements cross-fitted one-step estimation based on the efficient influence
#' function (EIF), the estimator analyzed in the accompanying paper. Also
#' supports an experimental TMLE-style estimator with no supporting theory.
#'
#' @param data Data frame or tibble with columns:
#'   - Treatment columns named A1, A2, ..., Am (binary 0/1)
#'   - Covariate columns named X1, X2, ..., Xp (numeric)
#'   - Outcome column Y (numeric)
#' @param delta Numeric vector of length m specifying target marginal treatment
#'   probabilities as odds ratios. Each element must be positive. Values > 1
#'   increase treatment probability; values < 1 decrease it. E.g., delta = c(1.5, 2.0)
#'   increases P(A1=1) by factor 1.5 and P(A2=1) by factor 2.0.
#' @param treatment_cols Character vector of treatment column names (required).
#'   Length defines m; order must match `delta` and the intended bit order.
#' @param covariate_cols Character vector of covariate column names (required).
#' @param outcome_col Character scalar, outcome column name. Default "Y".
#' @param divergence Character, divergence measure for Q_delta optimization.
#'   Options: "kl" (Kullback-Leibler, recommended), "reverse_kl", "chi2",
#'   "hellinger". Default: "kl". Different divergences have different stability
#'   and boundary behavior.
#' @param ps_method Character, propensity score estimation method.
#'   Options: "glm" (multinomial logistic), "randomforest", "xgboost",
#'   "binary_per_combination" (fits 2^m binary models). Default: "glm".
#' @param binary_method Character, base learner for binary_per_combination method.
#'   Options: "glm", "randomforest", "xgboost". Only used if ps_method = "binary_per_combination".
#' @param outcome_method Character, outcome model estimation method.
#'   Options: "glm", "gam", "randomforest", "superlearner". Default: "glm".
#' @param oracle_propensities Optional numeric matrix (n x 2^m) of known
#'   propensity scores. If supplied together with `oracle_outcomes`, model
#'   fitting is bypassed and these are used directly. Default: NULL.
#' @param oracle_outcomes Optional numeric matrix (n x 2^m) of known outcome
#'   means `E[Y|A=a,X]`. Must be supplied together with `oracle_propensities`.
#'   Default: NULL.
#' @param estimator Character, final estimator type.
#'   Options: "one_step" (recommended, EIF-based, the estimator analyzed in the
#'   paper's Theorem 1), "tmle" (EXPERIMENTAL -- targeting step implemented but
#'   no supporting theory; not analyzed in the paper), "plugin" (deprecated, no
#'   EIF correction). Default: "one_step".
#' @param cross_fit Logical, use cross-fitting for sample splitting to reduce
#'   bias. Recommended for small samples or flexible nuisance estimation.
#'   Default: FALSE.
#' @param n_folds Integer, number of folds for cross-fitting. Only used if
#'   cross_fit = TRUE. Default: 5.
#' @param confidence_level Numeric in (0, 1), confidence level for intervals.
#'   Default: 0.95 (95% confidence intervals).
#' @param bootstrap Logical, use bootstrap for inference instead of asymptotic
#'   variance. More robust but slower. Default: FALSE.
#' @param n_bootstrap Integer, number of bootstrap resamples. Only used if
#'   bootstrap = TRUE. Default: 1000.
#' @param threshold_type Character, active set selection method for handling
#'   boundary propensities. Options:
#'   \itemize{
#'     \item "none" (default): Standard estimator, all treatments always active
#'     \item "vanishing": Vanishing threshold ε_n → 0 with fixed delta (asymptotically optimal)
#'     \item "adaptive": Fixed threshold ε with adaptive delta_k(X) (localized margin)
#'   }
#' @param gamma Numeric, threshold rate parameter for vanishing threshold.
#'   `ε_n = C*n^{-γ}`. Only used if threshold_type = "vanishing".
#'   Default: NULL (auto-selects 0.3).
#' @param epsilon Numeric, fixed threshold for adaptive estimator.
#'   Only used if threshold_type = "adaptive". Default: 0.1.
#' @param scaling_method Character, delta scaling method for adaptive estimator.
#'   Options:
#'   \itemize{
#'     \item "proportional" (smooth scaling)
#'     \item "fixed" (fixed fraction)
#'     \item "zero" (binary)
#'   }
#'   Only used if threshold_type = "adaptive". Default: "proportional".
#' @param inactive_fraction Numeric, fraction of target delta for inactive treatments
#'   with fixed scaling. Only used if threshold_type = "adaptive" and
#'   scaling_method = "fixed". Default: 0.5.
#' @param min_active Integer, minimum number of active treatments required per
#'   observation. Only used if threshold_type != "none". Default: 1.
#' @param contrast Character, estimand specification. Options:
#'   \itemize{
#'     \item "absolute" (default): Estimate `E[Y^{Q_delta}]`, the counterfactual mean
#'     \item "difference": Estimate `E[Y^{Q_delta}] - E[Y]`, the difference from observed mean
#'   }
#'   The "difference" option provides correct inference for the causal contrast via
#'   modified EIF. Default: "absolute".
#' @param prob_infeasible Numeric in `[0, 1]`, maximum proportion of infeasible observations
#'   allowed before erroring. If the proportion of observations with infeasible delta
#'   is less than this threshold, estimation continues using only feasible observations
#'   with a warning. Default: 0 (strict mode: error if even one observation is infeasible).
#'   Example: prob_infeasible = 0.05 allows up to 5% infeasible observations.
#'
#' @param max_qdelta_failure_rate Numeric in `[0, 1]`: the largest share of observations whose
#'   Q_delta solve may fail and fall back to Q_delta = P0 (no shift for that observation, a bias
#'   toward the null). Default `0`: the first failure stops the estimate with the observation index
#'   and the solver's message. A positive value tolerates failures up to `rate * n` observations
#'   (each warned about and counted), then stops. `1` always falls back. The vanishing-threshold
#'   estimator never falls back and always stops.
#' @param qdelta_tol Numeric or `NULL`. Solver tolerance for the per-observation Q_delta problem;
#'   `NULL` (default) uses `QDELTA_TOLERANCE` for every sample size (it no longer relaxes above
#'   100,000 observations, which made the same stratum's estimate depend on whether it crossed that
#'   size).
#' @return List containing:
#' \describe{
#'   \item{estimate}{Numeric, point estimate of `E[Y^{Q_delta}]` or `E[Y^{Q_delta}] - E[Y]`}
#'   \item{se}{Numeric, standard error from EIF or bootstrap}
#'   \item{ci_lower}{Numeric, lower (1-alpha)/2 confidence bound}
#'   \item{ci_upper}{Numeric, upper (1-alpha)/2 confidence bound}
#'   \item{p_value}{Numeric, two-sided p-value for H0: psi_delta = 0}
#'   \item{eif_values}{Numeric vector (length n), centered EIF values}
#'   \item{qdelta_results}{List (length n), per-observation Q_delta solver results}
#'   \item{convergence_info}{List with converged (logical vector), iterations (integer vector)}
#'   \item{convergence_summary}{List with n_converged, convergence_rate, mean_iterations}
#'   \item{nuisance_estimates}{List with propensity_scores (n x 2^m matrix), outcome_means (n x 2^m matrix), models}
#'   \item{contrast}{Character, estimand type ("absolute" or "difference")}
#'   \item{n_original}{Integer, original sample size before filtering infeasible observations}
#'   \item{n_infeasible_excluded}{Integer, number of observations excluded due to infeasibility}
#'   \item{prop_infeasible_excluded}{Numeric, proportion of observations excluded (0 if none)}
#'   \item{prob_infeasible_threshold}{Numeric, threshold used for infeasibility filtering}
#' }
#'
#' @details
#' The function implements the estimation procedure described in the paper:
#' 1. Estimates nuisance functions: propensity scores P(A|X) and outcome models `E[Y|A,X]`
#' 2. For each observation, solves for Q_delta distribution using constrained optimization
#' 3. Computes the one-step estimate from the EIF (or the experimental TMLE-style estimator)
#' 4. Provides valid inference via influence function variance or bootstrap
#'
#' **Threshold types (handling boundary propensities):**
#'
#' - **"none"** (default): Standard estimator assuming all treatments active everywhere.
#'   Use when propensities are well away from 0 and 1.
#'
#' - **"vanishing"**: Vanishing threshold ε_n → 0 with active sets `K_{ε_n}(X)`.
#'   Asymptotically optimal under global margin condition. Automatically selects
#'   rates (gamma, alpha) or accepts user values. Use for large samples with
#'   theory-driven analysis.
#'
#' - **"adaptive"**: Fixed threshold ε with observation-specific delta_k(X).
#'   Uses localized margin condition. Simpler than vanishing threshold with
#'   potentially better finite sample stability. User controls epsilon (e.g., 0.1)
#'   and scaling_method. Use for finite samples or when wanting direct control.
#'
#' **When to use threshold types:**
#' - Well-behaved propensities: threshold_type = "none"
#' - Some boundary propensities, large n: threshold_type = "vanishing"
#' - Some boundary propensities, finite n: threshold_type = "adaptive"
#' - Very extreme propensities: Consider data quality or eligibility stratification
#'
#' **Recommended settings:**
#' - For most applications: estimator = "one_step", divergence = "kl", threshold_type = "none"
#' - For boundary propensities: threshold_type = "adaptive", epsilon = 0.1
#' - For high-dimensional treatments (m >= 5): ps_method = "binary_per_combination"
#' - For flexible estimation: use randomForest/xgboost methods with cross_fit = TRUE
#'
#' **Random number generation:** The caller controls RNG state. Set seed before
#' calling if using cross-fitting or bootstrap: set.seed(123).
#'
#' **Feasibility checking:** The function automatically checks if delta is feasible
#' for all observations (when threshold_type = "none"). If infeasible, suggests
#' alternatives via suggest_feasible_delta(). When using threshold types, active
#' set selection handles boundary issues automatically.
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate example data with 2 treatments
#' set.seed(123)
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 2, p = 3)
#' tcols <- paste0("A", 1:2)
#' xcols <- paste0("X", 1:3)
#'
#' # Basic estimation with default settings (one-step, GLM)
#' result <- estimate_qdelta_mean(
#'   data = sim_data$data,
#'   delta = c(1.5, 2.0),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   divergence = "kl"
#' )
#'
#' # Print results
#' print(result$estimate)
#' print(result$se)
#' print(c(result$ci_lower, result$ci_upper))
#'
#' # Check convergence
#' print(result$convergence_summary$convergence_rate)  # Should be close to 1
#'
#' # With flexible estimation and cross-fitting
#' result2 <- estimate_qdelta_mean(
#'   data = sim_data$data,
#'   delta = c(1.5, 2.0),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   ps_method = "randomforest",
#'   outcome_method = "randomforest",
#'   cross_fit = TRUE,
#'   n_folds = 5
#' )
#'
#' # High-dimensional treatments (m = 5)
#' sim_highdim <- simulate_vector_treatment_data(n = 1000, m = 5, p = 10)
#' result3 <- estimate_qdelta_mean(
#'   data = sim_highdim$data,
#'   delta = rep(1.2, 5),
#'   treatment_cols = paste0("A", 1:5),
#'   covariate_cols = paste0("X", 1:10),
#'   ps_method = "binary_per_combination",
#'   binary_method = "xgboost"
#' )
#'
#' # Vanishing threshold estimator (robust to boundary propensities)
#' result_vt <- estimate_qdelta_mean(
#'   data = sim_data$data,
#'   delta = c(1.5, 2.0),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   threshold_type = "vanishing",
#'   gamma = 0.3   # Optional: auto-selected if NULL
#' )
#'
#' # Adaptive estimator (fixed epsilon, adaptive delta)
#' result_adaptive <- estimate_qdelta_mean(
#'   data = sim_data$data,
#'   delta = c(1.5, 2.0),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   threshold_type = "adaptive",
#'   epsilon = 0.1,
#'   scaling_method = "proportional"
#' )
#'
#' # Compare all three approaches
#' cat("Standard:  ", result$estimate, "\n")
#' cat("Vanishing: ", result_vt$estimate, "\n")
#' cat("Adaptive:  ", result_adaptive$estimate, "\n")
#'
#' # Estimate difference from observed mean (causal contrast)
#' result_diff <- estimate_qdelta_mean(
#'   data = sim_data$data,
#'   delta = c(1.5, 2.0),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   contrast = "difference"
#' )
#'
#' # The estimate is now E[Y^{Q_delta}] - E[Y] with correct inference
#' cat("Difference from observed: ", result_diff$estimate, "\n")
#' cat("95% CI: [", result_diff$ci_lower, ", ", result_diff$ci_upper, "]\n")
#'
#' # Allow up to 5% infeasible observations
#' result_flexible <- estimate_qdelta_mean(
#'   data = sim_data$data,
#'   delta = c(1.8, 2.5),  # Potentially aggressive delta
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   prob_infeasible = 0.05  # Continue if < 5% infeasible
#' )
#'
#' # Check how many observations were excluded
#' cat("Excluded:", result_flexible$n_infeasible_excluded, "observations\n")
#' cat("Proportion:", result_flexible$prop_infeasible_excluded, "\n")
#' }
estimate_qdelta_mean <- function(data, delta, treatment_cols, covariate_cols,
                               outcome_col = "Y", divergence = "kl",
                               ps_method = "glm", binary_method = "glm",
                               outcome_method = "glm",
                               oracle_propensities = NULL,
                               oracle_outcomes = NULL,
                               estimator = "one_step", cross_fit = TRUE, n_folds = 5,
                               confidence_level = 0.95, bootstrap = FALSE,
                               n_bootstrap = 1000,
                               threshold_type = c("none", "vanishing", "adaptive"),
                               gamma = NULL,
                               epsilon = 0.1, scaling_method = "proportional",
                               inactive_fraction = 0.5, min_active = 1,
                               contrast = c("absolute", "difference"),
                               prob_infeasible = 0, xgb_params = NULL,
    max_qdelta_failure_rate = 0, qdelta_tol = NULL) {

  # Match arguments
  threshold_type <- match.arg(threshold_type)
  contrast <- match.arg(contrast)

  # Input validation
  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a positive numeric vector")
  }

  # Explicit column specifications (no auto-detection).
  validate_columns(data, treatment_cols, covariate_cols, outcome_col)
  m <- length(treatment_cols)

  if (length(delta) != m) {
    stop("delta length must equal number of treatments")
  }

  # Oracle mode: Use provided nuisances if supplied
  if (!is.null(oracle_propensities) && !is.null(oracle_outcomes)) {
    message("Using oracle nuisances (bypassing model fitting)")

    # Validate dimensions
    if (!is.matrix(oracle_propensities) ||
        nrow(oracle_propensities) != nrow(data) ||
        ncol(oracle_propensities) != 2^m) {
      stop("oracle_propensities must be n x 2^m matrix (n=", nrow(data), ", 2^m=", 2^m, ")")
    }
    if (!is.matrix(oracle_outcomes) ||
        nrow(oracle_outcomes) != nrow(data) ||
        ncol(oracle_outcomes) != 2^m) {
      stop("oracle_outcomes must be n x 2^m matrix (n=", nrow(data), ", 2^m=", 2^m, ")")
    }

    # Validate propensity constraints
    if (any(!is.finite(oracle_propensities))) {
      stop("oracle_propensities contains non-finite values")
    }
    ps_row_sums <- rowSums(oracle_propensities)
    if (!all(abs(ps_row_sums - 1) < PROPENSITY_SUM_TOL)) {
      max_deviation <- max(abs(ps_row_sums - 1))
      stop("oracle_propensities rows must sum to 1 (within tolerance). Max deviation: ",
           sprintf("%.2e", max_deviation))
    }
    if (any(!is.finite(oracle_outcomes))) {
      stop("oracle_outcomes contains non-finite values")
    }

    # Store as estimation results (bypass actual fitting)
    ps_est <- list(propensity_scores = oracle_propensities, models = NULL)
    outcome_est <- list(outcome_means = oracle_outcomes, models = NULL)

  } else if (!is.null(oracle_propensities) || !is.null(oracle_outcomes)) {
    stop("Must provide BOTH oracle_propensities AND oracle_outcomes, or neither")
  } else {
    # Standard mode: set flags to indicate we need to fit models later
    ps_est <- NULL
    outcome_est <- NULL
  }

  # Route to appropriate estimator based on threshold_type
  if (threshold_type == "vanishing") {
    # Vanishing threshold estimator
    return(estimate_qdelta_vanishing_threshold(
      data = data,
      delta = delta,
      treatment_cols = treatment_cols,
      covariate_cols = covariate_cols,
      outcome_col = outcome_col,
      divergence = divergence,
      ps_method = ps_method,
      binary_method = binary_method,
      outcome_method = outcome_method,
      estimator = estimator,
      cross_fit = cross_fit,
      n_folds = n_folds,
      confidence_level = confidence_level,
      bootstrap = bootstrap,
      n_bootstrap = n_bootstrap,
      gamma = gamma,
      min_active = min_active,
      max_qdelta_failure_rate = max_qdelta_failure_rate,
      qdelta_tol = qdelta_tol
    ))
  } else if (threshold_type == "adaptive") {
    # Adaptive estimator with fixed epsilon
    return(estimate_qdelta_adaptive(
      data = data,
      delta = delta,
      treatment_cols = treatment_cols,
      covariate_cols = covariate_cols,
      outcome_col = outcome_col,
      epsilon = epsilon,
      divergence = divergence,
      ps_method = ps_method,
      binary_method = binary_method,
      outcome_method = outcome_method,
      estimator = estimator,
      cross_fit = cross_fit,
      n_folds = n_folds,
      confidence_level = confidence_level,
      bootstrap = bootstrap,
      n_bootstrap = n_bootstrap,
      scaling_method = scaling_method,
      inactive_fraction = inactive_fraction,
      min_active = min_active,
      max_qdelta_failure_rate = max_qdelta_failure_rate,
      qdelta_tol = qdelta_tol
    ))
  }

  # Otherwise, threshold_type == "none", continue with standard estimator below

  # Generate treatment combinations
  vectors <- generate_binary_vectors(m)

  # Estimate nuisance parameters (if not in oracle mode)
  if (is.null(ps_est)) {
    ps_est <- estimate_propensity_scores(data,
                                        treatment_cols = treatment_cols,
                                        covariate_cols = covariate_cols,
                                        method = ps_method,
                                        binary_method = binary_method,
                                        cross_fit = cross_fit, n_folds = n_folds,
                                        xgb_params = xgb_params)
    gc(verbose = FALSE)  # Clear memory after PS estimation
  }

  if (is.null(outcome_est)) {
    outcome_est <- estimate_outcome_models(data,
                                          treatment_cols = treatment_cols,
                                          covariate_cols = covariate_cols,
                                          outcome_col = outcome_col,
                                          method = outcome_method,
                                          cross_fit = cross_fit, n_folds = n_folds,
                                          xgb_params = xgb_params)
    gc(verbose = FALSE)  # Clear memory after outcome estimation
  }

  # Validate propensity scores
  if (any(!is.finite(ps_est$propensity_scores))) {
    stop("Propensity scores contain non-finite values (NA, NaN, or Inf). ",
         "Check your propensity score estimation method and data quality.")
  }

  # Check propensity scores sum to 1 (within tolerance)
  ps_row_sums <- rowSums(ps_est$propensity_scores)
  if (!all(abs(ps_row_sums - 1) < PROPENSITY_SUM_TOL)) {
    max_deviation <- max(abs(ps_row_sums - 1))
    stop("Propensity scores do not sum to 1 (max deviation: ", sprintf("%.2e", max_deviation), "). ",
         "This indicates a bug in the propensity score estimation. ",
         "Check for normalization issues in method '", ps_method, "'.")
  }

  # Check for severe positivity violations (combinations with near-zero probability everywhere)
  min_ps_per_combo <- apply(ps_est$propensity_scores, 2, max)
  severe_violations <- sum(min_ps_per_combo < PROPENSITY_MIN_THRESHOLD)
  if (severe_violations > 0) {
    warning("Severe positivity violation: ", severe_violations, " treatment combinations have max probability < ",
            sprintf("%.2e", PROPENSITY_MIN_THRESHOLD), ". ",
            "Results may be unstable. Consider using fewer treatments or checking treatment overlap.")
  }

  # Validate outcome means
  if (any(!is.finite(outcome_est$outcome_means))) {
    stop("Outcome means contain non-finite values (NA, NaN, or Inf). ",
         "Check your outcome model estimation method and data quality.")
  }

  # Check feasibility at per-observation level
  # Use early stopping only if prob_infeasible == 0 (strict mode)
  use_early_stop <- (prob_infeasible == 0)

  feasibility_analysis <- analyze_delta_feasibility_all_obs(
    ps_est$propensity_scores, delta,
    early_stop = use_early_stop
  )

  # Track infeasible observations for filtering
  feasible_obs <- NULL
  n_infeasible <- nrow(data) - feasibility_analysis$n_feasible
  prop_infeasible <- n_infeasible / nrow(data)

  if (feasibility_analysis$n_feasible < nrow(data)) {
    # Some observations are infeasible

    # Determine whether to error or continue based on prob_infeasible threshold
    if (prop_infeasible >= prob_infeasible) {
      # Too many infeasible observations - error as before

      # Check if early stopped
      if (!is.null(feasibility_analysis$early_stopped) && feasibility_analysis$early_stopped) {
        # Early stop: only checked first few observations
        stop(sprintf(
          "Delta infeasible for observation %d (first infeasible found).\n
Only checked %d / %d observations before stopping early.
All %d preceding observations were feasible.\n
Proportion infeasible: %.3f >= threshold %.3f\n
To find a feasible delta, run:
  application/analysis/00a_check_delta_feasibility.R\n
This script uses binary search to find the maximum feasible delta
without checking all %d observations.\n
Alternatively, use a smaller delta (all components < %.2f) or increase prob_infeasible threshold.",
          feasibility_analysis$first_infeasible,
          feasibility_analysis$first_infeasible,
          nrow(data),
          feasibility_analysis$n_feasible,
          prop_infeasible,
          prob_infeasible,
          nrow(data),
          min(delta) * 0.9  # Suggest reducing delta by 10%
        ), call. = FALSE)
      }

      # Full check (early_stop=FALSE): checked all observations
      # Try to find feasible alternative
      suggested_result <- suggest_feasible_delta(ps_est$propensity_scores, delta)

      if (suggested_result$n_obs_feasible == nrow(data)) {
        stop(sprintf(
          "Delta infeasible for %d out of %d observations (%.1f%% >= threshold %.1f%%).\n
Suggested feasible delta: %s\n
Use this delta, modify your requested delta, or increase prob_infeasible threshold.",
          n_infeasible,
          nrow(data),
          100 * prop_infeasible,
          100 * prob_infeasible,
          paste(round(suggested_result$suggested_delta, 4), collapse = ", ")
        ), call. = FALSE)
      } else {
        stop(sprintf(
          "Delta infeasible for %d out of %d observations (%.1f%% >= threshold %.1f%%).\n
No universally feasible delta found. Consider using different delta values,
checking the support patterns in your propensity scores, or increase prob_infeasible threshold.",
          n_infeasible,
          nrow(data),
          100 * prop_infeasible,
          100 * prob_infeasible
        ), call. = FALSE)
      }
    } else {
      # Proportion of infeasible observations is below threshold - continue with warning
      warning(sprintf(
        "Delta infeasible for %d out of %d observations (%.2f%%).\n
This is below the threshold of %.2f%%, so continuing estimation with only feasible observations.\n
Results will be based on %d observations instead of %d.",
        n_infeasible,
        nrow(data),
        100 * prop_infeasible,
        100 * prob_infeasible,
        feasibility_analysis$n_feasible,
        nrow(data)
      ), call. = FALSE, immediate. = TRUE)

      # Identify feasible observations for filtering
      feasible_obs <- which(feasibility_analysis$feasible_per_obs)
    }
  }

  # Filter to feasible observations if needed
  n_original <- nrow(data)
  if (!is.null(feasible_obs)) {
    # Filter data, propensity scores, and outcome means to feasible observations only
    data <- data[feasible_obs, , drop = FALSE]
    ps_est$propensity_scores <- ps_est$propensity_scores[feasible_obs, , drop = FALSE]
    outcome_est$outcome_means <- outcome_est$outcome_means[feasible_obs, , drop = FALSE]
  }

  # Compute per-observation Q_delta
  qdelta_results <- vector("list", nrow(data))
  convergence_info <- list(
    converged = logical(nrow(data)),
    iterations = integer(nrow(data))
  )

  # One tolerance for every sample size (no relaxation above 100K observations; see qdelta_tol)
  adaptive_tol <- resolve_qdelta_tol(qdelta_tol)
  n_fallback <- 0L
  max_failure_rate <- resolve_max_failure_rate(max_qdelta_failure_rate)

  # Progress reporting for large datasets
  report_interval <- if (nrow(data) > QDELTA_LARGE_N_THRESHOLD) 10000 else Inf
  start_time <- Sys.time()
  timeout_seconds <- QDELTA_DEFAULT_TIMEOUT_SECS  # 5 minutes per observation

  for (i in 1:nrow(data)) {

    # Progress reporting
    if (i %% report_interval == 0) {
      pct <- 100 * i / nrow(data)
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      rate <- i / elapsed
      remaining <- (nrow(data) - i) / rate
      message(sprintf("Q_delta: %d / %d (%.1f%%) | %.1f obs/sec | ETA: %.1f min",
                     i, nrow(data), pct, rate, remaining / 60))
    }

    # Compute with timeout protection
    obs_start_time <- Sys.time()
    qdelta_result_i <- tryCatch({
      result <- compute_qdelta_per_obs(
        ps_est$propensity_scores[i, ], delta, divergence, tol = adaptive_tol
      )

      # Check if this observation took too long (warning only, no interruption)
      elapsed_obs <- as.numeric(difftime(Sys.time(), obs_start_time, units = "secs"))
      if (elapsed_obs > timeout_seconds) {
        warning(sprintf("Observation %d took %.1f seconds (> %d second threshold). ",
                       i, elapsed_obs, timeout_seconds,
                       "Consider investigating this observation."), call. = FALSE)
      }

      result
    }, error = function(e) {
      n_fallback <<- n_fallback + 1L
      qdelta_fallback_or_stop(i, e, ps_est$propensity_scores[i, ], delta,
                              n_fallback, nrow(data), max_failure_rate)
    })

    qdelta_results[[i]] <- qdelta_result_i
    convergence_info$converged[i] <- qdelta_result_i$converged
    convergence_info$iterations[i] <- qdelta_result_i$iterations
  }

  if (n_fallback > 0L) {
    warning(sprintf("Q_delta fallback to P0 used for %d of %d observations (max_qdelta_failure_rate = %g).",
                    n_fallback, nrow(data), max_failure_rate), call. = FALSE)
  }

  # Compute convergence summary statistics
  n_converged <- sum(convergence_info$converged, na.rm = TRUE)
  convergence_rate <- mean(convergence_info$converged, na.rm = TRUE)
  converged_iters <- convergence_info$iterations[convergence_info$converged]
  mean_iterations <- if (length(converged_iters) > 0) mean(converged_iters, na.rm = TRUE) else NA
  max_iterations <- max(convergence_info$iterations, na.rm = TRUE)

  convergence_summary <- list(
    n_converged = n_converged,
    convergence_rate = convergence_rate,
    mean_iterations = mean_iterations,
    max_iterations = max_iterations,
    n_total = nrow(data)
  )

  # Enhanced convergence warning with summary
  n_not_converged <- sum(!convergence_info$converged, na.rm = TRUE)
  if (n_not_converged > 0 || convergence_rate < MIN_CONVERGENCE_RATE) {
    warning(sprintf(
      "Q_delta convergence issues detected:\n  - Converged: %d / %d (%.1f%%)\n  - Mean iterations (converged): %.1f\n  - Max iterations: %d",
      n_converged, nrow(data), 100 * convergence_rate, mean_iterations, max_iterations
    ), call. = FALSE)
  } else {
    message(sprintf("Q_delta solver: %d / %d converged (%.1f%%), mean iterations: %.1f",
                   n_converged, nrow(data), 100 * convergence_rate, mean_iterations))
  }

  gc(verbose = FALSE)  # Clear memory after Q_delta computation

  # Compute point estimate
  if (estimator == "plugin") {
    # DEPRECATED: Plugin estimator is deprecated. Use "one_step" instead.
    warning("The 'plugin' estimator is deprecated and will be removed in a future version. ",
            "Please use 'one_step' instead. ",
            "The plugin estimator does not include the efficient-influence-function ",
            "correction and may be biased if the outcome model is misspecified.")

    psi_est <- compute_psi_delta(data, delta,
                                treatment_cols = treatment_cols,
                                covariate_cols = covariate_cols,
                                outcome_col = outcome_col,
                                qdelta_results = qdelta_results,
                                outcome_means = outcome_est$outcome_means,
                                propensity_scores = ps_est$propensity_scores,
                                vectors = vectors, method = "outcome_regression", divergence = divergence)

    # For plugin estimator, we need to compute EIF for variance estimation
    # Extract X, A, Y from data (explicit column specifications)
    X_matrix <- if (length(covariate_cols) > 0) {
      as.matrix(data[, covariate_cols, drop = FALSE])
    } else {
      NULL
    }
    A_matrix <- as.matrix(data[, treatment_cols, drop = FALSE])
    Y_vector <- data[[outcome_col]]

    eif_result <- efficient_influence_function(X = X_matrix, A = A_matrix, Y = Y_vector,
                                              propensity_scores = ps_est$propensity_scores,
                                              delta = delta,
                                              outcome_means = outcome_est$outcome_means,
                                              divergence = divergence,
                                              qdelta_results = qdelta_results)
    eif_values <- eif_result$centered_eif
    gc(verbose = FALSE)  # Clear memory after EIF computation

  } else if (estimator == "one_step") {
    # One-step estimator: mean of EIF values
    # Formula: ψ_one-step = (1/n) Σᵢ EIFᵢ

    # Extract X, A, Y from data (explicit column specifications)
    X_matrix <- if (length(covariate_cols) > 0) {
      as.matrix(data[, covariate_cols, drop = FALSE])
    } else {
      NULL
    }
    A_matrix <- as.matrix(data[, treatment_cols, drop = FALSE])
    Y_vector <- data[[outcome_col]]

    # Compute EIF with pre-computed qdelta_results
    eif_result <- efficient_influence_function(
      X = X_matrix,
      A = A_matrix,
      Y = Y_vector,
      propensity_scores = ps_est$propensity_scores,
      delta = delta,
      outcome_means = outcome_est$outcome_means,
      divergence = divergence,
      qdelta_results = qdelta_results
    )

    # One-step estimate is the mean of EIF values
    psi_est <- mean(eif_result$eif_values)
    eif_values <- eif_result$centered_eif
    gc(verbose = FALSE)  # Clear memory after EIF computation

  } else if (estimator == "tmle") {
    warning("The 'tmle' estimator is EXPERIMENTAL and is not covered by the theory ",
            "in the accompanying paper. No asymptotic normality or efficiency result ",
            "is established for it, and its targeting step fluctuates only the outcome ",
            "model while the estimand also depends on the propensity score through the ",
            "divergence tilt. Use estimator = 'one_step' for inference you intend to report.",
            call. = FALSE)
    # TMLE estimator with full targeting step (experimental, no paper-backed theory)
    # Algorithm:
    # 1. Get initial outcome model μₐ(X)
    # 2. Perform targeting step to get μ*ₐ(X)
    # 3. Compute EIF with targeted outcome means
    # 4. Final estimate: mean of EIF values


    # Step 1: Initial outcome model (already have outcome_est$outcome_means)
    initial_outcome_means <- outcome_est$outcome_means

    # Step 2: Perform targeting step
    targeted_outcome_means <- target_outcome_model_tmle(
      data = data,
      outcome_means = initial_outcome_means,
      propensity_scores = ps_est$propensity_scores,
      qdelta_results = qdelta_results,
      vectors = vectors,
      treatment_cols = treatment_cols,
      outcome_col = outcome_col,
      divergence = divergence
    )

    # Step 3: Extract X, A, Y from data (explicit column specifications)
    X_matrix <- if (length(covariate_cols) > 0) {
      as.matrix(data[, covariate_cols, drop = FALSE])
    } else {
      NULL
    }
    A_matrix <- as.matrix(data[, treatment_cols, drop = FALSE])
    Y_vector <- data[[outcome_col]]

    # Step 4: Compute EIF with targeted outcome means
    eif_result <- efficient_influence_function(
      X = X_matrix,
      A = A_matrix,
      Y = Y_vector,
      propensity_scores = ps_est$propensity_scores,
      delta = delta,
      outcome_means = targeted_outcome_means,
      divergence = divergence,
      qdelta_results = qdelta_results
    )

    # Step 5: Final estimate is mean of EIF values
    psi_est <- mean(eif_result$eif_values)
    eif_values <- eif_result$centered_eif
    gc(verbose = FALSE)  # Clear memory after EIF computation

    # Store targeted outcome means in nuisance estimates
    outcome_est$outcome_means <- targeted_outcome_means
    outcome_est$targeted <- TRUE

  } else {
    stop("estimator must be one of: 'one_step' (recommended), 'tmle' (experimental, ",
         "no paper-backed theory), or 'plugin' (deprecated).")
  }

  # Apply contrast modification (difference from observed mean)
  if (contrast == "difference") {
    # Compute observed mean
    observed_mean <- mean(Y_vector)

    # Modify estimate: E[Y^{Q_delta}] - E[Y]
    psi_est <- psi_est - observed_mean

    # Modify EIF: EIF_{E[Y^{Q_delta}] - E[Y]} = EIF_{E[Y^{Q_delta}]} - (Y_i - E[Y])
    # This gives correct inference for the difference parameter
    eif_values <- eif_values - (Y_vector - observed_mean)
  }

  # Validate EIF values
  if (any(!is.finite(eif_values))) {
    n_nonfinite <- sum(!is.finite(eif_values))
    stop("EIF contains ", n_nonfinite, " non-finite values (NA, NaN, or Inf). ",
         "This indicates numerical instability. ",
         "Try: (1) increasing regularization, (2) checking for positivity violations, ",
         "or (3) using a different divergence measure.")
  }

  # Check if EIF is approximately centered (mean ≈ 0)
  eif_mean <- mean(eif_values)
  eif_sd <- sd(eif_values)
  if (abs(eif_mean) > EIF_CENTERING_THRESHOLD * eif_sd && eif_sd > NEAR_ZERO_THRESHOLD) {
    warning("EIF may not be properly centered (mean = ", sprintf("%.4f", eif_mean),
            ", sd = ", sprintf("%.4f", eif_sd), "). ",
            "This could indicate estimation issues. ",
            "EIF should have mean close to 0 for valid inference.")
  }

  # Compute standard error and confidence interval
  if (bootstrap) {
    # Bootstrap inference
    bootstrap_estimates <- numeric(n_bootstrap)

    for (b in 1:n_bootstrap) {
      # Bootstrap sample
      boot_idx <- sample(1:nrow(data), nrow(data), replace = TRUE)
      boot_data <- data[boot_idx, ]

      # Re-estimate nuisance parameters for bootstrap sample
      boot_ps_est <- estimate_propensity_scores(boot_data,
                                               treatment_cols = treatment_cols,
                                               covariate_cols = covariate_cols,
                                               method = ps_method)
      boot_outcome_est <- estimate_outcome_models(boot_data,
                                                 treatment_cols = treatment_cols,
                                                 covariate_cols = covariate_cols,
                                                 outcome_col = outcome_col,
                                                 method = outcome_method)

      # Compute per-observation Q_delta for bootstrap sample
      boot_qdelta_results <- vector("list", nrow(boot_data))
      for (i in 1:nrow(boot_data)) {
        boot_qdelta_results[[i]] <- compute_qdelta_per_obs(
          boot_ps_est$propensity_scores[i, ], delta, divergence, tol = 1e-8
        )
      }

      bootstrap_estimates[b] <- compute_psi_delta(boot_data, delta,
                                                 treatment_cols = treatment_cols,
                                                 covariate_cols = covariate_cols,
                                                 outcome_col = outcome_col,
                                                 qdelta_results = boot_qdelta_results,
                                                 outcome_means = boot_outcome_est$outcome_means,
                                                 propensity_scores = boot_ps_est$propensity_scores,
                                                 m = m, method = "doubly_robust", divergence = divergence)
    }

    se <- sd(bootstrap_estimates)
    alpha <- 1 - confidence_level
    ci_lower <- quantile(bootstrap_estimates, alpha/2)
    ci_upper <- quantile(bootstrap_estimates, 1 - alpha/2)

  } else {
    # Asymptotic inference using EIF
    # Check for non-finite values in EIF before computing SE
    if (any(!is.finite(eif_values))) {
      n_bad <- sum(!is.finite(eif_values))
      warning(sprintf(
        "Found %d non-finite EIF values (NA/NaN/Inf) out of %d observations. ",
        n_bad, length(eif_values),
        "These will be excluded from SE computation. ",
        "This may indicate numerical issues in Q_delta solver or EIF computation."
      ), call. = FALSE)
      eif_values <- eif_values[is.finite(eif_values)]
      if (length(eif_values) == 0) {
        stop("All EIF values are non-finite. Cannot compute standard error.",
             call. = FALSE)
      }
    }

    se <- sqrt(var(eif_values) / nrow(data))
    alpha <- 1 - confidence_level
    z_alpha <- qnorm(1 - alpha/2)
    ci_lower <- psi_est - z_alpha * se
    ci_upper <- psi_est + z_alpha * se
  }

  # Compute p-value for testing H0: psi_delta = 0
  z_stat <- psi_est / se
  p_value <- 2 * (1 - pnorm(abs(z_stat)))

  return(list(
    estimate = psi_est,
    std_error = se,
    se = se,  # Backward compatibility alias for std_error
    variance = se^2,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    confidence_level = confidence_level,
    p_value = p_value,
    converged = convergence_summary$convergence_rate >= 0.95,  # Top-level flag for simulation compatibility
    eif_values = eif_values,
    centered_eif = eif_values,
    qdelta_results = qdelta_results,
    convergence_info = convergence_info,
    convergence_summary = convergence_summary,
    nuisance_estimates = list(
      propensity_scores = ps_est,
      outcome_models = outcome_est
    ),
    delta = delta,
    divergence = divergence,
    estimator = estimator,
    contrast = contrast,
    n = nrow(data),
    m = m,
    # Infeasibility filtering information
    n_original = n_original,
    n_infeasible_excluded = if (!is.null(feasible_obs)) n_original - nrow(data) else 0,
    prop_infeasible_excluded = if (!is.null(feasible_obs)) (n_original - nrow(data)) / n_original else 0,
    prob_infeasible_threshold = prob_infeasible
  ))
}

#' Identify eligibility patterns in data
#'
#' @description
#' Groups observations by unique eligibility patterns across treatments.
#' Each pattern is defined by a binary vector indicating which treatments
#' are eligible for that group of observations.
#'
#' @param data Data frame with treatment and eligibility columns
#' @param eligibility_cols Character vector of eligibility column names.
#'   If NULL, auto-detects columns named "elig_A1", "elig_A2", etc.
#'
#' @return List containing:
#' \describe{
#'   \item{patterns}{Data frame with:}
#'     \describe{
#'       \item{pattern_id}{Integer, sequential pattern identifier}
#'       \item{pattern_string}{Character, human-readable pattern (e.g., "10" for A1 eligible, A2 not)}
#'       \item{eligibility}{List column, binary vector indicating eligible treatments}
#'       \item{n_obs}{Integer, number of observations with this pattern}
#'       \item{prevalence}{Numeric, proportion of sample with this pattern}
#'       \item{n_eligible}{Integer, number of eligible treatments in this pattern}
#'       \item{obs_indices}{List column, indices of observations with this pattern}
#'     }
#'   \item{pattern_assignments}{Integer vector (length n), pattern ID for each observation}
#'   \item{m}{Integer, total number of treatments}
#' }
#'
#' @keywords internal
identify_eligibility_patterns <- function(data, eligibility_cols) {

  # Eligibility columns must be provided explicitly (no auto-detection).
  if (is.null(eligibility_cols) || length(eligibility_cols) == 0) {
    stop("`eligibility_cols` must be provided (e.g. c('elig_A1', 'elig_A2', ...)).",
         call. = FALSE)
  }

  # Validate eligibility columns exist
  missing_cols <- setdiff(eligibility_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Eligibility columns not found in data: ", paste(missing_cols, collapse = ", "))
  }

  m <- length(eligibility_cols)
  n <- nrow(data)

  # Extract eligibility matrix (n × m)
  elig_matrix <- as.matrix(data[, eligibility_cols, drop = FALSE])

  # Validate eligibility values are binary
  if (!all(elig_matrix %in% c(0, 1))) {
    stop("Eligibility columns must contain only 0 and 1")
  }

  # Create pattern string for each observation
  # E.g., "101" means A1 eligible, A2 not, A3 eligible
  pattern_strings <- apply(elig_matrix, 1, function(row) paste(row, collapse = ""))

  # Get unique patterns
  unique_patterns <- unique(pattern_strings)
  n_patterns <- length(unique_patterns)

  # Build patterns data frame
  patterns_list <- vector("list", n_patterns)

  for (i in seq_along(unique_patterns)) {
    pattern_str <- unique_patterns[i]

    # Find observations with this pattern
    obs_indices <- which(pattern_strings == pattern_str)
    n_obs_pattern <- length(obs_indices)

    # Extract eligibility vector for this pattern
    elig_vec <- as.integer(strsplit(pattern_str, "")[[1]])
    n_eligible <- sum(elig_vec)

    patterns_list[[i]] <- list(
      pattern_id = i,
      pattern_string = pattern_str,
      eligibility = elig_vec,
      n_obs = n_obs_pattern,
      prevalence = n_obs_pattern / n,
      n_eligible = n_eligible,
      obs_indices = obs_indices
    )
  }

  # Convert to data frame (preserving list columns)
  patterns_df <- do.call(rbind, lapply(patterns_list, function(x) {
    data.frame(
      pattern_id = x$pattern_id,
      pattern_string = x$pattern_string,
      n_obs = x$n_obs,
      prevalence = x$prevalence,
      n_eligible = x$n_eligible,
      stringsAsFactors = FALSE
    )
  }))

  # Add list columns separately (data.frame drops them)
  patterns_df$eligibility <- lapply(patterns_list, function(x) x$eligibility)
  patterns_df$obs_indices <- lapply(patterns_list, function(x) x$obs_indices)

  # Sort by prevalence (descending)
  patterns_df <- patterns_df[order(-patterns_df$prevalence), ]

  # Reassign pattern IDs after sorting
  patterns_df$pattern_id <- seq_len(nrow(patterns_df))

  # Create pattern assignment vector for all observations
  pattern_assignments <- integer(n)
  for (i in seq_len(nrow(patterns_df))) {
    obs_idx <- patterns_df$obs_indices[[i]]
    pattern_assignments[obs_idx] <- patterns_df$pattern_id[i]
  }

  return(list(
    patterns = patterns_df,
    pattern_assignments = pattern_assignments,
    m = m
  ))
}

#' Validate feasibility of stratification approach
#'
#' @description
#' Checks whether eligibility-stratified estimation is feasible for the
#' given data. Validates that there are enough patterns and observations
#' to support stratified estimation.
#'
#' @param patterns Data frame from identify_eligibility_patterns()
#' @param min_stratum_size Integer, minimum observations per pattern to consider feasible
#' @param min_coverage Numeric, minimum proportion of sample that must be in feasible patterns
#'
#' @return List containing:
#' \describe{
#'   \item{feasible}{Logical, whether stratification is feasible}
#'   \item{reason}{Character, reason for infeasibility (or "Feasible" if feasible)}
#'   \item{n_feasible_patterns}{Integer, number of patterns with n >= min_stratum_size}
#'   \item{coverage}{Numeric, proportion of sample in feasible patterns}
#'   \item{feasible_pattern_ids}{Integer vector, IDs of feasible patterns}
#' }
#'
#' @keywords internal
validate_stratification_feasibility <- function(patterns, min_stratum_size = 100, min_coverage = 0.5) {

  # Check: At least 2 patterns total
  if (nrow(patterns) < 2) {
    return(list(
      feasible = FALSE,
      reason = paste0("Only ", nrow(patterns), " eligibility pattern(s) found. Need at least 2 for stratification."),
      n_feasible_patterns = nrow(patterns),
      coverage = sum(patterns$prevalence),
      feasible_pattern_ids = patterns$pattern_id
    ))
  }

  # Check: Each pattern has at least one eligible treatment
  patterns_no_eligible <- patterns[patterns$n_eligible == 0, ]
  if (nrow(patterns_no_eligible) > 0) {
    return(list(
      feasible = FALSE,
      reason = paste0(nrow(patterns_no_eligible), " pattern(s) have no eligible treatments."),
      n_feasible_patterns = NA,
      coverage = NA,
      feasible_pattern_ids = integer(0)
    ))
  }

  # Filter to patterns with sufficient sample size
  feasible_patterns <- patterns[patterns$n_obs >= min_stratum_size, ]
  n_feasible <- nrow(feasible_patterns)

  # Check: At least 2 patterns with sufficient sample size
  if (n_feasible < 2) {
    return(list(
      feasible = FALSE,
      reason = paste0("Only ", n_feasible, " pattern(s) with n >= ", min_stratum_size, ". Need at least 2."),
      n_feasible_patterns = n_feasible,
      coverage = sum(feasible_patterns$prevalence),
      feasible_pattern_ids = feasible_patterns$pattern_id
    ))
  }

  # Check: Feasible patterns cover sufficient proportion of sample
  coverage <- sum(feasible_patterns$prevalence)
  if (coverage < min_coverage) {
    return(list(
      feasible = FALSE,
      reason = sprintf("Feasible patterns cover only %.1f%% of sample (< %.1f%% threshold).",
                      100 * coverage, 100 * min_coverage),
      n_feasible_patterns = n_feasible,
      coverage = coverage,
      feasible_pattern_ids = feasible_patterns$pattern_id
    ))
  }

  # All checks passed
  return(list(
    feasible = TRUE,
    reason = "Feasible",
    n_feasible_patterns = n_feasible,
    coverage = coverage,
    feasible_pattern_ids = feasible_patterns$pattern_id
  ))
}

#' Estimate effect within a single eligibility pattern
#'
#' @description
#' Estimates ψ_k(δ) within a single eligibility pattern by subsetting data,
#' renaming treatments to consecutive indices, and calling estimate_qdelta_mean().
#'
#' @param data Data frame, full dataset
#' @param pattern_info List, single row from patterns data frame (converted to list)
#' @param delta Numeric vector (length m), target deltas for all treatments
#' @param covariate_cols Character vector of covariate column names (required).
#' @param outcome_col Character scalar, outcome column name. Default "Y".
#' @param ps_method Character, propensity score method
#' @param outcome_method Character, outcome model method
#' @param estimator Character, estimator type
#' @param cross_fit Logical, use cross-fitting
#' @param n_folds Integer, number of folds
#' @param divergence Character, divergence measure
#' @param threshold_type Character, active-set selection for handling near-boundary
#'   marginal propensities within this pattern. One of "none" (default; standard
#'   estimator, no thresholding), "vanishing" (vanishing threshold ε_n → 0), or
#'   "adaptive" (fixed threshold ε with adaptive delta). Forwarded to
#'   \code{\link{estimate_qdelta_mean}}. Default "none" preserves prior behavior.
#' @param gamma Numeric, threshold rate for the vanishing estimator (`ε_n = C n^{-γ}`).
#'   Only used when threshold_type = "vanishing". Default NULL (auto-selected).
#' @param epsilon Numeric, fixed threshold for the adaptive estimator. Only used
#'   when threshold_type = "adaptive". Default 0.1.
#' @param verbose Logical, print progress
#' @param xgb_params Optional list of xgboost parameters, forwarded to the
#'   nuisance estimators.
#' @param nuisance_cache Optional \code{nuisance_cache} (see
#'   \code{\link{new_nuisance_cache}}). When supplied, this pattern's
#'   cross-fitted propensity and outcome matrices are taken from the cache if
#'   present and otherwise fitted once (with namespaced per-fit seeds) and
#'   stored. The matrices are handed to \code{\link{estimate_qdelta_mean}}
#'   through its existing pre-supplied-nuisance bypass
#'   (\code{oracle_propensities} / \code{oracle_outcomes}), so no new bypass is
#'   introduced. Only supported with \code{threshold_type = "none"}; the
#'   thresholded estimators fit their own nuisances and ignore pre-supplied
#'   ones, so passing a cache with them is an error rather than a silent no-op.
#'   Default NULL (no caching; behaviour identical to previous versions).
#' @param nuisance_tag Character scalar identifying the outcome (cost type) for
#'   cache keying. Outcome fits are cached per tag; propensity fits are shared
#'   across tags. Ignored when \code{nuisance_cache} is NULL.
#' @param pattern_pos Integer, 1-based position of this pattern in the caller's
#'   pattern loop. Required when \code{nuisance_cache} is supplied: it
#'   namespaces the per-fit seeds so they do not depend on how much RNG earlier
#'   patterns consumed. Ignored otherwise.
#'
#' @param max_qdelta_failure_rate Numeric in `[0, 1]`: the largest share of observations whose
#'   Q_delta solve may fail and fall back to Q_delta = P0 (no shift for that observation, a bias
#'   toward the null). Default `0`: the first failure stops the estimate with the observation index
#'   and the solver's message. A positive value tolerates failures up to `rate * n` observations
#'   (each warned about and counted), then stops. `1` always falls back. The vanishing-threshold
#'   estimator never falls back and always stops.
#' @param qdelta_tol Numeric or `NULL`. Solver tolerance for the per-observation Q_delta problem;
#'   `NULL` (default) uses `QDELTA_TOLERANCE` for every sample size (it no longer relaxes above
#'   100,000 observations, which made the same stratum's estimate depend on whether it crossed that
#'   size).
#' @return List containing:
#' \describe{
#'   \item{success}{Logical, whether estimation succeeded}
#'   \item{estimate}{Numeric, estimated ψ_k (or NA if failed)}
#'   \item{se}{Numeric, standard error (or NA if failed)}
#'   \item{eif_values_extended}{Numeric vector (length n), EIF extended to full dataset with zeros}
#'   \item{n_obs}{Integer, number of observations in pattern}
#'   \item{n_eligible}{Integer, number of eligible treatments}
#'   \item{treatment_mapping}{List with original and reduced indices}
#'   \item{error_message}{Character, error message if failed (or NA if success)}
#' }
#'
#' @section Relationship to Active Set Functions:
#' This function performs eligibility-based treatment subsetting at the data.frame
#' level (observation-specific eligibility), while the active set functions in
#' \code{\link{efficient_influence_function_active}} perform propensity-based
#' thresholding at the matrix level (combination-specific active sets).
#'
#' Both approaches use shared utilities from treatment_subset_utils.R:
#' \itemize{
#'   \item \code{\link{validate_active_set}} - validate treatment indices
#'   \item \code{\link{map_treatment_indices}} - track original→reduced mapping
#'   \item \code{\link{extend_eif_to_full}} - extend subset EIF to full dataset
#' }
#'
#' The treatment renaming (lines 1205-1219) operates on data.frames and is separate
#' from the matrix-level operations in \code{\link{compute_omega_reduced}}, which
#' are used internally during EIF computation.
#'
#' @keywords internal
estimate_pattern_effect <- function(data, pattern_info, delta,
                                    covariate_cols, outcome_col = "Y",
                                    ps_method = "glm", binary_method = "glm",
                                    outcome_method = "glm",
                                    estimator = "one_step", cross_fit = TRUE, n_folds = 5,
                                    divergence = "kl",
                                    threshold_type = c("none", "vanishing", "adaptive"),
                                    gamma = NULL, epsilon = 0.1,
                                    verbose = TRUE, xgb_params = NULL,
                                    nuisance_cache = NULL, nuisance_tag = NULL,
                                    pattern_pos = NULL,
    max_qdelta_failure_rate = 0, qdelta_tol = NULL) {

  threshold_type <- match.arg(threshold_type)
  n_full <- nrow(data)
  m <- length(delta)

  # Nuisance-fit cache (optional; see R/nuisance_cache.R). Validate up front so a
  # misuse fails before any fitting work. The thresholded estimators delegate to
  # estimate_qdelta_vanishing_threshold()/estimate_qdelta_adaptive(), which fit
  # their own nuisances and ignore pre-supplied ones -- a cache passed with them
  # would be silently useless, so refuse instead.
  use_nuisance_cache <- !is.null(nuisance_cache)
  if (use_nuisance_cache) {
    .check_nuisance_cache(nuisance_cache)
    if (threshold_type != "none") {
      stop("`nuisance_cache` is only supported with threshold_type = 'none'. ",
           "threshold_type = '", threshold_type, "' routes to a thresholded ",
           "estimator that fits its own nuisances and ignores pre-supplied ones, ",
           "so the cache would be silently ignored. Drop the cache or use ",
           "threshold_type = 'none'.", call. = FALSE)
    }
    if (is.null(pattern_pos)) {
      stop("`pattern_pos` must be supplied when `nuisance_cache` is used: it ",
           "namespaces the per-fit seeds so cached and freshly-fitted runs agree.",
           call. = FALSE)
    }
  }

  # Extract pattern information
  obs_indices <- pattern_info$obs_indices
  eligibility <- pattern_info$eligibility
  pattern_id <- pattern_info$pattern_id
  pattern_string <- pattern_info$pattern_string

  # Unlist eligibility if it's a list (from data frame list column)
  if (is.list(eligibility)) {
    eligibility <- unlist(eligibility)
  }
  if (is.list(obs_indices)) {
    obs_indices <- unlist(obs_indices)
  }

  # Identify eligible treatments (original indices)
  eligible_original <- which(eligibility == 1)
  m_k <- length(eligible_original)

  # Validate active set (using shared utility)
  validation <- validate_active_set(eligible_original, m = m, min_active = 0)

  if (!validation$valid) {
    return(list(
      success = FALSE,
      estimate = NA,
      se = NA,
      eif_values_extended = rep(NA, n_full),
      n_obs = length(obs_indices),
      n_eligible = m_k,
      treatment_mapping = list(original = integer(0), reduced = integer(0)),
      error_message = paste("Invalid active set:", paste(validation$issues, collapse = "; "))
    ))
  }

  # Check: Need at least one eligible treatment
  if (m_k == 0) {
    return(list(
      success = FALSE,
      estimate = NA,
      se = NA,
      eif_values_extended = rep(NA, n_full),
      n_obs = length(obs_indices),
      n_eligible = 0,
      treatment_mapping = list(original = integer(0), reduced = integer(0)),
      error_message = "No eligible treatments in pattern"
    ))
  }

  # Create treatment index mapping (using shared utility)
  treatment_mapping_df <- map_treatment_indices(eligible_original, m = m)
  treatment_mapping <- list(
    original = treatment_mapping_df$original_index,
    reduced = treatment_mapping_df$reduced_index
  )

  # Subset data to observations in this pattern
  data_k <- data[obs_indices, , drop = FALSE]

  # Extract eligible treatment columns and rename to A1, A2, ..., A_m_k
  original_treatment_cols <- paste0("A", eligible_original)

  # Treatment subsetting and renaming (data.frame level)
  # NOTE: This prepares data for estimate_qdelta_mean(), which expects consecutive
  # treatment names A1, A2, ..., Am_k. Matrix-level operations (e.g., compute_omega_reduced)
  # happen inside the EIF computation and use the utilities from treatment_subset_utils.R.

  # First, remove ALL original treatment columns (A1, A2, ..., Am)
  all_original_treatment_cols <- paste0("A", 1:m)
  data_k_renamed <- data_k[, !(names(data_k) %in% all_original_treatment_cols), drop = FALSE]

  # Then add back only the eligible treatments with consecutive renaming
  for (i in seq_along(eligible_original)) {
    original_col <- paste0("A", eligible_original[i])
    new_col <- paste0("A", i)
    data_k_renamed[[new_col]] <- data_k[[original_col]]
  }

  # Reorder columns: treatments first, then others
  treatment_cols_new <- paste0("A", 1:m_k)
  other_cols <- setdiff(names(data_k_renamed), treatment_cols_new)
  data_k_renamed <- data_k_renamed[, c(treatment_cols_new, other_cols)]

  # Extract corresponding delta values for eligible treatments
  delta_k <- delta[eligible_original]

  # Note: m_k = 1 is valid (scalar incremental effect on one treatment).
  # No minimum m requirement — estimate_qdelta_mean handles m >= 1.

  # Check: Sufficient sample size for estimation
  if (nrow(data_k_renamed) < max(50, n_folds * 10)) {
    return(list(
      success = FALSE,
      estimate = NA,
      se = NA,
      eif_values_extended = rep(NA, n_full),
      n_obs = nrow(data_k_renamed),
      n_eligible = m_k,
      treatment_mapping = treatment_mapping,
      error_message = paste0("Insufficient sample size (n=", nrow(data_k_renamed), ")")
    ))
  }

  # Check: Disable cross-fitting if too few observations per fold
  cross_fit_k <- cross_fit
  if (cross_fit && nrow(data_k_renamed) / n_folds < 50) {
    if (verbose) {
      message(sprintf("Pattern %d: Disabling cross-fitting (n=%d too small for %d folds)",
                     pattern_id, nrow(data_k_renamed), n_folds))
    }
    cross_fit_k <- FALSE
  }

  # Estimate within this pattern.
  # NOTE: the nuisance fetch/fit lives INSIDE this tryCatch so that a nuisance
  # failure is reported as a pattern failure (as before) rather than aborting the
  # whole scenario. Because the cache is an environment, any fit completed before
  # the error is still retained for later patterns/scenarios.
  result_k <- tryCatch({
    cached_nuisances <- if (use_nuisance_cache) {
      .pattern_nuisances(
        cache = nuisance_cache,
        tag = nuisance_tag,
        pattern_string = as.character(pattern_string),
        pattern_pos = pattern_pos,
        obs_indices = obs_indices,
        data_k = data_k_renamed,
        treatment_cols = treatment_cols_new,
        covariate_cols = covariate_cols,
        outcome_col = outcome_col,
        ps_method = ps_method,
        binary_method = binary_method,
        outcome_method = outcome_method,
        cross_fit = cross_fit_k,
        n_folds = n_folds,
        xgb_params = xgb_params
      )
    } else {
      list(propensity_scores = NULL, outcome_means = NULL)
    }

    fit_pattern <- function() {
      estimate_qdelta_mean(
        data = data_k_renamed,
        delta = delta_k,
        treatment_cols = treatment_cols_new,
        covariate_cols = covariate_cols,
        outcome_col = outcome_col,
        divergence = divergence,
        ps_method = ps_method,
        binary_method = binary_method,
        outcome_method = outcome_method,
        oracle_propensities = cached_nuisances$propensity_scores,
        oracle_outcomes = cached_nuisances$outcome_means,
        estimator = estimator,
        cross_fit = cross_fit_k,
        n_folds = n_folds,
        threshold_type = threshold_type,
        gamma = gamma,
        epsilon = epsilon,
        confidence_level = 0.95,
        bootstrap = FALSE,
        xgb_params = xgb_params,
        max_qdelta_failure_rate = max_qdelta_failure_rate,
        qdelta_tol = qdelta_tol
      )
    }

    if (use_nuisance_cache) {
      # Muffle only estimate_qdelta_mean's "oracle nuisances" notice: with the
      # cache it fires once per pattern per scenario (hundreds of lines in the
      # application log) and says nothing the cache diagnostics don't.
      withCallingHandlers(fit_pattern(), message = function(msg) {
        if (grepl("Using oracle nuisances", conditionMessage(msg), fixed = TRUE)) {
          invokeRestart("muffleMessage")
        }
      })
    } else {
      fit_pattern()
    }
  }, error = function(e) {
    # A cache misalignment is a programming/data error, not a pattern-level
    # estimation failure: rethrow so it aborts loudly instead of being recorded
    # as "this pattern failed" and averaged away.
    if (inherits(e, "nuisance_cache_mismatch")) stop(e)
    # A Q_delta solver failure beyond max_qdelta_failure_rate must stop the run, not remove
    # this pattern from the estimand (the stratified estimator would renormalise the weights).
    if (inherits(e, "qdelta_solver_failure")) stop(e)
    return(list(error = TRUE, message = as.character(e$message)))
  })

  # Check if estimation failed or result is invalid
  if (is.null(result_k) || (!is.null(result_k$error) && result_k$error)) {
    error_msg <- if (is.null(result_k)) {
      "estimate_qdelta_mean returned NULL"
    } else {
      as.character(result_k$message)
    }
    return(list(
      success = FALSE,
      estimate = NA,
      se = NA,
      eif_values_extended = rep(NA, n_full),
      n_obs = nrow(data_k_renamed),
      n_eligible = m_k,
      treatment_mapping = treatment_mapping,
      error_message = error_msg
    ))
  }

  # Validate result structure
  if (is.null(result_k$eif_values) || length(result_k$eif_values) == 0) {
    return(list(
      success = FALSE,
      estimate = NA,
      se = NA,
      eif_values_extended = rep(NA, n_full),
      n_obs = nrow(data_k_renamed),
      n_eligible = m_k,
      treatment_mapping = treatment_mapping,
      error_message = paste0("Invalid EIF values returned (length: ",
                            length(result_k$eif_values), ")")
    ))
  }

  # Extract EIF values for this pattern (length n_k)
  eif_values_k <- result_k$eif_values

  # Extend EIF to full dataset using shared utility
  eif_values_extended <- tryCatch({
    extend_eif_to_full(
      eif_reduced = eif_values_k,
      indices_in_full = obs_indices,
      n_full = n_full
    )
  }, error = function(e) {
    # If extension fails, return error details
    return(list(
      success = FALSE,
      estimate = NA,
      se = NA,
      eif_values_extended = rep(NA, n_full),
      n_obs = nrow(data_k_renamed),
      n_eligible = m_k,
      treatment_mapping = treatment_mapping,
      error_message = paste0("EIF extension failed: ", e$message,
                            " (eif_k length: ", length(eif_values_k),
                            ", obs_indices length: ", length(obs_indices),
                            ", n_full: ", n_full, ")")
    ))
  })

  # Check if extension failed (returned error list instead of vector)
  if (is.list(eif_values_extended) && !is.null(eif_values_extended$success)) {
    return(eif_values_extended)
  }

  # Return success
  return(list(
    success = TRUE,
    estimate = result_k$estimate,
    se = result_k$std_error,
    eif_values_extended = eif_values_extended,
    obs_indices = obs_indices,
    n_obs = nrow(data_k_renamed),
    n_eligible = m_k,
    treatment_mapping = treatment_mapping,
    convergence_summary = result_k$convergence_summary,
    error_message = NA
  ))
}

#' Combine per-pattern estimates into a stratified aggregate (core math)
#'
#' Single source of truth for the stratified aggregation, shared by
#' \code{aggregate_stratified_results()} (used inside estimation) and
#' \code{recombine_stratified_patterns()} (post-hoc, no re-estimation). Because
#' eligibility patterns partition the sample (each observation belongs to exactly
#' one pattern), the aggregate is a prevalence-weighted average of per-pattern
#' estimates. With prevalence weights \eqn{w_k = n_k / n_S}, the pooled estimator
#' is the plain mean of the observation-level influence values over the included
#' patterns, so its influence function for observation \eqn{i} in pattern
#' \eqn{k} is \eqn{(\varphi_i - \hat\psi_k)/c + (\hat\psi_k - \hat\psi)/c},
#' where \eqn{c} is the included share of the sample. The first term is the
#' within-pattern influence function rescaled by \eqn{w_k n / n_k} (it is NOT
#' \eqn{w_k \varphi_i}); the second is the between-pattern term, present because
#' pattern membership and hence the weights are random. For fixed user-supplied
#' weights only the rescaled within-pattern term applies.
#'
#' @param estimates Numeric vector, per-pattern point estimates (included only).
#' @param prevalences Numeric vector, per-pattern prevalence (fraction of the FULL
#'   stratified sample), aligned with \code{estimates}. Renormalized internally.
#' @param obs_eif Numeric vector length \code{n}, the UNWEIGHTED influence value
#'   for each observation from its own pattern's fit (0 or NA where the obs's
#'   pattern is not included). It is re-centered within pattern here, so centered
#'   or uncentered input are both accepted.
#' @param obs_pattern_pos Integer vector length \code{n}, index into
#'   \code{estimates}/\code{prevalences} giving each observation's pattern
#'   position, or NA if the observation's pattern is not among the included set.
#' @param n_full Integer, number of observations in the full stratified sample
#'   (denominator for the variance).
#' @param confidence_level Numeric, for the CI. Default 0.95.
#' @param prevalence_weights Logical. \code{TRUE} (default) when \code{prevalences}
#'   are the estimated pattern prevalences, so the between-pattern term is added.
#'   \code{FALSE} for fixed user-supplied weights.
#' @param coverage Numeric, share of the full sample in the included patterns.
#'   Defaults to \code{sum(prevalences)}; supply it explicitly when
#'   \code{prevalences} are custom weights.
#' @return List: estimate, se, ci_lower, ci_upper, p_value, eif_values (centered,
#'   length n_full), weights (aligned with estimates), coverage (sum of included
#'   prevalence).
#' @keywords internal
.combine_patterns <- function(estimates, prevalences, obs_eif, obs_pattern_pos,
                              n_full, confidence_level = 0.95,
                              prevalence_weights = TRUE, coverage = NULL) {
  # Renormalize prevalences over the INCLUDED patterns so weights sum to 1.
  total_prevalence <- sum(prevalences)
  if (is.null(coverage)) coverage <- total_prevalence
  if (total_prevalence <= 0) {
    stop("Cannot combine patterns: total prevalence of included patterns is 0.")
  }
  weights <- prevalences / total_prevalence

  # Point estimate: psi = sum_k w_k psi_k
  psi_agg <- sum(weights * estimates)

  # Aggregated influence function. For psi = sum_k w_k psi_k with w_k estimated
  # as n_k / n_S, psi is the plain mean of phi_i over the included patterns, so
  #   IF_i = (w_k n / n_k) (phi_i - psi_k)  [within, rescaled by 1/p_k]
  #        + (psi_k - psi) / c              [between; only for estimated weights]
  # with c = n_S / n. Weighting phi_i by w_k alone would understate the variance
  # by roughly w_k^2 per pattern.
  in_set <- !is.na(obs_pattern_pos)
  pos <- obs_pattern_pos[in_set]
  n_k <- tabulate(pos, nbins = length(estimates))
  if (any(n_k == 0)) {
    stop("Cannot combine patterns: an included pattern has no observations.")
  }
  # Re-center within pattern so uncentered inputs (threshold paths) are handled.
  mu_k <- as.numeric(rowsum(obs_eif[in_set], pos)) / n_k
  eif_agg <- rep(0, n_full)
  eif_agg[in_set] <- (weights * n_full / n_k)[pos] * (obs_eif[in_set] - mu_k[pos])
  if (prevalence_weights) {
    c_hat <- sum(n_k) / n_full
    eif_agg[in_set] <- eif_agg[in_set] + (estimates[pos] - psi_agg) / c_hat
  }

  eif_agg_centered <- eif_agg - mean(eif_agg)
  eif_var <- var(eif_agg_centered)
  se_agg <- sqrt(eif_var / n_full)

  alpha <- 1 - confidence_level
  z_alpha <- qnorm(1 - alpha / 2)
  ci_lower <- psi_agg - z_alpha * se_agg
  ci_upper <- psi_agg + z_alpha * se_agg
  z_stat <- psi_agg / se_agg
  p_value <- 2 * (1 - pnorm(abs(z_stat)))

  list(
    estimate = psi_agg,
    se = se_agg,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    eif_values = eif_agg_centered,
    weights = weights,
    coverage = coverage
  )
}

#' Aggregate stratified estimation results
#'
#' @description
#' Aggregates pattern-specific estimates via weighted averaging with delta
#' method variance estimation. Uses aggregated EIF for variance.
#'
#' @param pattern_results List, results from estimate_pattern_effect() for each pattern
#' @param patterns Data frame, pattern information from identify_eligibility_patterns()
#' @param confidence_level Numeric, confidence level for intervals
#'
#' @return List containing:
#' \describe{
#'   \item{estimate}{Numeric, aggregated estimate ψ = Σ p_k ψ_k}
#'   \item{se}{Numeric, standard error from aggregated EIF}
#'   \item{ci_lower}{Numeric, lower confidence bound}
#'   \item{ci_upper}{Numeric, upper confidence bound}
#'   \item{p_value}{Numeric, two-sided p-value for H0: ψ = 0}
#'   \item{eif_values}{Numeric vector (length n), aggregated EIF}
#'   \item{pattern_summaries}{Data frame with per-pattern estimates and weights}
#'   \item{n_patterns_used}{Integer, number of patterns with successful estimation}
#'   \item{coverage}{Numeric, proportion of sample in successfully estimated patterns}
#' }
#'
#' @keywords internal
aggregate_stratified_results <- function(pattern_results, patterns, confidence_level = 0.95) {

  # Extract successful patterns
  success_flags <- sapply(pattern_results, function(x) x$success)

  if (sum(success_flags) == 0) {
    # Collect error messages from failed patterns for diagnostics
    error_msgs <- vapply(pattern_results, function(x) {
      if (!is.null(x$error_message)) x$error_message else "Unknown error"
    }, FUN.VALUE = character(1))
    unique_errors <- unique(error_msgs)
    stop("All patterns failed estimation. Cannot aggregate results.\n",
         "  Unique errors (", length(unique_errors), "):\n",
         paste0("  - ", unique_errors, collapse = "\n"))
  }

  if (sum(success_flags) < nrow(patterns)) {
    # Loud, specific failure report. A silent "N of M failed" hid a one-line bug
    # (single-eligibility patterns failing under sequential PS) for a full
    # estimation cycle, because the scenario still "succeeded" with a quietly
    # smaller denominator. Report WHICH patterns, their size, why, and the
    # realized coverage lost — so this class of failure cannot hide again.
    failed_idx <- which(!success_flags)
    prevalence <- if (!is.null(patterns$prevalence)) patterns$prevalence else rep(NA_real_, nrow(patterns))
    lost_coverage <- sum(prevalence[failed_idx], na.rm = TRUE)
    detail <- vapply(failed_idx, function(i) {
      pr <- pattern_results[[i]]
      sprintf("    - pattern '%s' (n=%s, %s eligible): %s",
              if (!is.null(patterns$pattern_string)) patterns$pattern_string[i] else "?",
              if (!is.null(patterns$n_obs)) patterns$n_obs[i] else "?",
              if (!is.null(patterns$n_eligible)) patterns$n_eligible[i] else "?",
              if (!is.null(pr$error_message) && !is.na(pr$error_message)) pr$error_message else "Unknown error")
    }, FUN.VALUE = character(1))
    warning(sprintf(
      paste0("%d of %d eligibility patterns FAILED estimation and were EXCLUDED ",
             "from the aggregate.\n  Realized coverage dropped to %.1f%% of the ",
             "stratified sample (%.1f%% lost).\n  The reported estimate is over the ",
             "SURVIVING patterns only — treat it as biased if the lost patterns ",
             "differ systematically.\n  Failed patterns:\n%s"),
      length(failed_idx), nrow(patterns),
      100 * (1 - lost_coverage), 100 * lost_coverage,
      paste(detail, collapse = "\n")))
  }

  # Filter to successful patterns
  successful_indices <- which(success_flags)
  pattern_results_success <- pattern_results[successful_indices]
  patterns_success <- patterns[successful_indices, , drop = FALSE]

  estimates <- vapply(pattern_results_success, function(x) x$estimate, FUN.VALUE = numeric(1))
  prevalence_vec <- as.numeric(patterns_success$prevalence)
  eif_matrices <- lapply(pattern_results_success, function(x) x$eif_values_extended)
  n <- length(eif_matrices[[1]])

  # Build the partition representation the combine core expects: for each
  # observation, which included-pattern position it belongs to (NA if none), and
  # its own pattern's UNWEIGHTED influence value. Patterns partition the sample;
  # membership comes from obs_indices (a nonzero test would drop an observation
  # whose centered value is exactly 0).
  obs_pattern_pos <- rep(NA_integer_, n)
  obs_eif <- rep(0, n)
  for (i in seq_along(eif_matrices)) {
    members <- pattern_results_success[[i]]$obs_indices
    if (is.null(members)) {
      members <- which(!is.na(eif_matrices[[i]]) & eif_matrices[[i]] != 0)
    }
    obs_pattern_pos[members] <- i
    obs_eif[members] <- eif_matrices[[i]][members]
  }

  combined <- .combine_patterns(
    estimates = estimates,
    prevalences = prevalence_vec,
    obs_eif = obs_eif,
    obs_pattern_pos = obs_pattern_pos,
    n_full = n,
    confidence_level = confidence_level
  )

  # Per-pattern summary (as before)
  ses <- vapply(pattern_results_success, function(x) x$se, FUN.VALUE = numeric(1))
  pattern_summaries <- data.frame(
    pattern_id = as.integer(patterns_success$pattern_id),
    pattern_string = as.character(patterns_success$pattern_string),
    n_obs = as.integer(patterns_success$n_obs),
    prevalence_original = prevalence_vec,
    weight = combined$weights,
    estimate = estimates,
    se = ses,
    n_eligible = as.integer(patterns_success$n_eligible),
    stringsAsFactors = FALSE
  )

  # Compact, self-contained contributions object for POST-HOC recombination
  # (recombine_stratified_patterns): everything needed to re-aggregate under a
  # different inclusion threshold or reweighting, WITHOUT re-estimating. Stores
  # the per-observation pattern id and its own-pattern unweighted influence.
  pattern_contributions <- list(
    patterns = data.frame(
      pattern_id = as.integer(patterns_success$pattern_id),
      pattern_string = as.character(patterns_success$pattern_string),
      n_obs = as.integer(patterns_success$n_obs),
      prevalence = prevalence_vec,
      n_eligible = as.integer(patterns_success$n_eligible),
      estimate = estimates,
      se = ses,
      stringsAsFactors = FALSE
    ),
    obs_pattern_id = patterns_success$pattern_id[obs_pattern_pos],  # length n; NA if unassigned
    obs_eif = obs_eif,                                              # length n; unweighted
    n_full = n,
    confidence_level = confidence_level
  )

  return(list(
    estimate = combined$estimate,
    se = combined$se,
    ci_lower = combined$ci_lower,
    ci_upper = combined$ci_upper,
    p_value = combined$p_value,
    eif_values = combined$eif_values,
    pattern_summaries = pattern_summaries,
    pattern_contributions = pattern_contributions,
    n_patterns_used = sum(success_flags),
    coverage = combined$coverage
  ))
}

#' Recombine stratified pattern estimates post-hoc (no re-estimation)
#'
#' Re-aggregates a stratified estimate from the per-pattern contributions saved
#' during estimation, under a different inclusion threshold or a custom set of
#' patterns/weights. Because per-pattern estimation is the expensive step
#' (50--100 hrs on the full cohort) and combining is a weighted average plus an
#' influence-function variance, this lets you explore the sample-size threshold
#' (or any reweighting) in seconds without re-running the estimator.
#'
#' @details
#' Each observation belongs to exactly one eligibility pattern, so the aggregate
#' point estimate is the prevalence-weighted mean of per-pattern estimates and the
#' aggregate influence function places each observation's own (weighted) pattern
#' influence on it. Recombination re-selects which patterns are included, ren
#' ormalizes the prevalence weights over that set, and recomputes the estimate,
#' influence-function-based SE, and CI via the same core used during estimation
#' (\code{.combine_patterns()}), so results are identical to a fresh run at the
#' same threshold.
#'
#' \strong{Estimand caveat:} changing which patterns are included changes the
#' covered subpopulation, hence the target parameter --- this is a change of
#' estimand (who is represented), not merely of precision. Lowering the threshold
#' adds smaller/rarer patterns (broader coverage, noisier per-pattern fits);
#' raising it narrows coverage. The realized-coverage of the chosen set is
#' returned so the trade-off is explicit.
#'
#' @param pattern_contributions The \code{pattern_contributions} element returned
#'   by \code{estimate_qdelta_mean_stratified()} (via
#'   \code{aggregate_stratified_results()}). Contains per-pattern estimates and
#'   prevalences, the per-observation pattern id, and per-observation influence
#'   values.
#' @param min_stratum_size Integer. Include only patterns with at least this many
#'   observations. Default: 0 (include all estimated patterns).
#' @param include_patterns Optional integer vector of \code{pattern_id}s to
#'   include, overriding \code{min_stratum_size}. Use for arbitrary pattern sets.
#' @param weights Optional named numeric vector (names = \code{pattern_id}) giving
#'   custom, unnormalized pattern weights, overriding prevalence weighting. Only
#'   patterns in the included set are used; renormalized to sum to 1.
#' @param confidence_level Numeric. Default: the value stored at estimation time.
#' @return A list with \code{estimate}, \code{se}, \code{ci_lower},
#'   \code{ci_upper}, \code{p_value}, \code{n_patterns_used}, \code{coverage}
#'   (fraction of the full stratified sample represented), \code{eif_values}, and
#'   \code{patterns} (the included per-pattern summary with applied weights).
#' @examples
#' \dontrun{
#' fit <- estimate_qdelta_mean_stratified(data, delta, covariate_cols,
#'                                        eligibility_cols, min_stratum_size = 1000)
#' # Explore a lower threshold without re-estimating:
#' recombine_stratified_patterns(fit$pattern_contributions, min_stratum_size = 200)
#' }
#' @export
recombine_stratified_patterns <- function(pattern_contributions,
                                          min_stratum_size = 0,
                                          include_patterns = NULL,
                                          weights = NULL,
                                          confidence_level = NULL) {
  pc <- pattern_contributions
  if (is.null(pc) || is.null(pc$patterns) || is.null(pc$obs_eif)) {
    stop("`pattern_contributions` is missing or malformed. Pass the ",
         "$pattern_contributions element from a stratified fit.")
  }
  if (is.null(confidence_level)) {
    confidence_level <- if (!is.null(pc$confidence_level)) pc$confidence_level else 0.95
  }

  pats <- pc$patterns

  # Select included patterns.
  if (!is.null(include_patterns)) {
    keep <- pats$pattern_id %in% include_patterns
  } else {
    keep <- pats$n_obs >= min_stratum_size
  }
  if (sum(keep) < 1) {
    stop("No patterns satisfy the inclusion criterion ",
         "(min_stratum_size = ", min_stratum_size, ").")
  }
  kept <- pats[keep, , drop = FALSE]

  # Prevalence weights by default; custom weights if supplied.
  if (!is.null(weights)) {
    w <- weights[as.character(kept$pattern_id)]
    if (anyNA(w)) {
      stop("`weights` must supply a value (named by pattern_id) for every ",
           "included pattern; missing: ",
           paste(kept$pattern_id[is.na(w)], collapse = ", "))
    }
    prevalences <- as.numeric(w)
  } else {
    prevalences <- kept$prevalence
  }

  # Map each observation to a position within `kept` (NA if its pattern excluded).
  pos_by_id <- setNames(seq_len(nrow(kept)), kept$pattern_id)
  obs_pattern_pos <- unname(pos_by_id[as.character(pc$obs_pattern_id)])

  combined <- .combine_patterns(
    estimates = kept$estimate,
    prevalences = prevalences,
    obs_eif = pc$obs_eif,
    obs_pattern_pos = obs_pattern_pos,
    n_full = pc$n_full,
    confidence_level = confidence_level,
    prevalence_weights = is.null(weights),
    coverage = sum(kept$prevalence)
  )

  kept$weight <- combined$weights
  list(
    estimate = combined$estimate,
    se = combined$se,
    ci_lower = combined$ci_lower,
    ci_upper = combined$ci_upper,
    p_value = combined$p_value,
    eif_values = combined$eif_values,
    n_patterns_used = nrow(kept),
    coverage = combined$coverage,
    patterns = kept
  )
}

#' Contrast between two stratified fits (difference in counterfactual means)
#'
#' @description
#' Difference \eqn{\hat\psi_b - \hat\psi_a} between two stratified fits that
#' share the same sample, with an SE from the influence function of the
#' difference, \eqn{\mathrm{sd}(\mathrm{IF}_b - \mathrm{IF}_a)/\sqrt{n}}. This
#' accounts for the (typically large) positive covariance between the two
#' estimates, so a contrast against \eqn{\delta = 1} can be far more precise than
#' either level.
#'
#' @details
#' The covariance is only valid when the two fits used the SAME cross-fitting
#' splits (shared nuisance fits across delta scenarios). With independent splits
#' the influence values carry independent split noise and this SE is not valid.
#' The function checks that both fits cover the same observations and patterns
#' but cannot verify that the splits were shared; confirm it from the run
#' provenance (e.g. \code{nuisance_sharing} in the checkpoint signature).
#'
#' @param contributions_a,contributions_b \code{pattern_contributions} lists from
#'   two stratified fits on the same sample (\eqn{\psi_b - \psi_a} is returned).
#' @param min_stratum_size,include_patterns Passed to
#'   \code{recombine_stratified_patterns()}; applied identically to both fits.
#' @param confidence_level Numeric, default 0.95.
#' @return List: \code{estimate} (b minus a), \code{se}, \code{ci_lower},
#'   \code{ci_upper}, \code{z}, \code{p_value}, and the two level estimates
#'   \code{estimate_a}, \code{estimate_b}.
#' @seealso \code{\link{recombine_stratified_patterns}}
#' @export
stratified_contrast <- function(contributions_a, contributions_b,
                                min_stratum_size = 0, include_patterns = NULL,
                                confidence_level = 0.95) {
  if (!identical(contributions_a$n_full, contributions_b$n_full) ||
      !identical(as.integer(contributions_a$obs_pattern_id),
                 as.integer(contributions_b$obs_pattern_id))) {
    stop("The two fits must cover the same observations with the same ",
         "pattern assignment (same n_full and obs_pattern_id).")
  }
  ra <- recombine_stratified_patterns(contributions_a, min_stratum_size,
                                      include_patterns, confidence_level = confidence_level)
  rb <- recombine_stratified_patterns(contributions_b, min_stratum_size,
                                      include_patterns, confidence_level = confidence_level)
  d_if <- rb$eif_values - ra$eif_values
  est <- rb$estimate - ra$estimate
  se <- sqrt(var(d_if) / length(d_if))
  z_alpha <- qnorm(1 - (1 - confidence_level) / 2)
  z <- est / se
  list(
    estimate = est,
    se = se,
    ci_lower = est - z_alpha * se,
    ci_upper = est + z_alpha * se,
    z = z,
    p_value = 2 * (1 - pnorm(abs(z))),
    estimate_a = ra$estimate,
    estimate_b = rb$estimate
  )
}

#' Estimate counterfactual mean with eligibility stratification
#'
#' @description
#' Estimates `E[Y^{Q_delta}]` when treatments are only defined for subsets of
#' observations (e.g., medication adherence only measurable if prescribed).
#' Stratifies by eligibility pattern, estimates within each pattern, and
#' aggregates via weighted averaging with delta method inference.
#'
#' **When to use this function:**
#' - Treatments are only defined/eligible for some observations
#' - You have eligibility indicators (elig_A1, elig_A2, etc.)
#' - Want to estimate: "What if we improved treatment quality/adherence
#'   (by factor δ) conditional on current eligibility patterns?"
#'
#' **When NOT to use this function:**
#' - All treatments eligible for all observations → use estimate_qdelta_mean()
#' - Eligibility is rare and you want to include as covariate → use estimate_qdelta_mean()
#' - Interested in changing who gets eligible → different estimand (out of scope)
#'
#' @param data Data frame or tibble with columns:
#'   - Treatment columns: A1, A2, ..., Am (binary 0/1)
#'   - Eligibility columns: elig_A1, elig_A2, ..., elig_Am (binary 0/1)
#'   - Covariate columns: X1, X2, ..., Xp (numeric)
#'   - Outcome column: Y (numeric)
#' @param delta Numeric vector of length m, target marginal changes (odds ratios).
#'   Applied only to eligible treatments within each pattern.
#' @param covariate_cols Character vector of covariate column names (required).
#' @param eligibility_cols Character vector of eligibility column names (required).
#'   Length defines m; treatment columns are derived per pattern.
#' @param outcome_col Character scalar, outcome column name. Default "Y".
#' @param divergence Character, divergence measure. Options: "kl" (recommended),
#'   "reverse_kl", "chi2", "hellinger". Default: "kl".
#' @param min_stratum_size Integer, minimum observations per pattern to include
#'   in estimation. Patterns with fewer observations are excluded. Default: 100.
#' @param min_realized_coverage Numeric in `[0, 1]`, minimum fraction of the
#'   stratified sample that must be covered by SUCCESSFULLY-estimated patterns.
#'   Unlike the pre-estimation feasibility check (which only verifies identified
#'   patterns meet `min_stratum_size`), this reflects estimation failures: if
#'   patterns fail and realized coverage drops below this threshold, the function
#'   aborts rather than silently returning an estimate over a non-representative
#'   subset. Set to 0 to disable the gate. Default: 0.8.
#' @param max_patterns Integer, maximum number of patterns to use. If more patterns
#'   exist, uses only the top max_patterns by prevalence. NULL = use all feasible.
#'   Default: NULL.
#' @param aggregate_method Character, aggregation method. Currently only "weighted"
#'   is supported (weighted by pattern prevalence). Default: "weighted".
#' @param ps_method Character, propensity score estimation method.
#'   Options: "glm", "randomforest", "xgboost", "binary_per_combination".
#'   Default: "glm".
#' @param binary_method Character, base learner for binary_per_combination and
#'   sequential propensity methods. Options: "glm", "randomforest", "xgboost".
#'   Default: "glm".
#' @param outcome_method Character, outcome model estimation method.
#'   Options: "glm", "gam", "randomforest", "superlearner". Default: "glm".
#' @param estimator Character, final estimator type.
#'   Options: "one_step" (recommended, the estimator analyzed in the paper),
#'   "tmle" (EXPERIMENTAL, no supporting theory). Default: "one_step".
#' @param cross_fit Logical, use cross-fitting for nuisance estimation.
#'   Default: FALSE.
#' @param n_folds Integer, number of folds for cross-fitting. Only used if
#'   cross_fit = TRUE. Automatically disabled for patterns with insufficient
#'   sample size. Default: 5.
#' @param confidence_level Numeric in (0, 1), confidence level for intervals.
#'   Default: 0.95.
#' @param threshold_type Character, active-set selection for handling near-boundary
#'   marginal propensities within each eligibility pattern. One of "none" (default;
#'   standard estimator), "vanishing" (vanishing threshold ε_n → 0, targets ψ(δ)
#'   asymptotically), or "adaptive" (fixed threshold ε with adaptive delta, targets
#'   the localized-margin estimand). Applied per-coordinate on m_k(X) WITHIN each
#'   pattern (eligibility subsetting is separate from and prior to thresholding).
#'   Forwarded to \code{\link{estimate_pattern_effect}}. Default "none" preserves
#'   prior behavior.
#' @param gamma Numeric, threshold rate for the vanishing estimator (`ε_n = C n^{-γ}`).
#'   Only used when threshold_type = "vanishing". Default NULL (auto-selected).
#' @param epsilon Numeric, fixed threshold for the adaptive estimator. Only used
#'   when threshold_type = "adaptive". Default 0.1.
#' @param bootstrap Logical, use bootstrap for inference. Not currently supported
#'   for stratified estimation. Default: FALSE.
#' @param n_bootstrap Integer, number of bootstrap resamples. Not used.
#' @param verbose Logical, print progress messages. Default: TRUE.
#' @param nuisance_cache Optional \code{nuisance_cache} from
#'   \code{\link{new_nuisance_cache}}, shared across calls to avoid refitting
#'   nuisances that do not depend on \code{delta}. Within one cache, a pattern's
#'   propensity fit is computed once and reused by every subsequent call (any
#'   delta, any cost type); its outcome fit is reused by every call with the same
#'   \code{nuisance_tag} until \code{\link{.clear_outcomes}} is called. Fits use
#'   namespaced per-fit seeds and are RNG-neutral to the caller's stream, so
#'   results do not depend on which cells a resumed sweep recomputes. Only
#'   supported with \code{threshold_type = "none"} (error otherwise, rather than
#'   silently ignoring the cache). Default NULL (fit everything per call, as
#'   before).
#' @param nuisance_tag Character scalar identifying the outcome/cost type for
#'   cache keying (e.g. \code{"total"}). Ignored when \code{nuisance_cache} is
#'   NULL.
#'
#' @param max_qdelta_failure_rate Numeric in `[0, 1]`: the largest share of observations whose
#'   Q_delta solve may fail and fall back to Q_delta = P0 (no shift for that observation, a bias
#'   toward the null). Default `0`: the first failure stops the estimate with the observation index
#'   and the solver's message. A positive value tolerates failures up to `rate * n` observations
#'   (each warned about and counted), then stops. `1` always falls back. The vanishing-threshold
#'   estimator never falls back and always stops.
#' @param qdelta_tol Numeric or `NULL`. Solver tolerance for the per-observation Q_delta problem;
#'   `NULL` (default) uses `QDELTA_TOLERANCE` for every sample size (it no longer relaxes above
#'   100,000 observations, which made the same stratum's estimate depend on whether it crossed that
#'   size).
#' @return List with class "qdelta_stratified" containing:
#' \describe{
#'   \item{estimate}{Numeric, aggregated estimate ψ = Σ p_k ψ_k}
#'   \item{se}{Numeric, standard error from delta method}
#'   \item{ci_lower}{Numeric, lower confidence bound}
#'   \item{ci_upper}{Numeric, upper confidence bound}
#'   \item{p_value}{Numeric, two-sided p-value for H0: ψ = 0}
#'   \item{eif_values}{Numeric vector (length n), centered aggregated EIF}
#'   \item{pattern_summaries}{Data frame with per-pattern estimates and weights}
#'   \item{eligibility_patterns}{List from identify_eligibility_patterns()}
#'   \item{feasibility_check}{List from validate_stratification_feasibility()}
#'   \item{n_patterns_used}{Integer, number of patterns with successful estimation}
#'   \item{coverage}{Numeric, proportion of sample in successfully estimated patterns}
#'   \item{delta}{Numeric vector, input delta values}
#'   \item{divergence}{Character, divergence measure used}
#'   \item{estimator}{Character, estimator used}
#'   \item{n}{Integer, total sample size}
#'   \item{m}{Integer, number of treatments}
#' }
#'
#' @details
#' **Workflow:**
#' 1. Identifies unique eligibility patterns in the data
#' 2. Validates that stratification is feasible (enough patterns, enough observations)
#' 3. For each pattern k:
#'    - Subsets to observations with pattern k
#'    - Filters to eligible treatments only
#'    - Renames treatments to A1, ..., A_m_k (consecutive)
#'    - Calls estimate_qdelta_mean() with reduced dimension
#'    - Extends EIF to full dataset (zeros for other observations)
#' 4. Aggregates: ψ = Σ_k p_k ψ_k with delta method variance
#'
#' **Estimand interpretation:**
#' The estimated quantity is `E[Y^{Q_delta | Z}]` where Z = (Z_1, ..., Z_m)
#' are eligibility indicators. This answers: "What if we improved treatment
#' quality/adherence by factor δ_k for treatments k that are currently eligible,
#' without changing eligibility itself?"
#'
#' **Treatment renaming:**
#' Within each pattern, only eligible treatments are used. These are renamed
#' to consecutive indices A1, A2, ... for estimation, then results are mapped
#' back to original indices.
#'
#' **Variance estimation:**
#' The pooled estimator is the mean of the observation-level influence values,
#' so its influence function is phi_i - psi (see \code{.combine_patterns}). Each
#' pattern's influence function is rescaled by 1/p_k and a between-pattern term
#' is added, so the SE reflects both within- and between-pattern variability.
#'
#' @seealso
#' \code{\link{estimate_qdelta_mean}} for standard (non-stratified) estimation,
#' \code{\link{profile_eligibility_patterns}} for assessing stratification feasibility
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate example data with eligibility patterns
#' set.seed(123)
#' n <- 500; m <- 2; p <- 3
#'
#' # Simulate data
#' X <- matrix(rnorm(n * p), nrow = n)
#' colnames(X) <- paste0("X", 1:p)
#'
#' # Eligibility depends on covariates
#' elig_A1 <- rbinom(n, 1, plogis(X[,1]))
#' elig_A2 <- rbinom(n, 1, plogis(X[,2]))
#'
#' # Treatment only observed if eligible
#' A1 <- ifelse(elig_A1 == 1, rbinom(n, 1, 0.3), 0)
#' A2 <- ifelse(elig_A2 == 1, rbinom(n, 1, 0.4), 0)
#'
#' # Outcome depends on treatments and covariates
#' Y <- rnorm(n, mean = 2 + A1 + 1.5*A2 + 0.5*X[,1])
#'
#' # Combine into data frame
#' data <- data.frame(
#'   A1 = A1, A2 = A2,
#'   elig_A1 = elig_A1, elig_A2 = elig_A2,
#'   X, Y
#' )
#'
#' # Estimate with eligibility stratification
#' result <- estimate_qdelta_mean_stratified(
#'   data = data,
#'   delta = c(1.5, 2.0),
#'   covariate_cols = paste0("X", 1:p),
#'   eligibility_cols = c("elig_A1", "elig_A2"),
#'   divergence = "kl",
#'   min_stratum_size = 50
#' )
#'
#' # Print results
#' print(result$estimate)
#' print(result$se)
#' print(c(result$ci_lower, result$ci_upper))
#'
#' # Check pattern summaries
#' print(result$pattern_summaries)
#'
#' # Compare to non-stratified (eligibility as covariates)
#' result_std <- estimate_qdelta_mean(data, delta = c(1.5, 2.0),
#'   treatment_cols = c("A1", "A2"), covariate_cols = paste0("X", 1:p))
#' print(c(stratified = result$estimate, standard = result_std$estimate))
#' }
estimate_qdelta_mean_stratified <- function(data, delta,
                                           covariate_cols,
                                           eligibility_cols,
                                           outcome_col = "Y",
                                           divergence = "kl",
                                           min_stratum_size = 100,
                                           min_realized_coverage = 0.8,
                                           max_patterns = NULL,
                                           aggregate_method = "weighted",
                                           ps_method = "glm",
                                           binary_method = "glm",
                                           outcome_method = "glm",
                                           estimator = "one_step",
                                           cross_fit = TRUE,
                                           n_folds = 5,
                                           confidence_level = 0.95,
                                           threshold_type = c("none", "vanishing", "adaptive"),
                                           gamma = NULL,
                                           epsilon = 0.1,
                                           bootstrap = FALSE,
                                           n_bootstrap = 1000,
                                           verbose = TRUE,
                                           xgb_params = NULL,
                                           nuisance_cache = NULL,
                                           nuisance_tag = NULL,
    max_qdelta_failure_rate = 0, qdelta_tol = NULL) {

  threshold_type <- match.arg(threshold_type)

  # Nuisance-fit cache: validate here too so the failure is immediate and names
  # this function, instead of surfacing from the first pattern.
  if (!is.null(nuisance_cache)) {
    .check_nuisance_cache(nuisance_cache)
    if (threshold_type != "none") {
      stop("`nuisance_cache` is only supported with threshold_type = 'none'. ",
           "threshold_type = '", threshold_type, "' routes each pattern to a ",
           "thresholded estimator that fits its own nuisances and ignores ",
           "pre-supplied ones, so the cache would be silently ignored.",
           call. = FALSE)
    }
  }

  # Input validation
  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a positive numeric vector")
  }

  if (bootstrap) {
    stop("Bootstrap inference not currently supported for stratified estimation. ",
         "Use bootstrap = FALSE (delta method inference via aggregated EIF).")
  }

  if (aggregate_method != "weighted") {
    stop("Currently only aggregate_method = 'weighted' is supported.")
  }

  # Explicit column specifications (no auto-detection). Treatment columns are
  # derived per eligibility pattern (a pattern-specific active subset), so only
  # covariate/eligibility/outcome columns are validated up front here. The
  # per-pattern treatment columns are validated inside estimate_pattern_effect.
  if (is.null(eligibility_cols) || length(eligibility_cols) == 0) {
    stop("`eligibility_cols` must be provided (e.g. c('elig_A1', 'elig_A2', ...)).",
         call. = FALSE)
  }
  if (is.null(covariate_cols) || length(covariate_cols) == 0) {
    stop("`covariate_cols` must be provided (e.g. c('X1', 'X2', ...)).", call. = FALSE)
  }
  missing_cols <- setdiff(c(covariate_cols, outcome_col, eligibility_cols), names(data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in `data`: ", paste(missing_cols, collapse = ", "), ".",
         call. = FALSE)
  }

  # Step 1: Identify eligibility patterns
  if (verbose) {
    message("Step 1: Identifying eligibility patterns...")
  }

  pattern_info <- identify_eligibility_patterns(data, eligibility_cols)
  m <- pattern_info$m

  if (length(delta) != m) {
    stop("delta length (", length(delta), ") must equal number of treatments (", m, ")")
  }

  if (verbose) {
    message(sprintf("  Found %d unique eligibility patterns", nrow(pattern_info$patterns)))
  }

  # Step 2: Validate feasibility
  if (verbose) {
    message("Step 2: Validating stratification feasibility...")
  }

  feasibility_check <- validate_stratification_feasibility(
    pattern_info$patterns,
    min_stratum_size = min_stratum_size,
    min_coverage = 0.5
  )

  if (!feasibility_check$feasible) {
    stop("Stratification not feasible: ", feasibility_check$reason, "\n",
         "Consider: (1) reducing min_stratum_size, (2) using non-stratified estimation ",
         "with eligibility as covariates (estimate_qdelta_mean), or (3) collapsing patterns.")
  }

  if (verbose) {
    message(sprintf("  Feasible: %d patterns with n >= %d (%.1f%% coverage)",
                   feasibility_check$n_feasible_patterns, min_stratum_size,
                   100 * feasibility_check$coverage))
  }

  # Step 3: Filter to feasible patterns
  feasible_pattern_ids <- feasibility_check$feasible_pattern_ids
  patterns_to_use <- pattern_info$patterns[pattern_info$patterns$pattern_id %in% feasible_pattern_ids, ]

  # Optional: Limit to top max_patterns by prevalence
  if (!is.null(max_patterns) && nrow(patterns_to_use) > max_patterns) {
    if (verbose) {
      message(sprintf("  Limiting to top %d patterns by prevalence", max_patterns))
    }
    patterns_to_use <- patterns_to_use[1:max_patterns, ]
  }

  n_patterns_to_estimate <- nrow(patterns_to_use)

  if (verbose) {
    message(sprintf("Step 3: Estimating within %d patterns...", n_patterns_to_estimate))
  }

  # Step 4: Estimate within each pattern
  pattern_results <- vector("list", n_patterns_to_estimate)

  # Per-phase timing. Nuisance-fit seconds are read off the cache's cumulative
  # counters (fits happen inside .pattern_nuisances); everything else in the
  # pattern loop is Q_delta + EIF work. Without a cache, nuisance fitting happens
  # inside estimate_qdelta_mean and is not separable here, so it is reported NA.
  pattern_loop_start <- Sys.time()
  cache_secs_before <- if (!is.null(nuisance_cache)) {
    nuisance_cache$ps_fit_seconds + nuisance_cache$outcome_fit_seconds
  } else NA_real_

  for (i in seq_len(n_patterns_to_estimate)) {
    pattern_i_df <- patterns_to_use[i, ]

    # Convert data frame row to proper list (handle list columns correctly)
    pattern_i <- list(
      pattern_id = pattern_i_df$pattern_id,
      pattern_string = pattern_i_df$pattern_string,
      eligibility = pattern_i_df$eligibility[[1]],  # Extract vector from list column
      n_obs = pattern_i_df$n_obs,
      prevalence = pattern_i_df$prevalence,
      n_eligible = pattern_i_df$n_eligible,
      obs_indices = pattern_i_df$obs_indices[[1]]  # Extract vector from list column
    )

    if (verbose) {
      message(sprintf("  Pattern %d/%d: '%s' (n=%d, %.1f%%, %d eligible treatments)",
                     i, n_patterns_to_estimate, pattern_i$pattern_string,
                     pattern_i$n_obs, 100 * pattern_i$prevalence, pattern_i$n_eligible))
    }

    pattern_results[[i]] <- estimate_pattern_effect(
      data = data,
      pattern_info = pattern_i,
      delta = delta,
      covariate_cols = covariate_cols,
      outcome_col = outcome_col,
      ps_method = ps_method,
      binary_method = binary_method,
      outcome_method = outcome_method,
      estimator = estimator,
      cross_fit = cross_fit,
      n_folds = n_folds,
      divergence = divergence,
      threshold_type = threshold_type,
      gamma = gamma,
      epsilon = epsilon,
      verbose = verbose,
      xgb_params = xgb_params,
      nuisance_cache = nuisance_cache,
      nuisance_tag = nuisance_tag,
      pattern_pos = i,
      max_qdelta_failure_rate = max_qdelta_failure_rate,
      qdelta_tol = qdelta_tol
    )

    if (!pattern_results[[i]]$success) {
      warning(sprintf("Pattern %d ('%s') estimation failed: %s",
                     i, pattern_i$pattern_string, pattern_results[[i]]$error_message))
    } else {
      if (verbose) {
        message(sprintf("    Estimate: %.4f (SE: %.4f)",
                       pattern_results[[i]]$estimate, pattern_results[[i]]$se))
      }
    }
  }

  # Close out per-phase timing for this call (see pattern_loop_start above).
  pattern_loop_seconds <- as.numeric(difftime(Sys.time(), pattern_loop_start,
                                             units = "secs"))
  nuisance_fit_seconds <- if (!is.null(nuisance_cache)) {
    (nuisance_cache$ps_fit_seconds + nuisance_cache$outcome_fit_seconds) -
      cache_secs_before
  } else NA_real_
  nuisance_timing <- list(
    pattern_loop_seconds = pattern_loop_seconds,
    nuisance_fit_seconds = nuisance_fit_seconds,
    # Everything in the loop that is not nuisance fitting: Q_delta solving, EIF,
    # (optional) TMLE targeting, aggregation bookkeeping.
    qdelta_eif_seconds = if (is.na(nuisance_fit_seconds)) NA_real_ else
      pattern_loop_seconds - nuisance_fit_seconds,
    cache_used = !is.null(nuisance_cache),
    cache_stats = if (!is.null(nuisance_cache)) nuisance_cache_stats(nuisance_cache) else NULL
  )
  if (verbose && !is.null(nuisance_cache)) {
    message(sprintf("  Pattern loop: %.1fs total = %.1fs nuisance fitting + %.1fs Q_delta/EIF",
                    pattern_loop_seconds, nuisance_fit_seconds,
                    nuisance_timing$qdelta_eif_seconds))
  }

  # Step 5: Aggregate results
  if (verbose) {
    message("Step 4: Aggregating results...")
  }

  aggregated_results <- aggregate_stratified_results(
    pattern_results = pattern_results,
    patterns = patterns_to_use,
    confidence_level = confidence_level
  )

  # Realized coverage: fraction of the FULL stratified sample actually represented
  # by successfully-estimated patterns. Distinct from the pre-estimation feasibility
  # coverage (which only checks that identified patterns meet min_stratum_size) --
  # this reflects estimation FAILURES too. Guards against silently reporting an
  # estimate computed over a small, non-representative slice of the cohort.
  succeeded <- vapply(pattern_results, function(x) isTRUE(x$success), logical(1))
  n_estimated <- sum(vapply(pattern_results[succeeded], function(x) x$n_obs, numeric(1)))
  n_stratified <- sum(patterns_to_use$n_obs)
  realized_coverage <- if (n_stratified > 0) n_estimated / n_stratified else 0

  if (verbose) {
    message(sprintf("  Aggregated estimate: %.4f (SE: %.4f)",
                   aggregated_results$estimate, aggregated_results$se))
    message(sprintf("  95%% CI: [%.4f, %.4f]",
                   aggregated_results$ci_lower, aggregated_results$ci_upper))
    message(sprintf("  Used %d/%d patterns (%.1f%% of stratified sample)",
                   aggregated_results$n_patterns_used, n_patterns_to_estimate,
                   100 * realized_coverage))
  }

  # Realized-coverage gate: abort rather than silently return an estimate over a
  # non-representative subset. This is the guardrail that would have surfaced the
  # single-eligibility failure immediately (46.6% of the cohort was dropped).
  if (realized_coverage < min_realized_coverage) {
    dropped <- patterns_to_use[!succeeded, , drop = FALSE]
    dropped_desc <- if (nrow(dropped) > 0) {
      paste(sprintf("'%s' (n=%s, %s eligible)",
                    dropped$pattern_string, dropped$n_obs, dropped$n_eligible),
            collapse = ", ")
    } else "(none identified)"
    stop(sprintf(
      paste0("Realized pattern coverage %.1f%% is below min_realized_coverage %.1f%%: ",
             "%d of %d patterns failed estimation, dropping %.1f%% of the stratified ",
             "sample.\n  The aggregate would be computed over a non-representative subset.\n",
             "  Dropped patterns: %s\n  Inspect the per-pattern failure warnings above, ",
             "fix the cause, or lower min_realized_coverage deliberately if the drop is acceptable."),
      100 * realized_coverage, 100 * min_realized_coverage,
      sum(!succeeded), length(succeeded), 100 * (1 - realized_coverage), dropped_desc))
  }

  # Return structured output
  result <- list(
    estimate = aggregated_results$estimate,
    se = aggregated_results$se,
    ci_lower = aggregated_results$ci_lower,
    ci_upper = aggregated_results$ci_upper,
    p_value = aggregated_results$p_value,
    eif_values = aggregated_results$eif_values,
    pattern_summaries = aggregated_results$pattern_summaries,
    pattern_results = pattern_results,
    eligibility_patterns = pattern_info,
    feasibility_check = feasibility_check,
    n_patterns_used = aggregated_results$n_patterns_used,
    coverage = aggregated_results$coverage,
    realized_coverage = realized_coverage,
    n_patterns_failed = sum(!succeeded),
    pattern_contributions = aggregated_results$pattern_contributions,
    nuisance_timing = nuisance_timing,
    delta = delta,
    divergence = divergence,
    estimator = estimator,
    threshold_type = threshold_type,
    n = nrow(data),
    m = m
  )

  class(result) <- c("qdelta_stratified", "list")

  return(result)
}


#' Estimate Vector Incremental Effects with Vanishing Threshold
#'
#' Implements vanishing threshold estimator with threshold-based active set
#' selection. Uses ε_n → 0 with active sets `K_{ε_n}(X) = {k: m_k(X) ∈ [ε_n, 1-ε_n]}`
#' to satisfy margin conditions required for asymptotic normality.
#'
#' @inheritParams estimate_qdelta_mean
#' @param gamma Numeric, threshold rate parameter. `ε_n = C*n^{-γ}`. If NULL
#'   (default), automatically selects γ = 0.3. The asymptotic-normality condition
#'   γ > 1/(2(1+β)) (β = margin exponent of the data law) is not enforced: β is
#'   not estimable in practice and does not affect whether the estimate is
#'   computable. Larger γ gives a faster-vanishing threshold.
#' @param min_active Integer, minimum number of active treatments required per
#'   observation. Observations with fewer active treatments contribute zero to EIF.
#'   Default: 1.
#' @param eps_constant Numeric, constant C in `ε_n = C*n^{-γ}`. Default: 1.
#'
#' @param max_qdelta_failure_rate Numeric in `[0, 1]`: the largest share of observations whose
#'   Q_delta solve may fail and fall back to Q_delta = P0 (no shift for that observation, a bias
#'   toward the null). Default `0`: the first failure stops the estimate with the observation index
#'   and the solver's message. A positive value tolerates failures up to `rate * n` observations
#'   (each warned about and counted), then stops. `1` always falls back. The vanishing-threshold
#'   estimator never falls back and always stops.
#' @param qdelta_tol Numeric or `NULL`. Solver tolerance for the per-observation Q_delta problem;
#'   `NULL` (default) uses `QDELTA_TOLERANCE` for every sample size (it no longer relaxes above
#'   100,000 observations, which made the same stratum's estimate depend on whether it crossed that
#'   size).
#' @return List with same structure as \code{\link{estimate_qdelta_mean}}, plus:
#' \describe{
#'   \item{threshold_info}{List with epsilon_n, gamma, min_active, eps_constant}
#'   \item{active_set_summary}{Summary of active sets (from \code{\link{summarize_active_sets}})}
#'   \item{n_insufficient_active}{Number of observations with insufficient active treatments}
#' }
#'
#' @details
#' **Vanishing threshold approach:**
#'
#' The standard estimator assumes all treatments are active (m_k(X) bounded away from 0 and 1).
#' When treatments have marginal propensities near boundaries, this assumption fails.
#'
#' Solution: Only use "active" treatments where `m_k(X) ∈ [ε_n, 1-ε_n]` with ε_n → 0.
#' As n → ∞:
#' - ε_n → 0 (more treatments become active)
#' - But ε_n decreases slowly enough to maintain margin condition
#'
#' **Rate conditions (from paper theory):**
#'
#' For asymptotic normality, need:
#' 1. γ > 1/(2β) where β is nuisance rate
#' 2. `α ≥ max{1/4, 1/2 - γ(β-1)}` where α is PS rate
#'
#' **Automatic selection (if gamma/alpha not provided):**
#' - γ = 0.3 (safe for β ≥ 2, i.e., parametric nuisances)
#' - α = 0.4 (safe for most settings with γ = 0.3)
#'
#' **Insufficient active treatments:**
#'
#' Observations with < min_active treatments in `[ε_n, 1-ε_n]` contribute zero
#' to the EIF. The user should verify this fraction is acceptable or adjust ε_n.
#'
#' **Comparison to standard estimator:**
#' - Standard: Uses all treatments, assumes positivity everywhere
#' - Vanishing threshold: Adapts to propensity patterns, more robust
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data with some near-boundary propensities
#' sim_data <- simulate_vector_treatment_data(n = 1000, m = 3, p = 2, seed = 123)
#' tcols <- paste0("A", 1:3)
#' xcols <- paste0("X", 1:2)
#'
#' # Standard estimator (may be unstable with boundary propensities)
#' result_std <- estimate_qdelta_mean(
#'   data = sim_data$data,
#'   delta = c(0.1, 0.1, 0.1),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   divergence = "kl"
#' )
#'
#' # Vanishing threshold estimator (automatic rate selection)
#' result_vt <- estimate_qdelta_vanishing_threshold(
#'   data = sim_data$data,
#'   delta = c(0.1, 0.1, 0.1),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   divergence = "kl"
#' )
#'
#' # Check active set diagnostics
#' print(result_vt$active_set_summary)
#' print(result_vt$n_insufficient_active)
#'
#' # Custom threshold rate (slower-vanishing threshold)
#' result_vt_custom <- estimate_qdelta_vanishing_threshold(
#'   data = sim_data$data,
#'   delta = c(0.1, 0.1, 0.1),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   divergence = "kl",
#'   gamma = 0.25   # Slower threshold vanishing
#' )
#' }
estimate_qdelta_vanishing_threshold <- function(data, delta,
                                                 treatment_cols, covariate_cols,
                                                 outcome_col = "Y", divergence = "kl",
                                                 ps_method = "glm", binary_method = "glm",
                                                 outcome_method = "glm",
                                                 estimator = "one_step",
                                                 cross_fit = TRUE, n_folds = 5,
                                                 confidence_level = 0.95,
                                                 bootstrap = FALSE, n_bootstrap = 1000,
                                                 gamma = NULL,
                                                 min_active = 1, eps_constant = 1,
    max_qdelta_failure_rate = 0, qdelta_tol = NULL) {

  # Input validation
  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a positive numeric vector")
  }

  # Explicit column specifications (no auto-detection).
  validate_columns(data, treatment_cols, covariate_cols, outcome_col)
  m <- length(treatment_cols)

  if (length(delta) != m) {
    stop("delta length must equal number of treatments")
  }

  # Threshold rate. Only gamma (and n) determine the threshold epsilon_n; the
  # asymptotic-normality condition gamma > 1/(2(1+beta)) depends on the margin
  # exponent beta of the data-generating law, which is not estimable in practice
  # and does not affect whether the estimate is computable. We therefore do not
  # gate estimation on it (see quality_reports/plans/2026-07-04_remove-rate-guard.md).
  if (is.null(gamma)) {
    gamma <- 0.3  # Reasonable default threshold decay rate
    message("Using automatic gamma = 0.3")
  }

  # Estimate nuisance parameters (same as standard estimator)
  ps_est <- estimate_propensity_scores(data,
                                       treatment_cols = treatment_cols,
                                       covariate_cols = covariate_cols,
                                       method = ps_method,
                                       binary_method = binary_method,
                                       cross_fit = cross_fit, n_folds = n_folds)
  gc(verbose = FALSE)

  outcome_est <- estimate_outcome_models(data,
                                         treatment_cols = treatment_cols,
                                         covariate_cols = covariate_cols,
                                         outcome_col = outcome_col,
                                         method = outcome_method,
                                         cross_fit = cross_fit, n_folds = n_folds)
  gc(verbose = FALSE)

  # Validate propensity scores
  if (any(!is.finite(ps_est$propensity_scores))) {
    stop("Propensity scores contain non-finite values (NA, NaN, or Inf). ",
         "Check your propensity score estimation method and data quality.")
  }

  ps_row_sums <- rowSums(ps_est$propensity_scores)
  if (!all(abs(ps_row_sums - 1) < PROPENSITY_SUM_TOL)) {
    max_deviation <- max(abs(ps_row_sums - 1))
    stop("Propensity scores do not sum to 1 (max deviation: ", sprintf("%.2e", max_deviation), "). ",
         "Check propensity score estimation.")
  }

  # Validate outcome means
  if (any(!is.finite(outcome_est$outcome_means))) {
    stop("Outcome means contain non-finite values. Check outcome model estimation.")
  }

  # Extract marginal propensities (needed for active set computation)
  A_mat <- as.matrix(data[, treatment_cols, drop = FALSE])
  marginal_probs <- extract_marginal_probabilities(ps_est$propensity_scores, m = m)

  # Compute vanishing threshold
  n <- nrow(data)
  epsilon_n <- compute_vanishing_threshold(n, gamma, constant = eps_constant)
  message(sprintf("Computed vanishing threshold: epsilon_n = %.4e (n=%d, gamma=%.2f)",
                 epsilon_n, n, gamma))

  # Compute active sets
  active_sets <- compute_active_sets(marginal_probs, epsilon = epsilon_n, min_active = min_active)

  # Summarize active sets
  active_summary <- summarize_active_sets(active_sets)
  message(sprintf("Active set summary: %.1f%% obs have all %d treatments active, median = %d active",
                 100 * active_summary$pct_all_active, m, active_summary$median_n_active))

  # Warn if many observations have insufficient active treatments
  n_insufficient <- sum(!active_sets$sufficient_active)
  if (n_insufficient > 0) {
    pct_insufficient <- 100 * n_insufficient / n
    if (pct_insufficient > 10) {
      warning(sprintf("%.1f%% of observations have insufficient active treatments (< %d).\n",
                     pct_insufficient, min_active,
                     "Consider increasing epsilon_n (reduce gamma) or decreasing min_active."))
    } else {
      message(sprintf("%.1f%% of observations have insufficient active treatments", pct_insufficient))
    }
  }

  # Compute per-observation Q_delta (same as standard estimator)
  qdelta_results <- vector("list", n)

  # One tolerance for every sample size (see qdelta_tol)
  adaptive_tol <- resolve_qdelta_tol(qdelta_tol)

  # Progress reporting for large datasets
  report_interval <- if (n > QDELTA_LARGE_N_THRESHOLD) 10000 else Inf
  start_time <- Sys.time()

  for (i in 1:n) {
    if (i %% report_interval == 0) {
      pct <- 100 * i / n
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      rate <- i / elapsed
      remaining <- (n - i) / rate
      message(sprintf("Q_delta: %d / %d (%.1f%%) | %.1f obs/sec | ETA: %.1f min",
                     i, n, pct, rate, remaining / 60))
    }

    qdelta_results[[i]] <- compute_qdelta_per_obs(
      ps_est$propensity_scores[i, ], delta, divergence, tol = adaptive_tol
    )
  }

  # Compute EIF with active sets (key difference from standard estimator)
  message("Computing EIF with active set selection...")
  eif_result <- efficient_influence_function_active(
    A = A_mat,
    Y = data[[outcome_col]],
    propensity_scores = ps_est$propensity_scores,
    outcome_means = outcome_est$outcome_means,
    delta = delta,
    divergence = divergence,
    qdelta_results = qdelta_results,
    active_sets = active_sets
  )

  # Estimate and inference (same as standard estimator)
  estimate <- mean(eif_result$eif_values)
  variance <- var(eif_result$centered_eif) / n
  std_error <- sqrt(variance)

  # Confidence interval
  z_alpha <- qnorm(1 - (1 - confidence_level) / 2)
  ci_lower <- estimate - z_alpha * std_error
  ci_upper <- estimate + z_alpha * std_error

  # Bootstrap (if requested)
  if (bootstrap) {
    message(sprintf("Computing bootstrap CI (%d resamples)...", n_bootstrap))
    boot_estimates <- numeric(n_bootstrap)

    for (b in 1:n_bootstrap) {
      boot_idx <- sample(1:n, size = n, replace = TRUE)
      boot_eif <- eif_result$eif_values[boot_idx]
      boot_estimates[b] <- mean(boot_eif)
    }

    boot_ci_lower <- quantile(boot_estimates, probs = (1 - confidence_level) / 2)
    boot_ci_upper <- quantile(boot_estimates, probs = 1 - (1 - confidence_level) / 2)
    bootstrap_ci <- c(lower = boot_ci_lower, upper = boot_ci_upper)
  } else {
    bootstrap_ci <- NULL
  }

  # Compile results
  result <- list(
    estimate = estimate,
    std_error = std_error,
    se = std_error,  # Backward compatibility alias
    variance = variance,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    confidence_level = confidence_level,
    converged = eif_result$convergence_summary$convergence_rate >= 0.95,  # Top-level flag
    eif_values = eif_result$eif_values,
    centered_eif = eif_result$centered_eif,
    plugin = eif_result$plugin,
    convergence_info = eif_result$convergence_info,
    convergence_summary = eif_result$convergence_summary,
    bootstrap_ci = bootstrap_ci,
    delta = delta,
    divergence = divergence,
    estimator = estimator,
    n = n,
    m = m,
    # Vanishing threshold specific components
    threshold_info = list(
      epsilon_n = epsilon_n,
      gamma = gamma,
      min_active = min_active,
      eps_constant = eps_constant
    ),
    active_sets = active_sets,
    active_set_summary = active_summary,
    n_insufficient_active = n_insufficient
  )

  class(result) <- c("qdelta_vanishing_threshold", "list")

  return(result)
}


#' Estimate Vector Incremental Effects with Adaptive Delta
#'
#' Implements adaptive estimator with fixed threshold ε and observation-specific
#' delta_k(X). Uses active sets `K_ε(X) = {k: m_k(X) ∈ [ε, 1-ε]}` but adapts delta
#' based on propensity patterns to satisfy localized margin conditions.
#'
#' @inheritParams estimate_qdelta_mean
#' @param epsilon Numeric, fixed threshold for active set computation. Treatments
#'   with marginal propensities in `[ε, 1-ε]` are active. Typical values: 0.05-0.15.
#' @param delta Numeric vector (length m) of target delta values. These are adapted
#'   per observation based on propensity patterns.
#' @param scaling_method Character, how to scale delta for inactive treatments:
#'   \itemize{
#'     \item "proportional" (default): Scale proportionally to distance from boundary
#'     \item "fixed": Use fixed fraction of target delta
#'     \item "zero": Set to near-zero (exclude from constraint)
#'   }
#' @param inactive_fraction Numeric, fraction of target delta for inactive treatments
#'   when scaling_method = "fixed". Default: 0.5.
#' @param min_active Integer, minimum number of active treatments required per
#'   observation. Default: 1.
#'
#' @param max_qdelta_failure_rate Numeric in `[0, 1]`: the largest share of observations whose
#'   Q_delta solve may fail and fall back to Q_delta = P0 (no shift for that observation, a bias
#'   toward the null). Default `0`: the first failure stops the estimate with the observation index
#'   and the solver's message. A positive value tolerates failures up to `rate * n` observations
#'   (each warned about and counted), then stops. `1` always falls back. The vanishing-threshold
#'   estimator never falls back and always stops.
#' @param qdelta_tol Numeric or `NULL`. Solver tolerance for the per-observation Q_delta problem;
#'   `NULL` (default) uses `QDELTA_TOLERANCE` for every sample size (it no longer relaxes above
#'   100,000 observations, which made the same stratum's estimate depend on whether it crossed that
#'   size).
#' @return List with same structure as \code{\link{estimate_qdelta_mean}}, plus:
#' \describe{
#'   \item{epsilon}{The fixed threshold value}
#'   \item{delta_target}{The target delta values provided by user}
#'   \item{delta_adaptive}{Matrix (n x m) of observation-specific deltas}
#'   \item{scaling_method}{The delta scaling method used}
#'   \item{active_sets}{Active set information}
#'   \item{active_set_summary}{Summary of active sets}
#'   \item{n_insufficient_active}{Number of observations with insufficient active}
#' }
#'
#' @details
#' **Adaptive delta approach:**
#'
#' Unlike the vanishing threshold estimator (ε_n → 0, fixed δ), the adaptive
#' estimator uses:
#' - Fixed threshold ε (chosen by user)
#' - Adaptive δ_k(X) that varies by observation
#'
#' **Delta adaptation:**
#' For each observation:
#' - Active treatments (`m_k(X) ∈ [ε, 1-ε]`): Use target delta
#' - Inactive treatments: Scale delta based on distance from boundary
#'
#' This satisfies a **localized margin condition** (weaker than global condition
#' required by vanishing threshold).
#'
#' **Scaling methods:**
#' - **proportional**: δ_k(X) = δ_target * min(m_k/ε, (1-m_k)/ε, 1)
#'   Smooth scaling based on propensity distance from boundary
#' - **fixed**: δ_k(X) = δ_target * inactive_fraction for inactive treatments
#'   Simple fixed scaling
#' - **zero**: δ_k(X) ≈ 0 for inactive treatments (binary active/inactive)
#'
#' **Comparison to vanishing threshold:**
#' - Vanishing: ε_n → 0, needs global margin, asymptotically efficient
#' - Adaptive: Fixed ε, needs localized margin, potentially more stable with finite n
#'
#' **Choosing epsilon:**
#' - Smaller ε → more active treatments → closer to standard estimator
#' - Larger ε → fewer active treatments → more aggressive boundary handling
#' - Typical range: 0.05-0.15
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data with boundary propensities
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 3, p = 2, seed = 123)
#' tcols <- paste0("A", 1:3)
#' xcols <- paste0("X", 1:2)
#'
#' # Adaptive estimator with fixed epsilon = 0.1
#' result_adaptive <- estimate_qdelta_adaptive(
#'   data = sim_data$data,
#'   delta = c(0.1, 0.1, 0.1),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   epsilon = 0.1,
#'   divergence = "kl"
#' )
#'
#' # Check adaptive deltas
#' summary(result_adaptive$delta_adaptive)
#'
#' # Check active set diagnostics
#' print(result_adaptive$active_set_summary)
#'
#' # Compare to vanishing threshold
#' result_vt <- estimate_qdelta_vanishing_threshold(
#'   data = sim_data$data,
#'   delta = c(0.1, 0.1, 0.1),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   divergence = "kl"
#' )
#'
#' # Custom scaling method
#' result_adaptive_fixed <- estimate_qdelta_adaptive(
#'   data = sim_data$data,
#'   delta = c(0.1, 0.1, 0.1),
#'   treatment_cols = tcols,
#'   covariate_cols = xcols,
#'   epsilon = 0.1,
#'   divergence = "kl",
#'   scaling_method = "fixed",
#'   inactive_fraction = 0.3
#' )
#' }
estimate_qdelta_adaptive <- function(data, delta,
                                     treatment_cols, covariate_cols,
                                     outcome_col = "Y", epsilon = 0.1, divergence = "kl",
                                     ps_method = "glm", binary_method = "glm",
                                     outcome_method = "glm",
                                     estimator = "one_step",
                                     cross_fit = TRUE, n_folds = 5,
                                     confidence_level = 0.95,
                                     bootstrap = FALSE, n_bootstrap = 1000,
                                     scaling_method = "proportional",
                                     inactive_fraction = 0.5,
                                     min_active = 1,
    max_qdelta_failure_rate = 0, qdelta_tol = NULL) {

  # Input validation
  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a positive numeric vector")
  }

  if (!is.numeric(epsilon) || length(epsilon) != 1 || epsilon <= 0 || epsilon >= 0.5) {
    stop("epsilon must be a single number in (0, 0.5)")
  }

  if (!scaling_method %in% c("proportional", "fixed", "zero")) {
    stop("scaling_method must be 'proportional', 'fixed', or 'zero'")
  }

  # Explicit column specifications (no auto-detection).
  validate_columns(data, treatment_cols, covariate_cols, outcome_col)
  m <- length(treatment_cols)

  if (length(delta) != m) {
    stop("delta length must equal number of treatments")
  }

  message(sprintf("Adaptive estimator: fixed epsilon = %.3f, scaling_method = '%s'",
                 epsilon, scaling_method))

  # Estimate nuisance parameters
  ps_est <- estimate_propensity_scores(data,
                                       treatment_cols = treatment_cols,
                                       covariate_cols = covariate_cols,
                                       method = ps_method,
                                       binary_method = binary_method,
                                       cross_fit = cross_fit, n_folds = n_folds)
  gc(verbose = FALSE)

  outcome_est <- estimate_outcome_models(data,
                                         treatment_cols = treatment_cols,
                                         covariate_cols = covariate_cols,
                                         outcome_col = outcome_col,
                                         method = outcome_method,
                                         cross_fit = cross_fit, n_folds = n_folds)
  gc(verbose = FALSE)

  # Validate propensity scores
  if (any(!is.finite(ps_est$propensity_scores))) {
    stop("Propensity scores contain non-finite values")
  }

  ps_row_sums <- rowSums(ps_est$propensity_scores)
  if (!all(abs(ps_row_sums - 1) < PROPENSITY_SUM_TOL)) {
    max_deviation <- max(abs(ps_row_sums - 1))
    stop("Propensity scores do not sum to 1 (max deviation: ", sprintf("%.2e", max_deviation), ")")
  }

  # Validate outcome means
  if (any(!is.finite(outcome_est$outcome_means))) {
    stop("Outcome means contain non-finite values")
  }

  # Extract marginal propensities
  A_mat <- as.matrix(data[, treatment_cols, drop = FALSE])
  marginal_probs <- extract_marginal_probabilities(ps_est$propensity_scores, m = m)

  # Compute active sets with fixed epsilon
  n <- nrow(data)
  active_sets <- compute_active_sets(marginal_probs, epsilon = epsilon, min_active = min_active)

  # Summarize active sets
  active_summary <- summarize_active_sets(active_sets)
  message(sprintf("Active sets: %.1f%% obs have all %d treatments active, median = %d active",
                 100 * active_summary$pct_all_active, m, active_summary$median_n_active))

  # Warn if many observations have insufficient active treatments
  n_insufficient <- sum(!active_sets$sufficient_active)
  if (n_insufficient > 0) {
    pct_insufficient <- 100 * n_insufficient / n
    if (pct_insufficient > 10) {
      warning(sprintf("%.1f%% of observations have insufficient active treatments (< %d).\n",
                     pct_insufficient, min_active,
                     "Consider decreasing epsilon or min_active."))
    } else {
      message(sprintf("%.1f%% of observations have insufficient active treatments", pct_insufficient))
    }
  }

  # Compute adaptive deltas
  message("Computing observation-specific adaptive deltas...")
  delta_adaptive <- compute_adaptive_deltas(
    marginal_probs = marginal_probs,
    delta_target = delta,
    active_sets = active_sets,
    epsilon = epsilon,
    scaling_method = scaling_method,
    inactive_fraction = inactive_fraction
  )

  # Diagnostics: how much did deltas change?
  delta_changes <- delta_adaptive - matrix(delta, nrow = n, ncol = m, byrow = TRUE)
  mean_abs_change <- mean(abs(delta_changes))
  pct_unchanged <- 100 * mean(rowSums(abs(delta_changes) < 1e-8) == m)
  message(sprintf("Delta adaptation: %.1f%% obs unchanged, mean |change| = %.4f",
                 pct_unchanged, mean_abs_change))

  # Compute per-observation Q_delta with adaptive deltas
  qdelta_results <- vector("list", n)

  # One tolerance for every sample size (see qdelta_tol)
  adaptive_tol <- resolve_qdelta_tol(qdelta_tol)
  n_fallback <- 0L
  max_failure_rate <- resolve_max_failure_rate(max_qdelta_failure_rate)

  # Progress reporting
  report_interval <- if (n > QDELTA_LARGE_N_THRESHOLD) 10000 else Inf
  start_time <- Sys.time()

  for (i in 1:n) {
    if (i %% report_interval == 0) {
      pct <- 100 * i / n
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      rate <- i / elapsed
      remaining <- (n - i) / rate
      message(sprintf("Q_delta: %d / %d (%.1f%%) | %.1f obs/sec | ETA: %.1f min",
                     i, n, pct, rate, remaining / 60))
    }

    # Use observation-specific delta. A solver failure stops the estimate unless
    # more than max_qdelta_failure_rate of observations fail; otherwise the observation falls back to Q_delta = P0
    # (the intervention cannot move that observation: a bias toward the null that is counted
    # and warned about below).
    qdelta_results[[i]] <- tryCatch(
      compute_qdelta_per_obs(
        ps_est$propensity_scores[i, ],
        delta_adaptive[i, ],
        divergence = divergence,
        tol = adaptive_tol
      ),
      error = function(e) {
        n_fallback <<- n_fallback + 1L
        qdelta_fallback_or_stop(i, e, ps_est$propensity_scores[i, ], delta_adaptive[i, ],
                                n_fallback, n, max_failure_rate)
      }
    )
  }

  if (n_fallback > 0L) {
    warning(sprintf("Q_delta fallback to P0 used for %d of %d observations (max_qdelta_failure_rate = %g).",
                    n_fallback, n, max_failure_rate), call. = FALSE)
  }

  # Compute EIF with active sets
  message("Computing EIF with active set selection...")
  eif_result <- efficient_influence_function_active(
    A = A_mat,
    Y = data[[outcome_col]],
    propensity_scores = ps_est$propensity_scores,
    outcome_means = outcome_est$outcome_means,
    delta = delta,  # Use target delta for EIF (not adaptive)
    divergence = divergence,
    qdelta_results = qdelta_results,
    active_sets = active_sets
  )

  # Estimate and inference
  estimate <- mean(eif_result$eif_values)
  variance <- var(eif_result$centered_eif) / n
  std_error <- sqrt(variance)

  # Confidence interval
  z_alpha <- qnorm(1 - (1 - confidence_level) / 2)
  ci_lower <- estimate - z_alpha * std_error
  ci_upper <- estimate + z_alpha * std_error

  # Bootstrap (if requested)
  if (bootstrap) {
    message(sprintf("Computing bootstrap CI (%d resamples)...", n_bootstrap))
    boot_estimates <- numeric(n_bootstrap)

    for (b in 1:n_bootstrap) {
      boot_idx <- sample(1:n, size = n, replace = TRUE)
      boot_eif <- eif_result$eif_values[boot_idx]
      boot_estimates[b] <- mean(boot_eif)
    }

    boot_ci_lower <- quantile(boot_estimates, probs = (1 - confidence_level) / 2)
    boot_ci_upper <- quantile(boot_estimates, probs = 1 - (1 - confidence_level) / 2)
    bootstrap_ci <- c(lower = boot_ci_lower, upper = boot_ci_upper)
  } else {
    bootstrap_ci <- NULL
  }

  # Compile results
  result <- list(
    estimate = estimate,
    std_error = std_error,
    se = std_error,  # Backward-compat alias for std_error (matches estimate_qdelta_mean;
                     # consumers such as sims/study2/run_simulation.R read `se`).
    variance = variance,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    confidence_level = confidence_level,
    converged = eif_result$convergence_summary$convergence_rate >= 0.95,  # Top-level flag
    eif_values = eif_result$eif_values,
    centered_eif = eif_result$centered_eif,
    plugin = eif_result$plugin,
    convergence_info = eif_result$convergence_info,
    convergence_summary = eif_result$convergence_summary,
    bootstrap_ci = bootstrap_ci,
    delta_target = delta,
    delta_adaptive = delta_adaptive,
    divergence = divergence,
    estimator = estimator,
    n = n,
    m = m,
    # Adaptive estimator specific components
    epsilon = epsilon,
    scaling_method = scaling_method,
    inactive_fraction = if (scaling_method == "fixed") inactive_fraction else NULL,
    active_sets = active_sets,
    active_set_summary = active_summary,
    n_insufficient_active = n_insufficient,
    delta_adaptation_summary = list(
      mean_abs_change = mean_abs_change,
      pct_unchanged = pct_unchanged,
      min_delta = apply(delta_adaptive, 2, min),
      max_delta = apply(delta_adaptive, 2, max)
    )
  )

  class(result) <- c("qdelta_adaptive", "list")

  return(result)
}


#' Resolve the per-observation Q_delta solver tolerance (internal)
#'
#' @param qdelta_tol Numeric scalar, or `NULL` for the package default `QDELTA_TOLERANCE`.
#' @return Positive numeric scalar.
#' @keywords internal
resolve_qdelta_tol <- function(qdelta_tol = NULL) {
  if (is.null(qdelta_tol)) {
    return(QDELTA_TOLERANCE)
  }
  if (!is.numeric(qdelta_tol) || length(qdelta_tol) != 1L || !is.finite(qdelta_tol) ||
      qdelta_tol <= 0) {
    stop("qdelta_tol must be NULL or a single positive number")
  }
  qdelta_tol
}

#' Validate the tolerated Q_delta failure rate (internal)
#'
#' @param max_qdelta_failure_rate Numeric scalar in `[0, 1]`.
#' @return The validated scalar.
#' @keywords internal
resolve_max_failure_rate <- function(max_qdelta_failure_rate = 0) {
  if (!is.numeric(max_qdelta_failure_rate) || length(max_qdelta_failure_rate) != 1L ||
      !is.finite(max_qdelta_failure_rate) || max_qdelta_failure_rate < 0 ||
      max_qdelta_failure_rate > 1) {
    stop("max_qdelta_failure_rate must be a single number in [0, 1]")
  }
  max_qdelta_failure_rate
}

#' Handle a failed per-observation Q_delta solve (internal)
#'
#' Called once per failure, with `n_fallback` already counting this failure. If the failures so far
#' exceed `max_failure_rate * n` the run stops with a classed `qdelta_solver_failure` condition
#' (the default rate of 0 stops at the first failure); otherwise returns the degenerate solution
#' Q_delta = P0 (lambda = 0, the observation is left unshifted, a bias toward the null) marked
#' `converged = FALSE`, after a warning.
#'
#' @param i Observation index.
#' @param e The condition raised by the solver.
#' @param ps_row Numeric vector of length 2^m, the observation's propensity scores.
#' @param delta_i Numeric vector of length m, the delta used for this observation.
#' @param n_fallback Integer, failures so far including this one.
#' @param n Integer, observations in this fit.
#' @param max_failure_rate Numeric in `[0, 1]`, from `max_qdelta_failure_rate`.
#' @return A `compute_qdelta_per_obs()`-shaped list (only when the failure is tolerated).
#' @keywords internal
qdelta_fallback_or_stop <- function(i, e, ps_row, delta_i, n_fallback, n, max_failure_rate) {
  if (n_fallback > max_failure_rate * n) {
    msg <- if (max_failure_rate == 0) {
      sprintf(paste0("Q_delta solver failed for observation %d: %s\n",
                     "Set max_qdelta_failure_rate > 0 to tolerate failures (each falls back to ",
                     "Q_delta = P0, biasing the estimate toward the null)."),
              i, conditionMessage(e))
    } else {
      sprintf(paste0("Q_delta solver failures exceeded max_qdelta_failure_rate: %d of %d ",
                     "observations (%.2f%% > %.2f%%). Latest: observation %d: %s"),
              n_fallback, n, 100 * n_fallback / n, 100 * max_failure_rate, i, conditionMessage(e))
    }
    # Classed so estimate_pattern_effect() rethrows it instead of recording "this pattern failed"
    # and dropping the pattern from the stratified estimate.
    stop(structure(class = c("qdelta_solver_failure", "error", "condition"),
                   list(message = msg, call = NULL)))
  }
  warning(sprintf("Observation %d failed with error: %s. Using fallback (P0).",
                  i, conditionMessage(e)), call. = FALSE)
  list(
    q_delta = ps_row,
    lambda = rep(0, length(delta_i) + 1),
    converged = FALSE,
    constraint_norm = NA,
    iterations = NA,
    delta = delta_i
  )
}
