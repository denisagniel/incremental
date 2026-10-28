#' @file data_generation.R
#' @title Data Generation for Vector Treatments
#' @description
#' Simulation utilities for generating synthetic data with multiple binary
#' treatments, covariates, and outcomes. Provides true propensity scores and
#' outcome means for validation and testing.
#' @author Denis Agniel

#' Simulate vector treatment data
#'
#' Generates synthetic data with multiple binary treatments, covariates,
#' and outcomes following a specified data generating process.
#'
#' @param n Integer, sample size
#' @param m Integer, number of binary treatments
#' @param p Integer, number of covariates
#' @param propensity_formula Character or function, propensity score model specification
#' @param outcome_formula Character or function, outcome model specification
#' @param treatment_effects Numeric vector, true treatment effects
#' @param interaction_effects Numeric matrix, treatment interaction effects
#' @param covariate_effects Numeric vector, covariate effects on outcome
#' @param propensity_coef Numeric vector of length 3, coefficients for propensity model (coef for X1, X2, previous treatment)
#' @param noise_sd Numeric, standard deviation of outcome noise
#' @param missing_prob Numeric, probability of missing outcomes (0 = no missing)
#' @param seed Integer, random seed for reproducibility
#' @return List containing:
#'   - data: tibble with columns X1, ..., Xp, A1, ..., Am, Y, R (missing indicator)
#'   - true_propensity: matrix (n x 2^m) of true propensity scores for each treatment combination.
#'     Columns are named with format "A1=x_A2=y_..." indicating which treatment combination each column represents.
#'   - true_outcome: matrix (n x 2^m) of true outcome means for each treatment combination.
#'     Columns are named with format "A1=x_A2=y_..." indicating which treatment combination each column represents.
#'   - dgp_params: parameters used in data generation
#' @export
#' @examples
#' # Simple example with 2 treatments, 2 covariates
#' sim_data <- simulate_vector_treatment_data(
#'   n = 1000, m = 2, p = 2,
#'   treatment_effects = c(0.5, 0.3),
#'   interaction_effects = matrix(c(0, 0.2, 0.2, 0), nrow = 2),
#'   seed = 123
#' )
#' head(sim_data$data)
simulate_vector_treatment_data <- function(n, m, p = 2,
                                          propensity_formula = "linear",
                                          outcome_formula = "linear",
                                          treatment_effects = rep(0.5, m),
                                          interaction_effects = NULL,
                                          covariate_effects = rep(0.3, p),
                                          propensity_coef = c(0.5, 0.3, 0.2),
                                          noise_sd = 1.0,
                                          missing_prob = 0.0,
                                          seed = NULL) {

  if (!is.null(seed)) withr::local_seed(seed)

  # Input validation
  if (length(treatment_effects) != m) {
    stop("treatment_effects length must equal m")
  }
  if (length(covariate_effects) != p) {
    stop("covariate_effects length must equal p")
  }
  if (is.null(interaction_effects)) {
    interaction_effects <- matrix(0, nrow = m, ncol = m)
  }
  if (nrow(interaction_effects) != m || ncol(interaction_effects) != m) {
    stop("interaction_effects must be m x m matrix")
  }

  # Generate covariates
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)

  # Generate treatment assignments
  A <- matrix(0, nrow = n, ncol = m)
  colnames(A) <- paste0("A", 1:m)

  # Simple propensity score model: logistic for each treatment (vectorized)
  for (k in 1:m) {
    # Linear combination of all covariates with some treatment dependence
    # Use all p covariates (not just first 2)
    if (p == 1) {
      covariate_term <- propensity_coef[1] * X[, 1]
    } else {
      # Use first min(p, 3) propensity coefficients for covariates
      n_coef_use <- min(p, length(propensity_coef))
      covariate_term <- rowSums(sweep(X[, 1:n_coef_use, drop = FALSE], 2, propensity_coef[1:n_coef_use], "*"))
    }
    treatment_dependence <- if (k > 1 && length(propensity_coef) > p) {
      propensity_coef[p + 1] * A[, k-1]
    } else {
      0
    }
    linear_pred <- covariate_term + treatment_dependence
    prob <- plogis(linear_pred)
    A[, k] <- rbinom(n, 1, prob)
  }
  # Note: This loop is kept as-is since each treatment depends on previous ones

  # Generate outcomes
  Y <- numeric(n)

  # Base outcome (intercept + covariate effects)
  base_outcome <- 1.0 + X %*% covariate_effects

  # Add treatment effects (vectorized)
  Y <- base_outcome + A %*% treatment_effects

  # Add interaction effects (vectorized)
  # Compute all pairwise interactions at once
  for (i in 1:m) {
    for (j in 1:m) {
      if (i != j && interaction_effects[i, j] != 0) {
        Y <- Y + interaction_effects[i, j] * A[, i] * A[, j]
      }
    }
  }
  # Note: Interaction loop kept as conditional logic prevents full vectorization

  # Add noise
  Y <- Y + rnorm(n, 0, noise_sd)

  # Generate missing outcome indicator
  R <- rbinom(n, 1, 1 - missing_prob)

  # Create data frame
  data <- tibble::tibble(
    X = X,
    A = A,
    Y = Y,
    R = R
  )

  # Flatten for easier use
  data_flat <- tibble::tibble(id = 1:n)
  for (i in 1:p) {
    data_flat[[paste0("X", i)]] <- X[, i]
  }
  for (i in 1:m) {
    data_flat[[paste0("A", i)]] <- A[, i]
  }
  data_flat$Y <- Y
  data_flat$R <- R

  # Compute true propensity scores and outcome means
  vectors <- generate_binary_vectors(m)
  true_propensity <- matrix(0, nrow = n, ncol = 2^m)
  true_outcome <- matrix(0, nrow = n, ncol = 2^m)
  
  # Create column names based on treatment combinations
  # Each column corresponds to a treatment combination (row of vectors)
  column_names <- character(2^m)
  for (j in 1:nrow(vectors)) {
    # Create a readable name like "A1=0_A2=0" or "A1=1_A2=1"
    treatment_str <- paste0("A", 1:m, "=", vectors[j, ], collapse = "_")
    column_names[j] <- treatment_str
  }
  
  # Add column names to matrices
  colnames(true_propensity) <- column_names
  colnames(true_outcome) <- column_names

  for (i in 1:n) {
    for (j in 1:nrow(vectors)) {
      # True propensity for this treatment combination
      prob <- 1.0
      for (k in 1:m) {
        # Use same covariate logic as data generation (lines 82-88)
        if (p == 1) {
          covariate_term <- propensity_coef[1] * X[i, 1]
        } else {
          # Use first min(p, length(propensity_coef)) covariates
          n_coef_use <- min(p, length(propensity_coef))
          covariate_term <- sum(X[i, 1:n_coef_use] * propensity_coef[1:n_coef_use])
        }
        # Treatment dependence with bounds checking
        treatment_dependence <- if (k > 1 && length(propensity_coef) > p) {
          propensity_coef[p + 1] * vectors[j, k-1]
        } else {
          0
        }
        linear_pred <- covariate_term + treatment_dependence
        prob <- prob * if (vectors[j, k] == 1) plogis(linear_pred) else (1 - plogis(linear_pred))
      }
      true_propensity[i, j] <- prob

      # True outcome mean for this treatment combination (partially vectorized)
      outcome_mean <- 1.0 + sum(X[i, ] * covariate_effects) + 
                      sum(treatment_effects * vectors[j, ])
      # Add interaction effects (vectorized where possible)
      interaction_sum <- 0
      for (ii in 1:m) {
        for (jj in 1:m) {
          if (ii != jj && interaction_effects[ii, jj] != 0) {
            interaction_sum <- interaction_sum + interaction_effects[ii, jj] * vectors[j, ii] * vectors[j, jj]
          }
        }
      }
      true_outcome[i, j] <- outcome_mean + interaction_sum
    }
  }

  return(list(
    data = data_flat,
    true_propensity = true_propensity,
    true_outcome = true_outcome,
    dgp_params = list(
      n = n, m = m, p = p,
      treatment_effects = treatment_effects,
      interaction_effects = interaction_effects,
      covariate_effects = covariate_effects,
      noise_sd = noise_sd,
      missing_prob = missing_prob
    )
  ))
}

#' Generate binary vectors for treatment combinations
#'

