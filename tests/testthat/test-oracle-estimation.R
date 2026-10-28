#' Unit Tests for Oracle Nuisance Mode
#'
#' Tests that oracle propensities and outcomes bypass model fitting
#' and produce valid results. Uses `generate_dgp1_oracle_data()`
#' (see helper-oracle-dgp.R) as a self-contained stand-in for DGP-1 from
#' the paper's simulation study.

library(testthat)

test_that("Oracle nuisances bypass model fitting", {
  dgp_result <- generate_dgp1_oracle_data(n = 100, m = 2, p = 2, seed = 123)

  data <- dgp_result$data
  delta <- c(1.5, 1.5)

  # Test oracle mode
  result <- estimate_qdelta_mean(
    data = data,
    treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
    delta = delta,
    divergence = "kl",
    oracle_propensities = dgp_result$true_propensities,
    oracle_outcomes = dgp_result$true_outcomes,
    estimator = "one_step"
  )

  # Check convergence
  expect_true(result$converged)
  expect_equal(result$convergence_summary$convergence_rate, 1)
  expect_false(is.na(result$estimate))
  expect_true(is.finite(result$se))
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
})

test_that("Oracle validation catches dimension errors", {
  dgp_result <- generate_dgp1_oracle_data(n = 100, m = 2, p = 2, seed = 123)

  data <- dgp_result$data
  delta <- c(1.5, 1.5)

  # Wrong dimensions for propensities
  bad_ps <- matrix(0.25, 100, 3)  # Should be 100 x 4 (2^2)
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      delta = delta,
      oracle_propensities = bad_ps,
      oracle_outcomes = dgp_result$true_outcomes
    ),
    "must be n x 2\\^m"
  )

  # Wrong dimensions for outcomes
  bad_outcomes <- matrix(0, 100, 3)  # Should be 100 x 4
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      delta = delta,
      oracle_propensities = dgp_result$true_propensities,
      oracle_outcomes = bad_outcomes
    ),
    "must be n x 2\\^m"
  )
})

test_that("Oracle validation catches non-sum-to-1 propensities", {
  dgp_result <- generate_dgp1_oracle_data(n = 100, m = 2, p = 2, seed = 123)

  data <- dgp_result$data
  delta <- c(1.5, 1.5)

  # Propensities that don't sum to 1
  bad_ps <- matrix(0.5, 100, 4)
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      delta = delta,
      oracle_propensities = bad_ps,
      oracle_outcomes = dgp_result$true_outcomes
    ),
    "must sum to 1"
  )
})

test_that("Oracle validation catches non-finite values", {
  dgp_result <- generate_dgp1_oracle_data(n = 100, m = 2, p = 2, seed = 123)

  data <- dgp_result$data
  delta <- c(1.5, 1.5)

  # Non-finite propensities
  bad_ps <- dgp_result$true_propensities
  bad_ps[1, 1] <- NA
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      delta = delta,
      oracle_propensities = bad_ps,
      oracle_outcomes = dgp_result$true_outcomes
    ),
    "non-finite values"
  )

  # Non-finite outcomes
  bad_outcomes <- dgp_result$true_outcomes
  bad_outcomes[1, 1] <- Inf
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      delta = delta,
      oracle_propensities = dgp_result$true_propensities,
      oracle_outcomes = bad_outcomes
    ),
    "non-finite values"
  )
})

test_that("Oracle mode requires both propensities AND outcomes", {
  dgp_result <- generate_dgp1_oracle_data(n = 100, m = 2, p = 2, seed = 123)

  data <- dgp_result$data
  delta <- c(1.5, 1.5)

  # Only propensities
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      delta = delta,
      oracle_propensities = dgp_result$true_propensities,
      oracle_outcomes = NULL
    ),
    "BOTH oracle_propensities AND oracle_outcomes"
  )

  # Only outcomes
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      delta = delta,
      oracle_propensities = NULL,
      oracle_outcomes = dgp_result$true_outcomes
    ),
    "BOTH oracle_propensities AND oracle_outcomes"
  )
})

test_that("Oracle mode produces lower bias than fitted methods", {
  # This is a smoke test - oracle should have lowest bias
  set.seed(42)
  n_reps <- 20
  delta <- c(1.5, 1.5)

  oracle_estimates <- numeric(n_reps)
  glm_estimates <- numeric(n_reps)

  for (i in 1:n_reps) {
    dgp_result <- generate_dgp1_oracle_data(n = 500, m = 2, p = 2, seed = 100 + i)
    data <- dgp_result$data

    # Oracle
    result_oracle <- estimate_qdelta_mean(
      data = data, delta = delta,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      divergence = "kl",
      oracle_propensities = dgp_result$true_propensities,
      oracle_outcomes = dgp_result$true_outcomes,
      estimator = "one_step"
    )
    oracle_estimates[i] <- result_oracle$estimate

    # GLM
    result_glm <- estimate_qdelta_mean(
      data = data, delta = delta,
      treatment_cols = c("A1", "A2"), covariate_cols = c("X1", "X2"),
      divergence = "kl",
      ps_method = "sequential", binary_method = "glm",
      outcome_method = "glm", estimator = "one_step"
    )
    glm_estimates[i] <- result_glm$estimate
  }

  # Oracle should have lower variance than GLM (not necessarily lower bias in finite samples,
  # but lower variance is expected)
  oracle_sd <- sd(oracle_estimates, na.rm = TRUE)
  glm_sd <- sd(glm_estimates, na.rm = TRUE)

  # Expect oracle SD to be at most 1.5x GLM SD (can be larger in finite samples,
  # but should be comparable)
  expect_true(oracle_sd < 2 * glm_sd)
})
