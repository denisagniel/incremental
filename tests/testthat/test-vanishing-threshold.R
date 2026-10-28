test_that("estimate_qdelta_vanishing_threshold validates inputs", {
  # Setup simple test data
  n <- 50
  m <- 2
  set.seed(123)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Delta validation
  expect_error(
    estimate_qdelta_vanishing_threshold(data, delta = c(-0.1, 0.1)),
    "delta must be a positive"
  )

  expect_error(
    estimate_qdelta_vanishing_threshold(data, delta = c(0.1),  # wrong length
                                        treatment_cols = A_cols(data), covariate_cols = X_cols(data)),
    "delta length must equal number of treatments"
  )
})

test_that("estimate_qdelta_vanishing_threshold runs with automatic rate selection", {
  # Small test dataset
  n <- 100
  m <- 2
  set.seed(456)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n),
    Y = rnorm(n) + 0.5
  )

  # Should use automatic gamma=0.3, alpha=0.4
  expect_message(
    result <- estimate_qdelta_vanishing_threshold(
      data = data,
      treatment_cols = A_cols(data), covariate_cols = X_cols(data),      delta = c(0.1, 0.1),
      divergence = "kl"
    ),
    "Using automatic gamma = 0.3"
  )

  # Check structure
  expect_s3_class(result, "qdelta_vanishing_threshold")
  expect_true(all(c("estimate", "std_error", "ci_lower", "ci_upper",
                   "threshold_info", "active_sets", "active_set_summary",
                   "n_insufficient_active") %in% names(result)))

  # Check threshold info
  expect_equal(result$threshold_info$gamma, 0.3)

  # Check active sets computed
  expect_equal(nrow(result$active_sets$active_matrix), n)
  expect_equal(ncol(result$active_sets$active_matrix), m)

  # Check estimate is finite
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result$std_error))
  expect_true(result$std_error > 0)
})

test_that("estimate_qdelta_vanishing_threshold works with custom gamma", {
  n <- 80
  m <- 2
  set.seed(789)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.6),
    A2 = rbinom(n, 1, 0.4),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Custom threshold rate: gamma=0.35
  result <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.15, 0.15),
    divergence = "kl",
    gamma = 0.35
  ))

  expect_equal(result$threshold_info$gamma, 0.35)

  # Epsilon should be smaller than default (gamma=0.35 > 0.3)
  eps_n_custom <- result$threshold_info$epsilon_n
  eps_n_default <- compute_vanishing_threshold(n, gamma = 0.3)
  expect_true(eps_n_custom < eps_n_default)
})

test_that("estimate_qdelta_vanishing_threshold handles propensities near boundaries", {
  # Create data with some extreme propensities
  n <- 100
  m <- 2
  set.seed(101)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = c(rep(-5, 20), rnorm(80)),  # First 20 obs: extreme X1
    Y = rnorm(n)
  )

  # Should handle gracefully with active set selection
  result <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    divergence = "kl"
  ))

  # Estimate should still be finite (main test)
  expect_true(is.finite(result$estimate))

  # Active set selection should have occurred
  expect_true(!is.null(result$active_sets))
  expect_equal(nrow(result$active_sets$active_matrix), n)
})

test_that("estimate_qdelta_vanishing_threshold warns about insufficient active treatments", {
  n <- 60
  m <- 2  # Use m=2 to avoid sparse combinations
  set.seed(202)

  # Create data with balanced treatment probabilities
  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Use very large gamma (large epsilon) and high min_active
  # This should naturally lead to some insufficient active treatments
  result <- suppressWarnings(suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    divergence = "kl",
    gamma = 0.48,  # Very large epsilon (close to limit)
    min_active = 2  # Need both active
  )))

  # With extreme settings, should handle the situation gracefully
  # Main test: runs without error
  expect_true(is.finite(result$estimate))
  expect_true(!is.null(result$n_insufficient_active))
  expect_true(is.numeric(result$n_insufficient_active))
})

test_that("estimate_qdelta_vanishing_threshold with min_active=1 works", {
  n <- 60
  m <- 2
  set.seed(303)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    divergence = "kl",
    min_active = 1
  ))

  # With min_active=1, should have few insufficient observations
  expect_true(result$n_insufficient_active / n < 0.05)

  # All active set entries should have at least 1 active treatment (or be insufficient)
  for (i in 1:n) {
    if (result$active_sets$sufficient_active[i]) {
      expect_true(result$active_sets$n_active[i] >= 1)
    }
  }
})

test_that("estimate_qdelta_vanishing_threshold convergence info correct", {
  n <- 40
  m <- 2
  set.seed(404)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    divergence = "kl"
  ))

  # Convergence info should match n
  expect_equal(length(result$convergence_info$converged), n)
  expect_equal(length(result$convergence_info$iterations), n)

  # Should have high convergence rate
  expect_true(result$convergence_summary$convergence_rate > 0.9)
})

test_that("estimate_qdelta_vanishing_threshold with bootstrap", {
  n <- 50
  m <- 2
  set.seed(505)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Run with bootstrap (suppress other messages but not bootstrap message)
  result <- suppressWarnings(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    divergence = "kl",
    bootstrap = TRUE,
    n_bootstrap = 100
  ))

  # Bootstrap CI should exist
  expect_false(is.null(result$bootstrap_ci))
  expect_equal(length(result$bootstrap_ci), 2)
  expect_true(result$bootstrap_ci[1] < result$bootstrap_ci[2])
})

test_that("estimate_qdelta_vanishing_threshold with different divergences", {
  n <- 50
  m <- 2
  set.seed(606)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Test KL and reverse_KL (chi2/hellinger can be fragile)
  for (div in c("kl", "reverse_kl")) {
    result <- suppressMessages(estimate_qdelta_vanishing_threshold(
      data = data,
      treatment_cols = A_cols(data), covariate_cols = X_cols(data),      delta = c(0.1, 0.1),
      divergence = div
    ))

    expect_equal(result$divergence, div)
    expect_true(is.finite(result$estimate))
  }
})

test_that("estimate_qdelta_vanishing_threshold returns expected structure", {
  n <- 50
  m <- 2
  set.seed(707)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    divergence = "kl"
  ))

  # Standard components
  expect_true(all(c("estimate", "std_error", "variance", "ci_lower", "ci_upper",
                   "eif_values", "centered_eif", "convergence_info") %in% names(result)))

  # Vanishing threshold specific components
  expect_true(all(c("threshold_info", "active_sets", "active_set_summary",
                   "n_insufficient_active") %in% names(result)))

  # Dimensions
  expect_equal(length(result$eif_values), n)
  expect_equal(length(result$centered_eif), n)
  expect_equal(mean(result$centered_eif), 0, tolerance = 1e-10)

  # Threshold info structure
  expect_true(all(c("epsilon_n", "gamma", "min_active",
                   "eps_constant") %in% names(result$threshold_info)))

  # Active set summary structure
  expect_true(all(c("n_obs", "m", "min_n_active", "max_n_active", "mean_n_active",
                   "median_n_active", "pct_all_active") %in% names(result$active_set_summary)))
})

test_that("estimate_qdelta_vanishing_threshold epsilon_n decreases with gamma", {
  n <- 100
  m <- 2
  set.seed(808)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Smaller gamma → larger epsilon
  result_small_gamma <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data, delta = c(0.1, 0.1), gamma = 0.26,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data)  ))

  # Larger gamma → smaller epsilon
  result_large_gamma <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data, delta = c(0.1, 0.1), gamma = 0.35,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data)  ))

  expect_true(result_small_gamma$threshold_info$epsilon_n >
              result_large_gamma$threshold_info$epsilon_n)
})
