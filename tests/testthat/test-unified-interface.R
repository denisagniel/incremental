test_that("estimate_qdelta_mean with threshold_type='none' works (backward compatible)", {
  # Test that default behavior unchanged
  n <- 100
  m <- 2
  set.seed(123)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n),
    Y = rnorm(n) + 0.5
  )

  # Default should be threshold_type = "none"
  result <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    divergence = "kl"
  ))

  # Check it's the standard class (not vanishing or adaptive)
  expect_false(inherits(result, "qdelta_vanishing_threshold"))
  expect_false(inherits(result, "qdelta_adaptive"))

  # Should have standard components but not threshold-specific ones
  expect_true("estimate" %in% names(result))
  expect_true("std_error" %in% names(result))
  expect_false("threshold_info" %in% names(result))  # Vanishing-specific
  expect_false("delta_adaptive" %in% names(result))  # Adaptive-specific

  # Estimate should be finite
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result$std_error))
})

test_that("estimate_qdelta_mean with threshold_type='vanishing' routes correctly", {
  n <- 100
  m <- 2
  set.seed(456)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    divergence = "kl",
    threshold_type = "vanishing"
  ))

  # Should be vanishing threshold class
  expect_s3_class(result, "qdelta_vanishing_threshold")

  # Should have vanishing-specific components
  expect_true("threshold_info" %in% names(result))
  expect_true("active_sets" %in% names(result))
  expect_true("active_set_summary" %in% names(result))

  # threshold_info should have expected fields
  expect_true(all(c("epsilon_n", "gamma") %in%
                  names(result$threshold_info)))

  # Estimate should be finite
  expect_true(is.finite(result$estimate))
})

test_that("estimate_qdelta_mean with threshold_type='adaptive' routes correctly", {
  n <- 100
  m <- 2
  set.seed(789)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    divergence = "kl",
    threshold_type = "adaptive",
    epsilon = 0.1
  ))

  # Should be adaptive class
  expect_s3_class(result, "qdelta_adaptive")

  # Should have adaptive-specific components
  expect_true("epsilon" %in% names(result))
  expect_true("delta_adaptive" %in% names(result))
  expect_true("scaling_method" %in% names(result))
  expect_true("delta_adaptation_summary" %in% names(result))

  # epsilon should match input
  expect_equal(result$epsilon, 0.1)

  # Estimate should be finite
  expect_true(is.finite(result$estimate))
})

test_that("estimate_qdelta_mean passes parameters correctly to vanishing threshold", {
  n <- 80
  m <- 2
  set.seed(101)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Custom gamma passes through to the vanishing-threshold estimator
  result <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.15, 0.15),
    divergence = "kl",
    threshold_type = "vanishing",
    gamma = 0.35
  ))

  # Check parameter passed through
  expect_equal(result$threshold_info$gamma, 0.35)
})

test_that("estimate_qdelta_mean passes parameters correctly to adaptive", {
  n <- 80
  m <- 2
  set.seed(202)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Custom epsilon and scaling method
  result <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.15, 0.15),
    divergence = "kl",
    threshold_type = "adaptive",
    epsilon = 0.15,
    scaling_method = "fixed",
    inactive_fraction = 0.3
  ))

  # Check parameters passed through
  expect_equal(result$epsilon, 0.15)
  expect_equal(result$scaling_method, "fixed")
  expect_equal(result$inactive_fraction, 0.3)
})

test_that("estimate_qdelta_mean threshold_type argument matching works", {
  n <- 50
  m <- 2
  set.seed(303)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Partial matching should work
  result_v <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    threshold_type = "van"  # Partial match
  ))
  expect_s3_class(result_v, "qdelta_vanishing_threshold")

  result_a <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    threshold_type = "adap"  # Partial match
  ))
  expect_s3_class(result_a, "qdelta_adaptive")

  # Invalid should error
  expect_error(
    estimate_qdelta_mean(
      data = data,
      treatment_cols = A_cols(data), covariate_cols = X_cols(data),
      delta = c(0.1, 0.1),
      threshold_type = "invalid"
    ),
    "should be one of"
  )
})

test_that("estimate_qdelta_mean with different threshold types produces similar estimates", {
  # With well-behaved data, all three should give similar results
  n <- 200
  m <- 2
  set.seed(404)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n),
    Y = rnorm(n) + 0.2 * data.frame(A1 = rbinom(n, 1, 0.5))$A1
  )

  # Standard
  result_none <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    divergence = "kl",
    threshold_type = "none"
  ))

  # Vanishing threshold
  result_vt <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    divergence = "kl",
    threshold_type = "vanishing"
  ))

  # Adaptive
  result_adaptive <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    divergence = "kl",
    threshold_type = "adaptive",
    epsilon = 0.1
  ))

  # All should be finite
  expect_true(is.finite(result_none$estimate))
  expect_true(is.finite(result_vt$estimate))
  expect_true(is.finite(result_adaptive$estimate))

  # With well-behaved data, estimates should be reasonably close
  # (Not exact due to different handling, but within ~2 SEs)
  mean_est <- mean(c(result_none$estimate, result_vt$estimate, result_adaptive$estimate))
  max_se <- max(result_none$std_error, result_vt$std_error, result_adaptive$std_error)

  expect_true(abs(result_none$estimate - mean_est) < 2 * max_se)
  expect_true(abs(result_vt$estimate - mean_est) < 2 * max_se)
  expect_true(abs(result_adaptive$estimate - mean_est) < 2 * max_se)
})

test_that("estimate_qdelta_mean with bootstrap and threshold types", {
  n <- 50
  m <- 2
  set.seed(505)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Vanishing with bootstrap
  result_vt_boot <- suppressWarnings(suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    threshold_type = "vanishing",
    bootstrap = TRUE,
    n_bootstrap = 50  # Small for speed
  )))

  expect_false(is.null(result_vt_boot$bootstrap_ci))

  # Adaptive with bootstrap
  result_adaptive_boot <- suppressWarnings(suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    threshold_type = "adaptive",
    epsilon = 0.1,
    bootstrap = TRUE,
    n_bootstrap = 50
  )))

  expect_false(is.null(result_adaptive_boot$bootstrap_ci))
})

test_that("estimate_qdelta_mean unified interface handles cross-fitting", {
  n <- 100
  m <- 2
  set.seed(606)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n),
    Y = rnorm(n)
  )

  # Vanishing with cross-fitting
  result_vt_cf <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    threshold_type = "vanishing",
    cross_fit = TRUE,
    n_folds = 3
  ))

  expect_true(is.finite(result_vt_cf$estimate))

  # Adaptive with cross-fitting
  result_adaptive_cf <- suppressMessages(estimate_qdelta_mean(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),
    delta =c(0.1, 0.1),
    threshold_type = "adaptive",
    epsilon = 0.1,
    cross_fit = TRUE,
    n_folds = 3
  ))

  expect_true(is.finite(result_adaptive_cf$estimate))
})

test_that("estimate_qdelta_mean documentation examples work", {
  # Test examples from documentation
  set.seed(123)

  # Skip if simulate function not available
  skip_if_not(exists("simulate_vector_treatment_data"))

  sim_data <- simulate_vector_treatment_data(n = 100, m = 2, p = 3)

  # Basic estimation
  result <- suppressMessages(estimate_qdelta_mean(
    data = sim_data$data,
    treatment_cols = A_cols(sim_data$data), covariate_cols = X_cols(sim_data$data),
    delta = c(1.5, 2.0),
    divergence = "kl"
  ))
  expect_true(is.finite(result$estimate))

  # Vanishing threshold
  result_vt <- suppressMessages(estimate_qdelta_mean(
    data = sim_data$data,
    treatment_cols = A_cols(sim_data$data), covariate_cols = X_cols(sim_data$data),
    delta = c(1.5, 2.0),
    threshold_type = "vanishing"
  ))
  expect_s3_class(result_vt, "qdelta_vanishing_threshold")

  # Adaptive
  result_adaptive <- suppressMessages(estimate_qdelta_mean(
    data = sim_data$data,
    treatment_cols = A_cols(sim_data$data), covariate_cols = X_cols(sim_data$data),
    delta = c(1.5, 2.0),
    threshold_type = "adaptive",
    epsilon = 0.1
  ))
  expect_s3_class(result_adaptive, "qdelta_adaptive")

  # All should produce valid estimates
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result_vt$estimate))
  expect_true(is.finite(result_adaptive$estimate))
})
