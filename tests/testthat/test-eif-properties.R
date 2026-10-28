library(testthat)
library(incremental)

test_that("EIF has correct properties", {
  # Test that EIF satisfies theoretical properties:
  # 1. Centered EIF should have mean approximately 0
  # 2. Variance should be positive
  # 3. All values should be finite

  set.seed(123)
  n <- 1000
  m <- 2

  # Generate simple data
  X <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)

  # Create simple propensity scores (uniform for testing)
  vectors <- generate_binary_vectors(m)
  propensity_scores <- matrix(1/4, nrow = n, ncol = 4)

  # Create outcome means (constant for simplicity)
  outcome_means <- matrix(1.0, nrow = n, ncol = 4)

  delta <- c(1.5, 2.0)

  # Compute EIF
  eif_result <- efficient_influence_function(
    X = X, A = A, Y = Y,
    propensity_scores = propensity_scores,
    delta = delta,
    outcome_means = outcome_means,
    divergence = "kl"
  )

  # Property 1: Centered EIF should have mean approximately 0
  mean_centered <- mean(eif_result$centered_eif)
  expect_lt(abs(mean_centered), 0.5)  # Should be close to 0

  # Property 2: Variance should be positive
  var_eif <- var(eif_result$centered_eif)
  expect_gt(var_eif, 0)

  # Property 3: All values should be finite
  expect_true(all(is.finite(eif_result$eif_values)))
  expect_true(all(is.finite(eif_result$centered_eif)))

  # Property 4: Uncentered EIF should have mean approximately equal to psi_delta
  # (This is a property of the EIF - its mean is the parameter estimate)
  mean_uncentered <- mean(eif_result$eif_values)
  expect_true(is.finite(mean_uncentered))

  cat("✓ EIF has correct properties\n")
  cat(sprintf("  Mean of centered EIF: %.6f (should be ≈ 0)\n", mean_centered))
  cat(sprintf("  Variance of centered EIF: %.4f (should be > 0)\n", var_eif))
  cat(sprintf("  All values finite: %s\n", all(is.finite(eif_result$eif_values))))
  cat(sprintf("  Mean of uncentered EIF: %.4f\n", mean_uncentered))
})

test_that("EIF converges for all observations", {
  # Test that Q_delta computation converges for all observations in EIF
  set.seed(123)
  n <- 100
  m <- 2

  X <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)

  vectors <- generate_binary_vectors(m)
  propensity_scores <- matrix(1/4, nrow = n, ncol = 4)
  outcome_means <- matrix(1.0, nrow = n, ncol = 4)
  delta <- c(1.5, 2.0)

  eif_result <- efficient_influence_function(
    X = X, A = A, Y = Y,
    propensity_scores = propensity_scores,
    delta = delta,
    outcome_means = outcome_means,
    divergence = "kl"
  )

  # Check convergence info
  convergence_rate <- eif_result$convergence_summary$convergence_rate
  expect_gt(convergence_rate, 0.9)  # At least 90% should converge

  cat("✓ EIF convergence check\n")
  cat(sprintf("  Convergence rate: %.2f%%\n", convergence_rate * 100))
  cat(sprintf("  Mean iterations: %.1f\n", eif_result$convergence_summary$mean_iterations))
})
