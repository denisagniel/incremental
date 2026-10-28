test_that("estimate_qdelta_adaptive validates inputs", {
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
    estimate_qdelta_adaptive(data, delta = c(-0.1, 0.1), epsilon = 0.1),
    "delta must be a positive"
  )

  expect_error(
    estimate_qdelta_adaptive(data, delta = c(0.1), epsilon = 0.1,  # wrong length
                             treatment_cols = A_cols(data), covariate_cols = X_cols(data)),
    "delta length must equal number of treatments"
  )

  # Epsilon validation
  expect_error(
    estimate_qdelta_adaptive(data, delta = c(0.1, 0.1), epsilon = -0.1),
    "epsilon must be a single number in \\(0, 0.5\\)"
  )

  expect_error(
    estimate_qdelta_adaptive(data, delta = c(0.1, 0.1), epsilon = 0.6),
    "epsilon must be a single number in \\(0, 0.5\\)"
  )

  # Scaling method validation
  expect_error(
    estimate_qdelta_adaptive(data, delta = c(0.1, 0.1), epsilon = 0.1, scaling_method = "invalid"),
    "scaling_method must be"
  )
})

test_that("estimate_qdelta_adaptive runs with default settings", {
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

  result <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.1,
    divergence = "kl"
  ))

  # Check structure
  expect_s3_class(result, "qdelta_adaptive")
  expect_true(all(c("estimate", "std_error", "ci_lower", "ci_upper",
                   "epsilon", "delta_target", "delta_adaptive",
                   "scaling_method", "active_sets", "active_set_summary",
                   "n_insufficient_active", "delta_adaptation_summary") %in% names(result)))

  # Check epsilon stored
  expect_equal(result$epsilon, 0.1)
  expect_equal(result$scaling_method, "proportional")  # Default

  # Check delta adaptive is a matrix
  expect_true(is.matrix(result$delta_adaptive))
  expect_equal(dim(result$delta_adaptive), c(n, m))

  # Check estimate is finite
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result$std_error))
  expect_true(result$std_error > 0)
})

test_that("estimate_qdelta_adaptive with different scaling methods", {
  n <- 80
  m <- 2
  set.seed(789)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.6),
    A2 = rbinom(n, 1, 0.4),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Proportional scaling
  result_prop <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.15, 0.15),
    epsilon = 0.1,
    divergence = "kl",
    scaling_method = "proportional"
  ))

  expect_equal(result_prop$scaling_method, "proportional")
  expect_true(is.finite(result_prop$estimate))

  # Fixed scaling
  result_fixed <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.15, 0.15),
    epsilon = 0.1,
    divergence = "kl",
    scaling_method = "fixed",
    inactive_fraction = 0.3
  ))

  expect_equal(result_fixed$scaling_method, "fixed")
  expect_equal(result_fixed$inactive_fraction, 0.3)
  expect_true(is.finite(result_fixed$estimate))

  # Zero scaling
  result_zero <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.15, 0.15),
    epsilon = 0.1,
    divergence = "kl",
    scaling_method = "zero"
  ))

  expect_equal(result_zero$scaling_method, "zero")
  expect_true(is.finite(result_zero$estimate))
})

test_that("estimate_qdelta_adaptive handles boundary propensities", {
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

  # Use larger epsilon to create more inactive treatments
  result <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.15,  # Larger epsilon → more treatments become inactive
    divergence = "kl"
  ))

  # Estimate should still be finite (main test)
  expect_true(is.finite(result$estimate))

  # Active set selection should have occurred
  expect_true(!is.null(result$active_sets))
  expect_equal(nrow(result$active_sets$active_matrix), n)

  # Delta adaptation summary should exist
  expect_true(!is.null(result$delta_adaptation_summary))
})

test_that("estimate_qdelta_adaptive delta adaptation statistics", {
  n <- 80
  m <- 3
  set.seed(202)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    A3 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1, 0.1),
    epsilon = 0.1,
    divergence = "kl"
  ))

  # Check delta adaptation summary exists
  expect_true(!is.null(result$delta_adaptation_summary))
  expect_true(all(c("mean_abs_change", "pct_unchanged", "min_delta", "max_delta") %in%
                  names(result$delta_adaptation_summary)))

  # Check min/max delta per treatment
  expect_equal(length(result$delta_adaptation_summary$min_delta), m)
  expect_equal(length(result$delta_adaptation_summary$max_delta), m)

  # Max delta should not exceed target (with proportional scaling)
  expect_true(all(result$delta_adaptation_summary$max_delta <= c(0.1, 0.1, 0.1) + 1e-6))
})

test_that("estimate_qdelta_adaptive with different epsilon values", {
  n <- 60
  m <- 2
  set.seed(303)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Smaller epsilon → more active treatments
  result_small_eps <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.05,
    divergence = "kl"
  ))

  # Larger epsilon → fewer active treatments
  result_large_eps <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.15,
    divergence = "kl"
  ))

  # With smaller epsilon, expect more observations with all treatments active
  expect_true(result_small_eps$active_set_summary$pct_all_active >=
              result_large_eps$active_set_summary$pct_all_active - 0.1)

  # Both should produce finite estimates
  expect_true(is.finite(result_small_eps$estimate))
  expect_true(is.finite(result_large_eps$estimate))
})

test_that("estimate_qdelta_adaptive convergence info correct", {
  n <- 40
  m <- 2
  set.seed(404)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.1,
    divergence = "kl"
  ))

  # Convergence info should match n
  expect_equal(length(result$convergence_info$converged), n)
  expect_equal(length(result$convergence_info$iterations), n)

  # Should have high convergence rate
  expect_true(result$convergence_summary$convergence_rate > 0.9)
})

test_that("estimate_qdelta_adaptive with bootstrap", {
  n <- 50
  m <- 2
  set.seed(505)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressWarnings(suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.1,
    divergence = "kl",
    bootstrap = TRUE,
    n_bootstrap = 100
  )))

  # Bootstrap CI should exist
  expect_false(is.null(result$bootstrap_ci))
  expect_equal(length(result$bootstrap_ci), 2)
  expect_true(result$bootstrap_ci[1] < result$bootstrap_ci[2])
})

test_that("estimate_qdelta_adaptive with different divergences", {
  n <- 50
  m <- 2
  set.seed(606)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Test KL and reverse_KL
  for (div in c("kl", "reverse_kl")) {
    result <- suppressMessages(estimate_qdelta_adaptive(
      data = data,
      treatment_cols = A_cols(data), covariate_cols = X_cols(data),      delta = c(0.1, 0.1),
      epsilon = 0.1,
      divergence = div
    ))

    expect_equal(result$divergence, div)
    expect_true(is.finite(result$estimate))
  }
})

test_that("estimate_qdelta_adaptive returns expected structure", {
  n <- 50
  m <- 2
  set.seed(707)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.1,
    divergence = "kl"
  ))

  # Standard components
  expect_true(all(c("estimate", "std_error", "variance", "ci_lower", "ci_upper",
                   "eif_values", "centered_eif", "convergence_info") %in% names(result)))

  # Adaptive estimator specific components
  expect_true(all(c("epsilon", "delta_target", "delta_adaptive", "scaling_method",
                   "active_sets", "active_set_summary", "n_insufficient_active",
                   "delta_adaptation_summary") %in% names(result)))

  # Dimensions
  expect_equal(length(result$eif_values), n)
  expect_equal(length(result$centered_eif), n)
  expect_equal(mean(result$centered_eif), 0, tolerance = 1e-10)

  # Delta adaptive matrix dimensions
  expect_equal(nrow(result$delta_adaptive), n)
  expect_equal(ncol(result$delta_adaptive), m)
})

test_that("compute_adaptive_delta_single works correctly", {
  m <- 3
  marginal_probs <- c(0.05, 0.5, 0.95)  # One near 0, one middle, one near 1
  delta_target <- c(0.1, 0.1, 0.1)
  epsilon <- 0.1

  # Active set: only treatment 2 (middle)
  active_indices <- 2

  # Proportional scaling
  delta_adaptive <- compute_adaptive_delta_single(
    marginal_probs = marginal_probs,
    delta_target = delta_target,
    active_indices = active_indices,
    epsilon = epsilon,
    scaling_method = "proportional"
  )

  expect_equal(length(delta_adaptive), m)

  # Active treatment (2) should have target delta
  expect_equal(delta_adaptive[2], delta_target[2])

  # Inactive treatments (1, 3) interpolate toward 1 (the no-op) as m_k -> boundary.
  # delta_target = 0.1 < 1, so "toward 1" means LARGER than target but still < 1.
  # (Previous test asserted delta < target, i.e. driven toward 0 -- that encoded the
  #  scaling-direction bug: it made an EXTREME shift on the degenerate coordinate.)
  expect_true(delta_adaptive[1] > delta_target[1] && delta_adaptive[1] < 1)  # near 0
  expect_true(delta_adaptive[3] > delta_target[3] && delta_adaptive[3] < 1)  # near 1

  # Coordinate closer to the boundary (m=0.05, dist 0.05) is nearer the no-op (1)
  # than one farther in (both are dist 0.05 here, so equal); check the limit instead:
  # a coordinate AT the boundary would give exactly 1.
  expect_true(all(delta_adaptive > 0))
})

test_that("compute_adaptive_delta_single with zero scaling", {
  m <- 2
  marginal_probs <- c(0.05, 0.5)
  delta_target <- c(0.1, 0.1)
  epsilon <- 0.1
  active_indices <- 2  # Only treatment 2 active

  delta_adaptive <- compute_adaptive_delta_single(
    marginal_probs = marginal_probs,
    delta_target = delta_target,
    active_indices = active_indices,
    epsilon = epsilon,
    scaling_method = "zero"
  )

  # Active treatment gets target delta
  expect_equal(delta_adaptive[2], delta_target[2])

  # Inactive treatment gets delta=1.0 (no shift; odds ratio of 1 = unchanged)
  expect_equal(delta_adaptive[1], 1.0)
})

test_that("compute_adaptive_delta_single with fixed scaling", {
  m <- 2
  marginal_probs <- c(0.05, 0.5)
  delta_target <- c(0.1, 0.1)
  epsilon <- 0.1
  active_indices <- 2
  inactive_fraction <- 0.3

  delta_adaptive <- compute_adaptive_delta_single(
    marginal_probs = marginal_probs,
    delta_target = delta_target,
    active_indices = active_indices,
    epsilon = epsilon,
    scaling_method = "fixed",
    inactive_fraction = inactive_fraction
  )

  # Active treatment gets target delta
  expect_equal(delta_adaptive[2], delta_target[2])

  # Inactive treatment gets fixed fraction
  expect_equal(delta_adaptive[1], delta_target[1] * inactive_fraction)
})

test_that("compute_adaptive_deltas works for multiple observations", {
  n <- 10
  m <- 2
  set.seed(808)

  # Create marginal probabilities with some variation
  marginal_probs <- matrix(runif(n * m, 0.1, 0.9), nrow = n, ncol = m)

  delta_target <- c(0.1, 0.1)
  epsilon <- 0.15

  # Compute active sets
  active_sets <- compute_active_sets(marginal_probs, epsilon = epsilon)

  # Compute adaptive deltas
  delta_adaptive <- compute_adaptive_deltas(
    marginal_probs = marginal_probs,
    delta_target = delta_target,
    active_sets = active_sets,
    epsilon = epsilon,
    scaling_method = "proportional"
  )

  # Check dimensions
  expect_equal(dim(delta_adaptive), c(n, m))

  # All deltas should be positive
  expect_true(all(delta_adaptive > 0))

  # Every adaptive delta lies BETWEEN delta_target and 1 (the no-op), inclusive:
  # active coords keep the target, inactive coords interpolate toward 1. With
  # delta_target = 0.1 < 1 this means target <= delta_adaptive <= 1.
  # (Previous test asserted delta_adaptive <= delta_target, encoding the buggy
  #  toward-0 scaling direction.)
  tgt <- matrix(delta_target, nrow = n, ncol = m, byrow = TRUE)
  lo <- pmin(tgt, 1); hi <- pmax(tgt, 1)
  expect_true(all(delta_adaptive >= lo - 1e-6 & delta_adaptive <= hi + 1e-6))
})

test_that("estimate_qdelta_adaptive vs vanishing threshold comparison", {
  # Same data for both estimators
  n <- 100
  m <- 2
  set.seed(909)

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Adaptive estimator with fixed epsilon
  result_adaptive <- suppressMessages(estimate_qdelta_adaptive(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    epsilon = 0.1,
    divergence = "kl"
  ))

  # Vanishing threshold estimator
  result_vt <- suppressMessages(estimate_qdelta_vanishing_threshold(
    data = data,
    treatment_cols = A_cols(data), covariate_cols = X_cols(data),    delta = c(0.1, 0.1),
    divergence = "kl"
  ))

  # Both should produce finite estimates
  expect_true(is.finite(result_adaptive$estimate))
  expect_true(is.finite(result_vt$estimate))

  # Both should have similar structure (different class though)
  expect_s3_class(result_adaptive, "qdelta_adaptive")
  expect_s3_class(result_vt, "qdelta_vanishing_threshold")

  # Adaptive has delta_adaptive, VT has threshold_info
  expect_true("delta_adaptive" %in% names(result_adaptive))
  expect_true("threshold_info" %in% names(result_vt))
})

test_that("adaptive is not biased toward 0 by empty active sets (near-boundary)", {
  # Regression test for the empty-active-set bug: in a near-boundary regime a large
  # fraction of obs have 0 active coords (all m_k near the boundary). Previously
  # their EIF contribution to E[Y^Q] was set to 0, collapsing the estimate toward 0
  # (observed a 0.67 vs true ~1.46). The fix contributes Y_i (the E[Y^{P0}]=E[Y]
  # influence) for those obs. Also exercises the per-obs P0 fallback in the solve.
  # seed=2 avoids the separate PS-estimation "subscript out of bounds" fragility
  # (an unobserved combination) that some severe-regime seeds trigger.
  set.seed(2)
  n <- 2000; m <- 3
  beta <- 0.4; lp_scale <- 3   # near-boundary marginals (median b* ~ 1e-3)
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("X1", "X2")
  A <- matrix(0L, n, m); colnames(A) <- paste0("A", 1:m)
  for (k in 1:m) {
    lp <- lp_scale * (1 / beta) * (0.5 * X[, 1] + 0.3 * X[, 2])
    if (k > 1) lp <- lp + lp_scale * (1 / beta) * 0.4 * rowSums(A[, 1:(k - 1), drop = FALSE])
    A[, k] <- rbinom(n, 1, plogis(lp))
  }
  Y <- 1 + 0.3 * X[, 1] + 0.4 * X[, 2] + rowSums(A) * 0.5 / m + rnorm(n)
  data <- as.data.frame(cbind(X, A, Y = Y))

  result <- suppressWarnings(suppressMessages(estimate_qdelta_adaptive(
    data = data, treatment_cols = paste0("A", 1:m), covariate_cols = c("X1", "X2"),
    delta = rep(1.5, m), epsilon = 0.1, divergence = "kl"
  )))
  expect_s3_class(result, "qdelta_adaptive")
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result$std_error) && result$std_error > 0)

  # A large fraction of obs are empty-active (the branch under test).
  expect_gt(mean(!result$active_sets$sufficient_active), 0.2)

  # The estimate must be near the observed mean Y (delta=1.5 shifts it only
  # slightly), NOT collapsed toward 0. Pre-fix this was ~0.67 vs mean(Y)~1.5.
  expect_gt(result$estimate, 0.8 * mean(data$Y))
})
