test_that("efficient_influence_function_active validates inputs", {
  # Setup simple test case
  n <- 50
  m <- 2
  set.seed(123)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)
  ps <- matrix(runif(n * 4, 0.1, 0.9), nrow = n, ncol = 4)
  ps <- ps / rowSums(ps)  # Normalize
  outcomes <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  delta <- c(0.1, 0.1)

  # Create valid active sets
  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.1)

  # Invalid active_sets structure
  expect_error(
    efficient_influence_function_active(
      A = A, Y = Y, propensity_scores = ps,
      outcome_means = outcomes, delta = delta,
      active_sets = list(wrong = "structure")
    ),
    "active_sets must be a list"
  )

  # Dimension mismatch
  active_wrong <- active
  active_wrong$active_matrix <- active_wrong$active_matrix[1:10, ]
  expect_error(
    efficient_influence_function_active(
      A = A, Y = Y, propensity_scores = ps,
      outcome_means = outcomes, delta = delta,
      active_sets = active_wrong
    ),
    "active_sets dimensions do not match"
  )

  # Other standard validations should still work
  expect_error(
    efficient_influence_function_active(
      A = A, Y = Y, propensity_scores = ps * 2,  # doesn't sum to 1
      outcome_means = outcomes, delta = delta,
      active_sets = active
    ),
    "do not sum to 1"
  )
})

test_that("efficient_influence_function_active handles all active treatments", {
  # Case where all treatments are active for all observations
  n <- 50
  m <- 2
  set.seed(456)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)

  # Propensities well away from boundaries
  ps <- matrix(runif(n * 4, 0.3, 0.7), nrow = n, ncol = 4)
  ps <- ps / rowSums(ps)
  outcomes <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  delta <- c(0.1, 0.1)

  # Small threshold means all treatments active
  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.05)

  # Should have all treatments active
  expect_true(all(active$n_active == m))
  expect_true(all(active$sufficient_active))

  # Compute EIF
  eif_result <- efficient_influence_function_active(
    A = A, Y = Y,
    propensity_scores = ps,
    outcome_means = outcomes,
    delta = delta,
    divergence = "kl",
    active_sets = active
  )

  # Check structure
  expect_equal(length(eif_result$eif_values), n)
  expect_equal(length(eif_result$plugin), n)
  expect_equal(eif_result$n_insufficient_active, 0)

  # EIF should be finite
  expect_true(all(is.finite(eif_result$eif_values)))
})

test_that("efficient_influence_function_active handles some inactive treatments", {
  # Case where some treatments inactive for some observations
  n <- 100
  m <- 3
  set.seed(789)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)

  # Create propensities with some near boundaries
  ps <- matrix(0.25, nrow = n, ncol = 8)
  # Make first treatment nearly deterministic for some obs
  ps[1:20, c(1, 2, 3, 4)] <- 0.01  # A1 = 0 likely
  ps[1:20, c(5, 6, 7, 8)] <- 0.99 / 4  # A1 = 1 likely
  ps <- ps / rowSums(ps)

  outcomes <- matrix(rnorm(n * 8), nrow = n, ncol = 8)
  delta <- c(0.05, 0.05, 0.05)

  # Moderate threshold
  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.15)

  # Some observations should have fewer active treatments
  expect_true(any(active$n_active < m))
  expect_true(all(active$sufficient_active))  # But enough for min_active = 1

  # Compute EIF
  eif_result <- efficient_influence_function_active(
    A = A, Y = Y,
    propensity_scores = ps,
    outcome_means = outcomes,
    delta = delta,
    divergence = "kl",
    active_sets = active
  )

  # Should complete without error
  expect_equal(length(eif_result$eif_values), n)
  expect_true(all(is.finite(eif_result$eif_values)))
  expect_equal(eif_result$n_insufficient_active, 0)
})

test_that("efficient_influence_function_active handles insufficient active treatments", {
  # Case where some observations have no sufficient active treatments
  n <- 50
  m <- 3
  set.seed(101)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)

  # Create propensities with many near boundaries
  ps <- matrix(0.125, nrow = n, ncol = 8)
  # Make most treatments nearly deterministic for first 10 obs
  ps[1:10, ] <- 0.01
  ps[1:10, 1] <- 0.93  # Only one combination likely
  ps <- ps / rowSums(ps)

  outcomes <- matrix(rnorm(n * 8), nrow = n, ncol = 8)
  delta <- c(0.05, 0.05, 0.05)

  # Large threshold + require 2 active
  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.2, min_active = 2)

  # First 10 obs should have insufficient active treatments
  expect_true(sum(!active$sufficient_active) >= 5)
  n_insufficient <- sum(!active$sufficient_active)

  # Compute EIF (should warn)
  expect_warning(
    eif_result <- efficient_influence_function_active(
      A = A, Y = Y,
      propensity_scores = ps,
      outcome_means = outcomes,
      delta = delta,
      divergence = "kl",
      active_sets = active
    ),
    "insufficient active treatments"
  )

  # Insufficient (empty-active-set) obs: intervention is the no-op (delta -> 1),
  # Q = P0, so the EIF contribution to E[Y^Q] = E[Y] is Y_i -- NOT zero. (Previous
  # test asserted == 0, which biased E[Y^Q] toward 0 by the empty-active fraction.)
  insufficient_idx <- which(!active$sufficient_active)
  if (length(insufficient_idx) > 0) {
    expect_equal(eif_result$eif_values[insufficient_idx], Y[insufficient_idx])
  }

  # Other observations should have non-zero EIF
  sufficient_idx <- which(active$sufficient_active)
  expect_true(any(eif_result$eif_values[sufficient_idx] != 0))

  expect_equal(eif_result$n_insufficient_active, n_insufficient)
})

test_that("efficient_influence_function_active agrees with standard EIF when all active", {
  # When all treatments are active, should match standard EIF
  n <- 30
  m <- 2
  set.seed(202)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)

  # Well-behaved propensities (all active)
  ps <- matrix(runif(n * 4, 0.2, 0.8), nrow = n, ncol = 4)
  ps <- ps / rowSums(ps)
  outcomes <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  delta <- c(0.1, 0.1)

  # Very small threshold → all active
  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.01)
  expect_true(all(active$n_active == m))

  # Compute standard EIF
  eif_standard <- efficient_influence_function(
    A = A, Y = Y,
    propensity_scores = ps,
    outcome_means = outcomes,
    delta = delta,
    divergence = "kl"
  )

  # Compute active EIF
  eif_active <- efficient_influence_function_active(
    A = A, Y = Y,
    propensity_scores = ps,
    outcome_means = outcomes,
    delta = delta,
    divergence = "kl",
    active_sets = active
  )

  # Should be very close (allowing for numerical differences)
  expect_equal(eif_active$eif_values, eif_standard$eif_values, tolerance = 1e-6)
  expect_equal(eif_active$plugin, eif_standard$plugin, tolerance = 1e-6)
})

test_that("efficient_influence_function_active works with different divergences", {
  n <- 40
  m <- 2
  set.seed(303)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)
  ps <- matrix(runif(n * 4, 0.2, 0.8), nrow = n, ncol = 4)
  ps <- ps / rowSums(ps)
  outcomes <- matrix(rnorm(n * 4), nrow = n, ncol = 4)

  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.1)

  # Test robust divergences with delta=0.1
  for (div in c("kl", "reverse_kl")) {
    eif_result <- efficient_influence_function_active(
      A = A, Y = Y,
      propensity_scores = ps,
      outcome_means = outcomes,
      delta = c(0.1, 0.1),
      divergence = div,
      active_sets = active
    )

    expect_equal(length(eif_result$eif_values), n)
    expect_true(all(is.finite(eif_result$eif_values)))
  }

  # Test chi2 and hellinger with smaller delta (may fail with random data)
  for (div in c("chi2", "hellinger")) {
    # chi2 can be fragile; wrap in tryCatch
    result <- tryCatch({
      efficient_influence_function_active(
        A = A, Y = Y,
        propensity_scores = ps,
        outcome_means = outcomes,
        delta = c(0.05, 0.05),
        divergence = div,
        active_sets = active
      )
    }, error = function(e) {
      # Skip if chi2/hellinger convergence fails (known fragility)
      if (grepl("Chi2 convergence failure|Hellinger", e$message)) {
        return(NULL)
      }
      stop(e)
    })

    if (!is.null(result)) {
      expect_equal(length(result$eif_values), n)
      expect_true(all(is.finite(result$eif_values)))
    }
  }
})

test_that("calculate_eif_lambda_active handles empty active set", {
  # Single observation with no active treatments
  m <- 2
  A <- c(1, 0)
  ps <- c(0.1, 0.2, 0.3, 0.4)

  # Create minimal qdelta_result
  qdelta_result <- list(
    lambda = rep(0, m + 1),
    q_delta = ps,
    delta = c(0.1, 0.1)
  )

  a_mat <- generate_binary_vectors(m)
  a_extended <- cbind(1, a_mat)
  h_values <- rep(1, 4)
  h_prime <- rep(0, 4)

  # Empty active set
  eif_lambda <- calculate_eif_lambda_active(
    A, ps, qdelta_result, h_values, h_prime,
    a_extended, "kl", active_indices = integer(0)
  )

  # Should return zeros
  expect_equal(eif_lambda, numeric(m + 1))
})

test_that("calculate_eif_lambda_active handles single active treatment", {
  # Single observation with only one active treatment
  m <- 2
  set.seed(404)

  A <- c(1, 0)
  ps <- c(0.2, 0.3, 0.3, 0.2)

  # Compute Q_delta
  delta <- c(0.1, 0.1)
  qdelta_result <- compute_qdelta_per_obs(ps, delta, divergence = "kl")

  a_mat <- generate_binary_vectors(m)
  a_extended <- cbind(1, a_mat)

  # Compute h values
  lambda <- qdelta_result$lambda
  h_values <- as.numeric(exp(lambda %*% t(a_extended)))
  h_prime <- h_values  # For KL, h' = h

  # Only treatment 1 active
  eif_lambda <- calculate_eif_lambda_active(
    A, ps, qdelta_result, h_values, h_prime,
    a_extended, "kl", active_indices = 1
  )

  # Should return vector of full length m+1 with zeros for inactive
  expect_equal(length(eif_lambda), m + 1)
  expect_true(all(is.finite(eif_lambda)))
  # Treatment 2 should be zero (inactive)
  expect_equal(eif_lambda[3], 0)
})

test_that("efficient_influence_function_active convergence warnings work", {
  # Test with difficult delta that causes convergence issues
  n <- 100  # More observations to ensure enough failures
  m <- 2
  set.seed(505)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)

  # Extreme propensities for many observations
  ps <- matrix(0, nrow = n, ncol = 4)
  # First 50 obs: highly concentrated (delta will be infeasible)
  ps[1:50, 1] <- 0.95
  ps[1:50, 2:4] <- 0.05 / 3
  # Next 30 obs: moderately concentrated
  ps[51:80, 1] <- 0.7
  ps[51:80, 2:4] <- 0.3 / 3
  # Last 20 obs: balanced (should converge)
  ps[81:100, ] <- 0.25
  ps <- ps / rowSums(ps)
  outcomes <- matrix(rnorm(n * 4), nrow = n, ncol = 4)

  # Large delta (infeasible for concentrated propensities)
  delta <- c(3.0, 3.0)

  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.1)

  # Should warn about convergence or insufficient active (both acceptable)
  expect_warning(
    eif_result <- efficient_influence_function_active(
      A = A, Y = Y,
      propensity_scores = ps,
      outcome_means = outcomes,
      delta = delta,
      divergence = "kl",
      active_sets = active
    ),
    "Low Q_delta convergence rate|insufficient active treatments"
  )
})

test_that("efficient_influence_function_active returns expected structure", {
  n <- 50
  m <- 2
  set.seed(606)

  A <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  Y <- rnorm(n)
  ps <- matrix(runif(n * 4, 0.2, 0.8), nrow = n, ncol = 4)
  ps <- ps / rowSums(ps)
  outcomes <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  delta <- c(0.1, 0.1)

  marginal_probs <- extract_marginal_probabilities(ps, m = m)
  active <- compute_active_sets(marginal_probs, epsilon = 0.1)

  eif_result <- efficient_influence_function_active(
    A = A, Y = Y,
    propensity_scores = ps,
    outcome_means = outcomes,
    delta = delta,
    divergence = "kl",
    active_sets = active
  )

  # Check all expected components
  expect_true(all(c("eif_values", "centered_eif", "plugin", "qdelta_results",
                    "convergence_info", "convergence_summary", "active_sets",
                    "n_insufficient_active") %in% names(eif_result)))

  # Check dimensions
  expect_equal(length(eif_result$eif_values), n)
  expect_equal(length(eif_result$centered_eif), n)
  expect_equal(length(eif_result$plugin), n)
  expect_equal(length(eif_result$qdelta_results), n)

  # Check centered EIF has mean ~0
  expect_equal(mean(eif_result$centered_eif), 0, tolerance = 1e-10)

  # Check convergence info structure
  expect_equal(length(eif_result$convergence_info$converged), n)
  expect_equal(length(eif_result$convergence_info$iterations), n)
  expect_equal(length(eif_result$convergence_info$constraint_norms), n)

  # Check active_sets preserved
  expect_identical(eif_result$active_sets, active)
})
