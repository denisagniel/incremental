test_that("compute_vanishing_threshold computes correctly", {
  # Basic computation: ε_n = n^{-γ}
  eps <- compute_vanishing_threshold(1000, gamma = 0.3)
  expect_equal(eps, 1000^(-0.3), tolerance = 1e-10)

  # With constant
  eps2 <- compute_vanishing_threshold(1000, gamma = 0.3, constant = 2)
  expect_equal(eps2, 2 * 1000^(-0.3), tolerance = 1e-10)

  # Threshold decreases with n
  eps_small <- compute_vanishing_threshold(100, gamma = 0.3)
  eps_large <- compute_vanishing_threshold(10000, gamma = 0.3)
  expect_true(eps_small > eps_large)

  # Threshold decreases faster with larger γ
  eps_slow <- compute_vanishing_threshold(1000, gamma = 0.2)
  eps_fast <- compute_vanishing_threshold(1000, gamma = 0.4)
  expect_true(eps_slow > eps_fast)
})

test_that("compute_vanishing_threshold validates inputs", {
  # Invalid n
  expect_error(compute_vanishing_threshold(-1, gamma = 0.3))
  expect_error(compute_vanishing_threshold(0, gamma = 0.3))
  expect_error(compute_vanishing_threshold(10.5, gamma = 0.3))  # non-integer
  expect_error(compute_vanishing_threshold(c(10, 20), gamma = 0.3))  # vector

  # Invalid gamma
  expect_error(compute_vanishing_threshold(100, gamma = 0))
  expect_error(compute_vanishing_threshold(100, gamma = -0.1))
  expect_error(compute_vanishing_threshold(100, gamma = c(0.2, 0.3)))

  # Invalid constant
  expect_error(compute_vanishing_threshold(100, gamma = 0.3, constant = 0))
  expect_error(compute_vanishing_threshold(100, gamma = 0.3, constant = -1))
})

test_that("compute_vanishing_threshold handles edge cases", {
  # Very small n with large γ → large ε (should be capped)
  eps <- compute_vanishing_threshold(2, gamma = 0.5)
  expect_equal(eps, 0.49)  # Capped at 0.49

  # Warning for ε >= 0.5 (2^(-0.8) = 0.574)
  expect_warning(
    eps2 <- compute_vanishing_threshold(2, gamma = 0.8),
    "Threshold.*>= 0.5"
  )
  expect_equal(eps2, 0.49)  # Also capped
})

test_that("compute_active_sets identifies active treatments correctly", {
  # Simple case: 3 obs, 2 treatments, ε = 0.2
  # Obs 1: m = [0.5, 0.5] → both active
  # Obs 2: m = [0.1, 0.5] → only treatment 2 active (0.1 < 0.2)
  # Obs 3: m = [0.9, 0.5] → only treatment 2 active (0.9 > 0.8)
  m_probs <- matrix(c(0.5, 0.1, 0.9,
                       0.5, 0.5, 0.5), nrow = 3, ncol = 2)

  active <- compute_active_sets(m_probs, epsilon = 0.2)

  # Check active matrix
  expect_equal(active$active_matrix[1, ], c(TRUE, TRUE))
  expect_equal(active$active_matrix[2, ], c(FALSE, TRUE))
  expect_equal(active$active_matrix[3, ], c(FALSE, TRUE))

  # Check counts
  expect_equal(active$n_active, c(2, 1, 1))

  # Check sufficient
  expect_true(all(active$sufficient_active))

  # Check indices
  expect_equal(active$active_indices[[1]], c(1, 2))
  expect_equal(active$active_indices[[2]], 2)
  expect_equal(active$active_indices[[3]], 2)
})

test_that("compute_active_sets handles boundaries correctly", {
  # Test exact boundaries
  m_probs <- matrix(c(0.05, 0.95, 0.5,   # Treatment 1
                       0.1, 0.9, 0.5),    # Treatment 2
                     nrow = 3, ncol = 2)
  eps <- 0.1

  active <- compute_active_sets(m_probs, epsilon = eps)

  # Obs 1: m1=0.05 < 0.1 (inactive), m2=0.1 = 0.1 (active, boundary)
  expect_equal(active$active_matrix[1, ], c(FALSE, TRUE))

  # Obs 2: m1=0.95 > 0.9 (inactive), m2=0.9 = 0.9 (active, boundary)
  expect_equal(active$active_matrix[2, ], c(FALSE, TRUE))

  # Obs 3: both at 0.5 (clearly active)
  expect_equal(active$active_matrix[3, ], c(TRUE, TRUE))
})

test_that("compute_active_sets handles min_active requirement", {
  # 3 obs, 3 treatments, require at least 2 active
  m_probs <- matrix(c(0.5, 0.5, 0.5,   # Treatment 1
                       0.5, 0.5, 0.05,  # Treatment 2
                       0.5, 0.05, 0.05), # Treatment 3
                     nrow = 3, ncol = 3)

  active <- compute_active_sets(m_probs, epsilon = 0.1, min_active = 2)

  # Obs 1: 3 active → sufficient
  expect_true(active$sufficient_active[1])
  expect_equal(active$n_active[1], 3)

  # Obs 2: 2 active → sufficient
  expect_true(active$sufficient_active[2])
  expect_equal(active$n_active[2], 2)

  # Obs 3: only 1 active → insufficient, mark all as inactive
  expect_false(active$sufficient_active[3])
  expect_equal(active$n_active[3], 0)
  expect_equal(active$active_matrix[3, ], c(FALSE, FALSE, FALSE))
  expect_equal(active$active_indices[[3]], integer(0))
})

test_that("compute_active_sets handles vector epsilon", {
  # Different threshold per observation
  m_probs <- matrix(c(0.15, 0.15,   # Treatment 1
                       0.5, 0.5),     # Treatment 2
                     nrow = 2, ncol = 2)

  # Obs 1: ε = 0.1, treatment 1 is active (0.15 >= 0.1)
  # Obs 2: ε = 0.2, treatment 1 is inactive (0.15 < 0.2)
  epsilon_vec <- c(0.1, 0.2)
  active <- compute_active_sets(m_probs, epsilon = epsilon_vec)

  expect_equal(active$active_matrix[1, ], c(TRUE, TRUE))
  expect_equal(active$active_matrix[2, ], c(FALSE, TRUE))
})

test_that("compute_active_sets validates inputs", {
  m_probs <- matrix(runif(20, 0.1, 0.9), nrow = 10, ncol = 2)

  # Invalid marginal_probs
  expect_error(compute_active_sets("not a matrix", epsilon = 0.1))
  expect_error(compute_active_sets(as.data.frame(m_probs), epsilon = 0.1))

  # Probabilities out of range
  m_bad <- m_probs
  m_bad[1, 1] <- 1.5
  expect_error(compute_active_sets(m_bad, epsilon = 0.1))

  m_bad2 <- m_probs
  m_bad2[2, 1] <- -0.1
  expect_error(compute_active_sets(m_bad2, epsilon = 0.1))

  # Invalid epsilon
  expect_error(compute_active_sets(m_probs, epsilon = 0))
  expect_error(compute_active_sets(m_probs, epsilon = 0.5))
  expect_error(compute_active_sets(m_probs, epsilon = 0.6))
  expect_error(compute_active_sets(m_probs, epsilon = c(0.1, 0.2, 0.3)))  # wrong length

  # Invalid min_active
  expect_error(compute_active_sets(m_probs, epsilon = 0.1, min_active = 0))
  expect_error(compute_active_sets(m_probs, epsilon = 0.1, min_active = 3))  # > m
  expect_error(compute_active_sets(m_probs, epsilon = 0.1, min_active = 1.5))  # non-integer
})


test_that("summarize_active_sets computes statistics correctly", {
  # Create test data: 10 obs, 3 treatments
  set.seed(123)
  m_probs <- matrix(runif(30, 0.1, 0.9), nrow = 10, ncol = 3)
  active <- compute_active_sets(m_probs, epsilon = 0.15, min_active = 1)

  summary <- summarize_active_sets(active)

  # Basic counts
  expect_equal(summary$n_obs, 10)
  expect_equal(summary$m, 3)
  expect_equal(summary$epsilon, 0.15)

  # n_active_summary should be summary stats
  expect_equal(length(summary$n_active_summary), 6)  # min, Q1, median, mean, Q3, max

  # Insufficient observations
  expect_true(summary$n_insufficient >= 0)
  expect_true(summary$n_insufficient <= 10)
  expect_equal(summary$prop_insufficient, summary$n_insufficient / 10)

  # By treatment
  expect_equal(length(summary$active_by_treatment), 3)
  expect_equal(length(summary$prop_active_by_treatment), 3)
  expect_true(all(summary$active_by_treatment >= 0))
  expect_true(all(summary$active_by_treatment <= 10))
})

test_that("summarize_active_sets includes propensity diagnostics", {
  set.seed(456)
  m_probs <- matrix(runif(30, 0.1, 0.9), nrow = 10, ncol = 3)
  active <- compute_active_sets(m_probs, epsilon = 0.15)

  # Without marginal_probs
  summary1 <- summarize_active_sets(active)
  expect_null(summary1$min_propensity)
  expect_null(summary1$max_propensity)

  # With marginal_probs
  summary2 <- summarize_active_sets(active, marginal_probs = m_probs)
  expect_equal(length(summary2$min_propensity), 3)
  expect_equal(length(summary2$max_propensity), 3)
  expect_true(all(summary2$min_propensity >= 0))
  expect_true(all(summary2$max_propensity <= 1))
  expect_true(all(summary2$min_propensity <= summary2$max_propensity))
})

test_that("active set functions integrate correctly", {
  # Full workflow test
  n <- 100
  m <- 4
  set.seed(789)

  # Simulate marginal propensities
  m_probs <- matrix(runif(n * m, 0.05, 0.95), nrow = n, ncol = m)

  # Compute vanishing threshold
  gamma <- 0.3
  eps_n <- compute_vanishing_threshold(n, gamma = gamma)

  # Compute active sets
  active <- compute_active_sets(m_probs, epsilon = eps_n, min_active = 2)

  # Summarize
  summary <- summarize_active_sets(active, marginal_probs = m_probs)

  # Sanity checks
  expect_equal(summary$n_obs, n)
  expect_equal(summary$m, m)
  expect_true(summary$prop_insufficient >= 0)
  expect_true(summary$prop_insufficient <= 1)
  expect_equal(length(active$active_indices), n)
})
