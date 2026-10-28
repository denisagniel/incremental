library(testthat)
library(incremental)

test_that("compute_qdelta_per_obs with KL divergence works correctly", {
  # Simple 2-treatment case
  vectors <- generate_binary_vectors(2)
  p0_dist <- c(0.3, 0.2, 0.1, 0.4)
  p0_marginals <- compute_marginals(p0_dist, vectors)
  delta <- c(1.5, 2.0)

  result <- compute_qdelta_per_obs(p0_dist, delta, "kl")

  # Check that result is a valid probability distribution
  expect_equal(sum(result$q_delta), 1, tolerance = 1e-6)
  expect_true(all(result$q_delta >= 0))

  # Check that marginals are approximately correct
  computed_marginals <- compute_marginals(result$q_delta, vectors)
  target_marginals <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)
  expect_equal(computed_marginals, target_marginals, tolerance = 1e-3)

  # Check convergence
  expect_true(result$converged)
})

test_that("compute_qdelta_per_obs with reverse KL divergence works", {
  vectors <- generate_binary_vectors(2)
  p0_dist <- c(0.3, 0.2, 0.1, 0.4)
  delta <- c(1.5, 2.0)

  result <- compute_qdelta_per_obs(p0_dist, delta, "reverse_kl")

  # Check that result is a valid probability distribution
  expect_equal(sum(result$q_delta), 1, tolerance = 1e-6)
  expect_true(all(result$q_delta >= 0))

  # Check convergence
  expect_true(result$converged)
})

test_that("compute_qdelta_per_obs with chi2 divergence works", {
  vectors <- generate_binary_vectors(2)
  p0_dist <- c(0.3, 0.2, 0.1, 0.4)
  delta <- c(1.5, 2.0)

  result <- compute_qdelta_per_obs(p0_dist, delta, "chi2")

  # Check that result is a valid probability distribution
  expect_equal(sum(result$q_delta), 1, tolerance = 1e-6)
  expect_true(all(result$q_delta >= 0))

  # Check convergence
  expect_true(result$converged)
})

test_that("compute_h_function works correctly", {
  # Test KL divergence case
  a_mat <- generate_binary_vectors(2)
  lambda <- c(0.5, 0.3, 0.7)  # lambda_0, lambda_1, lambda_2

  h_result <- compute_H_vectorized_general(lambda, cbind(1,a_mat), "kl")

  index <- which(apply(cbind(1, a_mat), 1, paste0, collapse = "") == "110")

  # For KL: h(a, lambda) = exp(lambda_0 + sum(a_k * lambda_k))
  expected_h <- exp(0.5 + 1 * 0.3 + 0 * 0.7)
  expect_equal(h_result$H_mat[index], expected_h, tolerance = 1e-10)

  # Check derivatives: dh/dlambda_0 = h, dh/dlambda_k = a_k * h
  expected_derivatives <- c(expected_h, 1 * expected_h, 0 * expected_h)
  expect_equal(h_result$H_prime, h_result$H_mat, tolerance = 1e-10)
})

test_that("m=1 case matches Kennedy 2019 result", {
  # Single treatment case should match original incremental propensity score
  vectors <- generate_binary_vectors(1)
  p0_dist <- c(0.7, 0.3)  # P(A=0) = 0.7, P(A=1) = 0.3
  p0_marginals <- compute_marginals(p0_dist, vectors)
  delta <- 2.0

  result <- compute_qdelta_per_obs(p0_dist, delta, "kl")

  # For m=1, the solution should be:
  # P(A=1) = delta * p / (1 + (delta-1) * p) = 2 * 0.3 / (1 + 1 * 0.3) = 0.6/1.3 ≈ 0.462
  computed_marginals <- compute_marginals(result$q_delta, vectors)
  expected_marginal <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)
  expect_equal(computed_marginals, expected_marginal, tolerance = 1e-3)
})

test_that("m=1 closed-form solver is exact and divergence-invariant", {
  # With a single eligible treatment, Q_delta is fully determined by the
  # marginal + normalization constraints, so it must be identical (exact,
  # not just approximate) across all divergences.
  vectors <- generate_binary_vectors(1)
  p0_dist <- c(0.65, 0.35)
  delta <- 1.7
  target <- delta * 0.35 / (1 + (delta - 1) * 0.35)

  q_by_div <- lapply(c("kl", "reverse_kl", "chi2", "hellinger"), function(div) {
    res <- compute_qdelta_per_obs(p0_dist, delta, div)
    # Closed form: no iterations, exact convergence
    expect_identical(res$iterations, 0L)
    expect_true(res$converged)
    expect_equal(res$constraint_norm, 0)
    # Exact marginal (machine precision, not 1e-3)
    expect_equal(sum(res$q_delta * vectors[, 1]), target, tolerance = 1e-12)
    expect_equal(sum(res$q_delta), 1, tolerance = 1e-12)
    res$q_delta
  })

  # Divergence-invariance: all q_delta vectors coincide exactly
  for (q in q_by_div[-1]) {
    expect_equal(q, q_by_div[[1]], tolerance = 1e-12)
  }
})

test_that("m=1 closed-form matches iterative Newton solver", {
  set.seed(20260721)
  for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
    for (trial in 1:4) {
      p1 <- runif(1, 0.15, 0.85)
      p0 <- c(1 - p1, p1)
      delta <- runif(1, 0.5, 2.0)

      cf <- compute_qdelta_per_obs(p0, delta, div)               # closed-form path
      it <- newton_qdelta(p0, delta, div, fix_lambda0 = FALSE)   # iterative reference

      expect_equal(cf$q_delta, it$q_delta_weights, tolerance = 1e-6)
      expect_equal(cf$lambda, it$lambda, tolerance = 1e-6)
    }
  }
})

test_that("m=1 closed-form lambda reproduces weights via EIF h-function", {
  # Downstream EIF code re-derives weights from lambda; the closed-form lambda
  # must round-trip through compute_H_vectorized_general().
  a_ext <- cbind(1, generate_binary_vectors(1))
  for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
    p0 <- c(0.6, 0.4)
    res <- compute_qdelta_per_obs(p0, 1.5, div)
    H <- compute_H_vectorized_general(matrix(res$lambda, nrow = 1), a_ext, div)
    q_recon <- as.numeric(p0 * H$H_mat)  # internal normalization sums to 1
    expect_equal(q_recon, res$q_delta, tolerance = 1e-10)
  }
})
