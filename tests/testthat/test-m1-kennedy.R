library(testthat)
library(incremental)

test_that("m=1 case matches Kennedy 2019 formula exactly", {
  # For m=1 with KL divergence, the solution should match the original
  # incremental propensity score intervention from Kennedy 2019
  
  vectors <- generate_binary_vectors(1)
  p0_dist <- c(0.7, 0.3)  # P(A=0)=0.7, P(A=1)=0.3
  p0_marginals <- compute_marginals(p0_dist, vectors)  # = 0.3
  delta <- 2.0
  
  result <- compute_qdelta_per_obs(p0_dist, delta, "kl")
  
  # For m=1 with KL, the solution has a known form:
  # q(1) = p_0(1) * exp(log(delta)) / Z = p_0(1) * delta / Z
  # q(0) = p_0(0) / Z
  # where Z = p_0(0) + p_0(1) * delta = 0.7 + 0.3 * 2 = 1.3
  
  # Expected distribution
  Z_expected <- p0_dist[1] + p0_dist[2] * delta  # 0.7 + 0.3*2 = 1.3
  expected_q0 <- p0_dist[1] / Z_expected  # 0.7 / 1.3 ≈ 0.5385
  expected_q1 <- p0_dist[2] * delta / Z_expected  # 0.3 * 2 / 1.3 ≈ 0.4615
  
  expect_equal(result$q_delta[1], expected_q0, tolerance = 1e-4)
  expect_equal(result$q_delta[2], expected_q1, tolerance = 1e-4)
  
  # Expected marginal: delta * p / (1 + (delta-1)*p) = 2*0.3/(1+1*0.3) = 0.6/1.3 ≈ 0.4615
  expected_marginal <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)
  computed_marginals <- compute_marginals(result$q_delta, vectors)
  expect_equal(computed_marginals, expected_marginal, tolerance = 1e-4)
  
  # For m=1, lambda_1 should be approximately log(delta) = log(2) ≈ 0.693
  # (lambda_0 is the normalizing constant)
  expect_equal(result$lambda[2], log(delta), tolerance = 1e-2)
  
  cat("✓ m=1 case matches Kennedy 2019 formula\n")
  cat(sprintf("  P_0(A=1) = %.4f\n", p0_marginals))
  cat(sprintf("  delta = %.4f\n", delta))
  cat(sprintf("  Expected Q(A=1) = %.4f\n", expected_marginal))
  cat(sprintf("  Computed Q(A=1) = %.4f\n", computed_marginals))
  cat(sprintf("  Expected lambda_1 = log(delta) = %.4f\n", log(delta)))
  cat(sprintf("  Computed lambda_1 = %.4f\n", result$lambda[2]))
  cat(sprintf("  Distribution: q(0)=%.4f, q(1)=%.4f\n", result$q_delta[1], result$q_delta[2]))
})

test_that("m=1 case works for different delta values", {
  vectors <- generate_binary_vectors(1)
  p0_dist <- c(0.6, 0.4)
  p0_marginals <- compute_marginals(p0_dist, vectors)
  
  delta_values <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  
  for (delta in delta_values) {
    result <- compute_qdelta_per_obs(p0_dist, delta, "kl")
    
    # Check distribution form
    Z_expected <- p0_dist[1] + p0_dist[2] * delta
    expected_q0 <- p0_dist[1] / Z_expected
    expected_q1 <- p0_dist[2] * delta / Z_expected
    
    expect_equal(result$q_delta[1], expected_q0, tolerance = 1e-4)
    expect_equal(result$q_delta[2], expected_q1, tolerance = 1e-4)
    
    # Check marginals
    expected_marginal <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)
    computed_marginals <- compute_marginals(result$q_delta, vectors)
    expect_equal(computed_marginals, expected_marginal, tolerance = 1e-4)
    
    cat(sprintf("✓ delta=%.1f: Q(A=1)=%.4f (expected=%.4f)\n", 
                delta, computed_marginals, expected_marginal))
  }
})


