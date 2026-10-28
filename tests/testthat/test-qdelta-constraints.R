library(testthat)
library(incremental)

test_that("Q_delta satisfies marginal constraints", {
  # Test that Q_delta distribution satisfies:
  # 1. Sums to 1 (valid probability distribution)
  # 2. All probabilities non-negative
  # 3. Marginal probabilities match target marginals (within tolerance)
  
  vectors <- generate_binary_vectors(2)
  p0_dist <- c(0.3, 0.2, 0.1, 0.4)
  p0_marginals <- compute_marginals(p0_dist, vectors)
  delta <- c(1.5, 2.0)
  
  result <- compute_qdelta_per_obs(p0_dist, delta, "kl")
  
  # Check 1: Distribution sums to 1
  expect_equal(sum(result$q_delta), 1, tolerance = 1e-6)
  
  # Check 2: All probabilities are non-negative
  expect_true(all(result$q_delta >= 0))
  expect_true(all(result$q_delta <= 1))
  
  # Check 3: Marginal probabilities match target marginals
  # Target marginals from odds scaling formula
  target_marginals <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)
  
  # Compute actual marginals from Q_delta distribution
  computed_marginals <- compute_marginals(result$q_delta, vectors)
  
  expect_equal(computed_marginals, target_marginals, tolerance = 1e-4)
  
  # Check 4: Convergence was successful
  expect_true(result$converged)
  
  cat("✓ Q_delta satisfies all constraints\n")
  cat(sprintf("  Distribution sums to: %.6f (should be 1.0)\n", sum(result$q_delta)))
  cat(sprintf("  All probabilities >= 0: %s\n", all(result$q_delta >= 0)))
  cat(sprintf("  Target marginals: (%.4f, %.4f)\n", target_marginals[1], target_marginals[2]))
  cat(sprintf("  Computed marginals: (%.4f, %.4f)\n", computed_marginals[1], computed_marginals[2]))
  cat(sprintf("  Max difference: %.6f (should be < 1e-4)\n", max(abs(computed_marginals - target_marginals))))
  cat(sprintf("  Convergence: %s\n", result$converged))
})

test_that("Q_delta works for different divergences", {
  vectors <- generate_binary_vectors(2)
  p0_dist <- c(0.3, 0.2, 0.1, 0.4)
  p0_marginals <- compute_marginals(p0_dist, vectors)
  delta <- c(1.5, 2.0)
  
  divergences <- c("kl", "reverse_kl", "chi2", "hellinger")
  
  for (div in divergences) {
    result <- compute_qdelta_per_obs(p0_dist, delta, div)
    
    # All should satisfy basic constraints
    expect_equal(sum(result$q_delta), 1, tolerance = 1e-5)
    expect_true(all(result$q_delta >= 0))
    expect_true(result$converged)
    
    # Check marginals match target (may be less accurate for some divergences)
    target_marginals <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)
    computed_marginals <- compute_marginals(result$q_delta, vectors)
    
    # Allow slightly larger tolerance for non-KL divergences
    tol <- ifelse(div == "kl", 1e-4, 1e-3)
    expect_equal(computed_marginals, target_marginals, tolerance = tol)
    
    cat(sprintf("✓ %s divergence: constraints satisfied (tolerance: %.0e)\n", div, tol))
  }
})


