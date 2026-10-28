library(testthat)
library(incremental)

test_that("Target marginals use correct odds scaling formula", {
  # Test that compute_target_marginals uses the odds scaling formula
  # from theory (main.tex line 81): 
  # Q_delta(A_k = 1) = delta_k * P_0(A_k = 1) / (1 + (delta_k - 1) * P_0(A_k = 1))
  
  # Simple test case
  p0_marginals <- c(0.3, 0.4)
  delta <- c(1.5, 2.0)
  
  # Expected target marginals using correct formula
  expected_target <- c(
    1.5 * 0.3 / (1 + (1.5 - 1) * 0.3),  # = 0.45 / (1 + 0.15) = 0.45 / 1.15 ≈ 0.3913
    2.0 * 0.4 / (1 + (2.0 - 1) * 0.4)   # = 0.80 / (1 + 0.40) = 0.80 / 1.40 ≈ 0.5714
  )
  
  # Create a propensity score vector that matches these marginals
  vectors <- generate_binary_vectors(2)
  # For marginals (0.3, 0.4), we need a distribution where:
  # P(A1=1) = 0.3, P(A2=1) = 0.4
  # One such distribution: (0,0)=0.4, (0,1)=0.3, (1,0)=0.1, (1,1)=0.2
  # Check: P(A1=1) = 0.1 + 0.2 = 0.3 ✓, P(A2=1) = 0.3 + 0.2 = 0.5 ✗
  # Better: (0,0)=0.4, (0,1)=0.3, (1,0)=0.0, (1,1)=0.3
  # Check: P(A1=1) = 0.0 + 0.3 = 0.3 ✓, P(A2=1) = 0.3 + 0.3 = 0.6 ✗
  # Let's use: (0,0)=0.5, (0,1)=0.2, (1,0)=0.2, (1,1)=0.1
  # Check: P(A1=1) = 0.2 + 0.1 = 0.3 ✓, P(A2=1) = 0.2 + 0.1 = 0.3 ✗
  # Actually, let's compute it properly:
  # We want P(A1=1) = p(1,0) + p(1,1) = 0.3
  # We want P(A2=1) = p(0,1) + p(1,1) = 0.4
  # One solution: p(0,0)=0.3, p(0,1)=0.3, p(1,0)=0.0, p(1,1)=0.4
  # Check: P(A1=1) = 0.0 + 0.4 = 0.4 ✗
  # Let's solve: p(0,0) + p(0,1) + p(1,0) + p(1,1) = 1
  # p(1,0) + p(1,1) = 0.3
  # p(0,1) + p(1,1) = 0.4
  # One solution: p(0,0)=0.3, p(0,1)=0.3, p(1,0)=0.0, p(1,1)=0.4
  # Wait, that gives P(A1=1)=0.4, not 0.3
  # Let's try: p(0,0)=0.4, p(0,1)=0.3, p(1,0)=0.0, p(1,1)=0.3
  # P(A1=1) = 0.0 + 0.3 = 0.3 ✓, P(A2=1) = 0.3 + 0.3 = 0.6 ✗
  # Let's use a simpler approach - just create any distribution and compute marginals from it
  p0_dist <- c(0.4, 0.2, 0.1, 0.3)
  computed_marginals <- compute_marginals(p0_dist, vectors)
  
  # Now test compute_target_marginals with this distribution
  result <- compute_target_marginals(p0_dist, vectors, delta)
  
  # The result should use the odds scaling formula on the computed marginals
  expected_from_computed <- delta * computed_marginals / (1 + (delta - 1) * computed_marginals)
  
  expect_equal(result, expected_from_computed, tolerance = 1e-6)
  
  # Also verify it's NOT using the wrong formula (addition)
  wrong_formula_result <- computed_marginals + delta
  expect_false(all(abs(result - wrong_formula_result) < 1e-6))
  
  cat("✓ Target marginals use correct odds scaling formula\n")
  cat(sprintf("  Computed marginals: (%.4f, %.4f)\n", computed_marginals[1], computed_marginals[2]))
  cat(sprintf("  Target marginals: (%.4f, %.4f)\n", result[1], result[2]))
  cat(sprintf("  Expected (odds scaling): (%.4f, %.4f)\n", expected_from_computed[1], expected_from_computed[2]))
  cat(sprintf("  Wrong formula (addition) would give: (%.4f, %.4f)\n", wrong_formula_result[1], wrong_formula_result[2]))
})



