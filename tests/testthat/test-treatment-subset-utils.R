test_that("create_reduced_treatment_matrix works correctly", {
  # Generate full treatment matrix for m=4
  vectors <- generate_binary_vectors(4)  # 16 x 4

  # Extract treatments 2 and 4
  active <- c(2, 4)
  reduced <- create_reduced_treatment_matrix(vectors, active)

  # Check dimensions
  expect_equal(nrow(reduced), 16)
  expect_equal(ncol(reduced), 2)

  # Check content matches
  expect_equal(reduced[, 1], vectors[, 2])
  expect_equal(reduced[, 2], vectors[, 4])

  # Test error handling
  expect_error(create_reduced_treatment_matrix(vectors, c(0, 2)))
  expect_error(create_reduced_treatment_matrix(vectors, c(2, 5)))
  expect_error(create_reduced_treatment_matrix(vectors, integer(0)))
})

test_that("create_a_bar_matrix adds intercept correctly", {
  vectors <- generate_binary_vectors(3)  # 8 x 3
  active <- c(1, 3)

  a_bar <- create_a_bar_matrix(vectors, active)

  # Check dimensions
  expect_equal(nrow(a_bar), 8)
  expect_equal(ncol(a_bar), 3)  # intercept + 2 active

  # Check intercept column
  expect_equal(a_bar[, 1], rep(1, 8))

  # Check treatment columns
  expect_equal(a_bar[, 2], vectors[, 1])
  expect_equal(a_bar[, 3], vectors[, 3])

  # Check column names
  expect_equal(colnames(a_bar), c("intercept", "A1", "A3"))
})

test_that("extract_marginal_probabilities computes correctly", {
  # Simple case: m=2, uniform distribution
  # P(A1=0, A2=0) = 0.25, P(A1=1, A2=0) = 0.25
  # P(A1=0, A2=1) = 0.25, P(A1=1, A2=1) = 0.25
  ps <- matrix(0.25, nrow = 10, ncol = 4)

  marginals <- extract_marginal_probabilities(ps, m = 2)

  # Each marginal should be 0.5
  expect_equal(dim(marginals), c(10, 2))
  expect_true(all(abs(marginals - 0.5) < 1e-10))

  # Test with skewed distribution
  # P(0,0)=0.4, P(1,0)=0.3, P(0,1)=0.2, P(1,1)=0.1
  ps2 <- matrix(c(0.4, 0.3, 0.2, 0.1), nrow = 1, ncol = 4)
  marginals2 <- extract_marginal_probabilities(ps2, m = 2)

  # m1 = P(A1=1) = P(1,0) + P(1,1) = 0.3 + 0.1 = 0.4
  # m2 = P(A2=1) = P(0,1) + P(1,1) = 0.2 + 0.1 = 0.3
  expect_equal(unname(marginals2[1, 1]), 0.4, tolerance = 1e-10)
  expect_equal(unname(marginals2[1, 2]), 0.3, tolerance = 1e-10)
})

test_that("compute_omega_reduced produces positive definite matrix", {
  # Simple test case: m=2, active = {1}
  # So we have intercept + A1 only (dimension 2)

  # Propensity scores for 4 combinations (must sum to 1)
  ps <- c(0.3, 0.2, 0.4, 0.1)

  # Lambda: intercept + A1 coefficient
  lambda <- c(0.5, 0.3)

  # Create a_bar for active set {1}
  vectors <- generate_binary_vectors(2)
  a_bar <- create_a_bar_matrix(vectors, c(1))

  # Compute omega (should be 2x2)
  omega <- compute_omega_reduced(ps, lambda, a_bar, "kl")

  # Check dimensions
  expect_equal(dim(omega), c(2, 2))

  # Check symmetric
  expect_equal(omega[1, 2], omega[2, 1], tolerance = 1e-10)

  # Check positive definite (all eigenvalues > 0)
  eigenvals <- eigen(omega, only.values = TRUE)$values
  expect_true(all(eigenvals > 0))
})

test_that("compute_omega_reduced works for different divergences", {
  ps <- c(0.25, 0.25, 0.25, 0.25)
  lambda <- c(0, 0)
  vectors <- generate_binary_vectors(2)
  a_bar <- create_a_bar_matrix(vectors, c(1))

  # All divergences should produce valid matrices
  for (div in c("kl", "chi2", "hellinger", "reverse_kl")) {
    # [2026-08-11] reverse_kl's domain is x < 0, and lambda = c(0, 0) puts every
    # linear term at exactly 0 -- i.e. ON the singularity, outside the domain. The
    # pre-bugfix code used the wrong conjugate (h'(z) = 1/(1-z)^2) which returned an
    # innocuous 1 there and silently hid that; the canonical divergence_h_pair()
    # correctly saturates and warns. This test's intent is "every divergence returns
    # a valid 2x2 finite matrix", not "reverse_kl accepts out-of-domain input", so
    # reverse_kl gets a lambda inside its own domain. Same assertions for all four.
    lambda_div <- if (div == "reverse_kl") c(-1, -0.2) else lambda
    omega <- compute_omega_reduced(ps, lambda_div, a_bar, div)
    expect_equal(dim(omega), c(2, 2))
    expect_true(all(is.finite(omega)))
  }
})

test_that("extend_eif_to_full places values correctly", {
  # Simple case: subset of 3 observations from 10
  eif_subset <- c(1.5, 2.0, -0.5)
  indices <- c(2, 5, 8)

  eif_full <- extend_eif_to_full(eif_subset, indices, 10)

  # Check length
  expect_equal(length(eif_full), 10)

  # Check values at correct positions
  expect_equal(eif_full[2], 1.5)
  expect_equal(eif_full[5], 2.0)
  expect_equal(eif_full[8], -0.5)

  # Check zeros elsewhere
  expect_equal(eif_full[c(1, 3, 4, 6, 7, 9, 10)], rep(0, 7))

  # Test error handling
  expect_error(extend_eif_to_full(eif_subset, c(1, 2), 10))  # length mismatch
  expect_error(extend_eif_to_full(eif_subset, c(2, 5, 11), 10))  # out of range
  expect_error(extend_eif_to_full(eif_subset, c(2, 5, 5), 10))  # duplicates
})

test_that("map_treatment_indices creates correct mapping", {
  active <- c(2, 5, 7)
  mapping <- map_treatment_indices(active, m = 10)

  expect_equal(nrow(mapping), 3)
  expect_equal(mapping$original_index, c(2, 5, 7))
  expect_equal(mapping$reduced_index, c(1, 2, 3))
  expect_equal(mapping$treatment_name, c("A2", "A5", "A7"))
})

test_that("validate_active_set catches errors", {
  # Valid case
  result <- validate_active_set(c(1, 3, 5), m = 10)
  expect_true(result$valid)
  expect_equal(result$n_active, 3)
  expect_equal(length(result$issues), 0)

  # Empty set
  result <- validate_active_set(integer(0), m = 10)
  expect_false(result$valid)
  expect_true(any(grepl("empty", result$issues, ignore.case = TRUE)))

  # Out of range
  result <- validate_active_set(c(1, 11), m = 10)
  expect_false(result$valid)
  expect_true(any(grepl("range", result$issues, ignore.case = TRUE)))

  # Duplicates
  result <- validate_active_set(c(1, 3, 3), m = 10)
  expect_false(result$valid)
  expect_true(any(grepl("duplicate", result$issues, ignore.case = TRUE)))

  # All active (warning but valid)
  result <- validate_active_set(1:10, m = 10)
  expect_true(result$valid)
  expect_true(length(result$warnings) > 0)
})
