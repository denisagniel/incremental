library(testthat)
library(incremental)

test_that("generate_binary_vectors works correctly", {
  # Test m = 2
  vectors_2 <- generate_binary_vectors(2)
  expect_equal(nrow(vectors_2), 4)
  expect_equal(ncol(vectors_2), 2)
  # Correct order: (0,0), (1,0), (0,1), (1,1)
  expect_equal(vectors_2, matrix(c(0,0,1,0,0,1,1,1), nrow=4, byrow=TRUE))
  
  # Test m = 3
  vectors_3 <- generate_binary_vectors(3)
  expect_equal(nrow(vectors_3), 8)
  expect_equal(ncol(vectors_3), 3)
  
  # Test edge case
  vectors_1 <- generate_binary_vectors(1)
  expect_equal(nrow(vectors_1), 2)
  expect_equal(ncol(vectors_1), 1)
  expect_equal(vectors_1, matrix(c(0,1), ncol=1))
})

test_that("compute_marginals works correctly", {
  vectors <- generate_binary_vectors(2)

  # Test with known distribution (using vectors)
  # Order: (0,0), (1,0), (0,1), (1,1)
  # joint_dist[1]=0.3 is P(0,0), joint_dist[2]=0.2 is P(1,0),
  # joint_dist[3]=0.1 is P(0,1), joint_dist[4]=0.4 is P(1,1)
  joint_dist <- c(0.3, 0.2, 0.1, 0.4)
  marginals <- compute_marginals(joint_dist, vectors)

  # Manual calculation with correct ordering
  expected_marginals <- c(
    0.2 + 0.4,  # P(A1 = 1) = P(1,0) + P(1,1) = joint_dist[2] + joint_dist[4]
    0.1 + 0.4   # P(A2 = 1) = P(0,1) + P(1,1) = joint_dist[3] + joint_dist[4]
  )
  
  expect_equal(marginals, expected_marginals, tolerance = 1e-10)
  
  # Test with m parameter instead of vectors
  marginals2 <- compute_marginals(joint_dist, m = 2)
  expect_equal(marginals2, expected_marginals, tolerance = 1e-10)
  
  # Test that both methods give same result
  expect_equal(marginals, marginals2, tolerance = 1e-10)
})

test_that("find_active_constraints works correctly", {
  # Test with mixed probabilities
  p0_marginals <- c(0.3, 0.0, 0.7, 1.0)
  active <- find_active_constraints(p0_marginals)
  expect_equal(active, c(1, 3))
  
  # Test with all active
  p0_marginals <- c(0.3, 0.7, 0.5)
  active <- find_active_constraints(p0_marginals)
  expect_equal(active, c(1, 2, 3))
  
  # Test with none active
  p0_marginals <- c(0.0, 1.0)
  active <- find_active_constraints(p0_marginals)
  expect_equal(active, integer(0))
})
