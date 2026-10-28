library(testthat)
library(incremental)

test_that("simulate_vector_treatment_data works correctly", {
  # Test basic functionality
  sim_data <- simulate_vector_treatment_data(n = 100, m = 2, p = 2, seed = 123)
  
  # Check data structure
  expect_equal(nrow(sim_data$data), 100)
  expect_true("X1" %in% names(sim_data$data))
  expect_true("X2" %in% names(sim_data$data))
  expect_true("A1" %in% names(sim_data$data))
  expect_true("A2" %in% names(sim_data$data))
  expect_true("Y" %in% names(sim_data$data))
  expect_true("R" %in% names(sim_data$data))
  
  # Check treatment values are binary
  expect_true(all(sim_data$data$A1 %in% c(0, 1)))
  expect_true(all(sim_data$data$A2 %in% c(0, 1)))
  
  # Check outcome is numeric
  expect_true(is.numeric(sim_data$data$Y))
  
  # Check missing indicator
  expect_true(all(sim_data$data$R %in% c(0, 1)))
})


test_that("data generation with missing outcomes works", {
  sim_data <- simulate_vector_treatment_data(n = 100, m = 2, p = 2, 
                                           missing_prob = 0.2, seed = 123)
  
  # Check that some outcomes are missing
  expect_true(any(sim_data$data$R == 0))
  expect_true(any(sim_data$data$R == 1))
  
  # Check that missing probability is approximately correct
  missing_rate <- mean(sim_data$data$R == 0)
  expect_true(abs(missing_rate - 0.2) < 0.1)
})
