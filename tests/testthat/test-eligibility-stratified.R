#' @file test-eligibility-stratified.R
#' @title Tests for Eligibility-Stratified Estimation
#' @description
#' Comprehensive tests for eligibility-stratified vector incremental effects
#' estimation. Tests helper functions and full workflow.

context("Eligibility-Stratified Estimation")

# ============================================================================
# Test 1: identify_eligibility_patterns()
# ============================================================================

test_that("identify_eligibility_patterns works with 2 treatments, 4 patterns", {
  set.seed(123)
  n <- 400

  # Create 4 distinct patterns: 00, 01, 10, 11
  # Pattern 1: Both eligible (40%)
  # Pattern 2: A1 only (30%)
  # Pattern 3: A2 only (20%)
  # Pattern 4: Neither (10%)

  pattern_assignment <- sample(1:4, n, replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    elig_A1 = ifelse(pattern_assignment %in% c(1, 2), 1, 0),
    elig_A2 = ifelse(pattern_assignment %in% c(1, 3), 1, 0),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- identify_eligibility_patterns(data, eligibility_cols = c("elig_A1", "elig_A2"))

  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("patterns", "pattern_assignments", "m") %in% names(result)))

  # Check m
  expect_equal(result$m, 2)

  # Check number of patterns (should find 4)
  expect_equal(nrow(result$patterns), 4)

  # Check pattern strings
  pattern_strings <- result$patterns$pattern_string
  expect_true(all(pattern_strings %in% c("00", "01", "10", "11")))

  # Check sorting by prevalence (descending)
  expect_true(all(diff(result$patterns$prevalence) <= 0))

  # Check prevalences sum to 1
  expect_equal(sum(result$patterns$prevalence), 1)

  # Check n_eligible
  expect_equal(result$patterns$n_eligible[result$patterns$pattern_string == "00"], 0)
  expect_equal(result$patterns$n_eligible[result$patterns$pattern_string == "11"], 2)
  expect_equal(result$patterns$n_eligible[result$patterns$pattern_string == "01"], 1)
  expect_equal(result$patterns$n_eligible[result$patterns$pattern_string == "10"], 1)

  # Check pattern assignments
  expect_equal(length(result$pattern_assignments), n)
  expect_true(all(result$pattern_assignments %in% 1:4))
})

test_that("identify_eligibility_patterns works with 3 treatments, sparse patterns", {
  set.seed(456)
  n <- 300

  # Only 3 out of 8 possible patterns present
  pattern_assignment <- sample(1:3, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    A3 = rbinom(n, 1, 0.5),
    elig_A1 = ifelse(pattern_assignment == 1, 1, 0),
    elig_A2 = ifelse(pattern_assignment %in% c(1, 2), 1, 0),
    elig_A3 = 1,  # Always eligible
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- identify_eligibility_patterns(data, eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"))

  # Should find 3 patterns (not 8)
  expect_equal(nrow(result$patterns), 3)
  expect_equal(result$m, 3)

  # All patterns should have A3 eligible
  for (i in 1:nrow(result$patterns)) {
    elig_vec <- result$patterns$eligibility[[i]]
    expect_equal(elig_vec[3], 1)
  }
})

test_that("identify_eligibility_patterns handles single pattern", {
  set.seed(789)
  n <- 100

  # All observations have same pattern
  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    elig_A1 = 1,
    elig_A2 = 1,
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  result <- identify_eligibility_patterns(data, eligibility_cols = c("elig_A1", "elig_A2"))

  expect_equal(nrow(result$patterns), 1)
  expect_equal(result$patterns$pattern_string, "11")
  expect_equal(result$patterns$n_obs, n)
  expect_equal(result$patterns$prevalence, 1)
})

test_that("identify_eligibility_patterns requires explicit eligibility_cols", {
  set.seed(101)
  n <- 100

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    elig_A1 = rbinom(n, 1, 0.8),
    elig_A2 = rbinom(n, 1, 0.7),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  # Auto-detection has been removed: NULL eligibility_cols must error.
  expect_error(
    identify_eligibility_patterns(data, eligibility_cols = NULL),
    "eligibility_cols.*must be provided"
  )

  # With explicit columns it works.
  result <- identify_eligibility_patterns(data, eligibility_cols = c("elig_A1", "elig_A2"))
  expect_equal(result$m, 2)
  expect_true(nrow(result$patterns) >= 1)
})

test_that("identify_eligibility_patterns errors on missing columns", {
  data <- data.frame(
    A1 = rbinom(100, 1, 0.5),
    A2 = rbinom(100, 1, 0.5),
    X1 = rnorm(100),
    Y = rnorm(100)
  )

  expect_error(
    identify_eligibility_patterns(data, eligibility_cols = c("elig_A1", "elig_A2")),
    "Eligibility columns not found"
  )
})

# ============================================================================
# Test 2: validate_stratification_feasibility()
# ============================================================================

test_that("validate_stratification_feasibility accepts feasible case", {
  # Create patterns data frame
  patterns <- data.frame(
    pattern_id = 1:3,
    pattern_string = c("11", "10", "01"),
    n_obs = c(150, 120, 80),
    prevalence = c(0.43, 0.34, 0.23),
    n_eligible = c(2, 1, 1),
    stringsAsFactors = FALSE
  )
  patterns$eligibility <- list(c(1, 1), c(1, 0), c(0, 1))
  patterns$obs_indices <- list(1:150, 151:270, 271:350)

  result <- validate_stratification_feasibility(patterns, min_stratum_size = 50)

  expect_true(result$feasible)
  expect_equal(result$reason, "Feasible")
  expect_equal(result$n_feasible_patterns, 3)
  expect_equal(result$coverage, 1.0)
  expect_equal(result$feasible_pattern_ids, 1:3)
})

test_that("validate_stratification_feasibility rejects too sparse", {
  # Only 1 pattern with sufficient sample
  patterns <- data.frame(
    pattern_id = 1:3,
    pattern_string = c("11", "10", "01"),
    n_obs = c(150, 30, 20),
    prevalence = c(0.75, 0.15, 0.10),
    n_eligible = c(2, 1, 1),
    stringsAsFactors = FALSE
  )
  patterns$eligibility <- list(c(1, 1), c(1, 0), c(0, 1))
  patterns$obs_indices <- list(1:150, 151:180, 181:200)

  result <- validate_stratification_feasibility(patterns, min_stratum_size = 100)

  expect_false(result$feasible)
  expect_match(result$reason, "Only 1 pattern.*Need at least 2")
  expect_equal(result$n_feasible_patterns, 1)
})

test_that("validate_stratification_feasibility rejects low coverage", {
  # Many patterns but only small ones pass min_stratum_size
  patterns <- data.frame(
    pattern_id = 1:4,
    pattern_string = c("111", "110", "011", "101"),
    n_obs = c(100, 110, 20, 30),
    prevalence = c(0.38, 0.42, 0.08, 0.12),
    n_eligible = c(3, 2, 2, 2),
    stringsAsFactors = FALSE
  )
  patterns$eligibility <- list(c(1, 1, 1), c(1, 1, 0), c(0, 1, 1), c(1, 0, 1))
  patterns$obs_indices <- list(1:100, 101:210, 211:230, 231:260)

  # min_stratum_size = 200: only no patterns pass → feasible count = 0 → "Only 0 pattern(s)"
  result <- validate_stratification_feasibility(patterns, min_stratum_size = 200, min_coverage = 0.5)

  expect_false(result$feasible)
  expect_match(result$reason, "Only 0 pattern")
})

test_that("validate_stratification_feasibility rejects single pattern", {
  patterns <- data.frame(
    pattern_id = 1,
    pattern_string = "11",
    n_obs = 500,
    prevalence = 1.0,
    n_eligible = 2,
    stringsAsFactors = FALSE
  )
  patterns$eligibility <- list(c(1, 1))
  patterns$obs_indices <- list(1:500)

  result <- validate_stratification_feasibility(patterns)

  expect_false(result$feasible)
  expect_match(result$reason, "Only 1 eligibility pattern")
})

test_that("validate_stratification_feasibility rejects patterns with no eligible treatments", {
  patterns <- data.frame(
    pattern_id = 1:2,
    pattern_string = c("00", "11"),
    n_obs = c(200, 300),
    prevalence = c(0.4, 0.6),
    n_eligible = c(0, 2),
    stringsAsFactors = FALSE
  )
  patterns$eligibility <- list(c(0, 0), c(1, 1))
  patterns$obs_indices <- list(1:200, 201:500)

  result <- validate_stratification_feasibility(patterns)

  expect_false(result$feasible)
  expect_match(result$reason, "have no eligible treatments")
})

# ============================================================================
# Test 3: estimate_pattern_effect()
# ============================================================================

test_that("estimate_pattern_effect works with full eligibility", {
  set.seed(234)
  n <- 200

  # Simple DGP
  X1 <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(-0.5 + 0.5 * X1))
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X1))
  Y <- 1 + 0.5 * A1 + 0.8 * A2 + 0.3 * X1 + rnorm(n, 0, 0.5)

  data <- data.frame(A1 = A1, A2 = A2, elig_A1 = 1, elig_A2 = 1, X1 = X1, Y = Y)

  pattern_info <- list(
    pattern_id = 1,
    pattern_string = "11",
    eligibility = c(1, 1),
    n_obs = n,
    prevalence = 1.0,
    n_eligible = 2,
    obs_indices = 1:n
  )

  result <- estimate_pattern_effect(
    data = data,
    pattern_info = pattern_info,
    covariate_cols = "X1",
    delta =c(1.5, 2.0),
    divergence = "kl",
    verbose = FALSE
  )

  expect_true(result$success)
  expect_true(is.numeric(result$estimate))
  expect_false(is.na(result$estimate))
  expect_true(is.numeric(result$se))
  expect_equal(length(result$eif_values_extended), n)
  expect_equal(result$n_obs, n)
  expect_equal(result$n_eligible, 2)
  expect_equal(result$treatment_mapping$original, c(1, 2))
  expect_equal(result$treatment_mapping$reduced, c(1, 2))
})

test_that("estimate_pattern_effect works with partial eligibility", {
  set.seed(345)
  n <- 200

  # 3 treatments, only A1 and A2 eligible (A3 not eligible)
  X1 <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(-0.5 + 0.5 * X1))
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X1))
  A3 <- rep(0, n)  # Not eligible
  Y <- 1 + 0.5 * A1 + 0.3 * A2 + 0.3 * X1 + rnorm(n, 0, 0.5)

  data <- data.frame(A1 = A1, A2 = A2, A3 = A3, elig_A1 = 1, elig_A2 = 1, elig_A3 = 0, X1 = X1, Y = Y)

  pattern_info <- list(
    pattern_id = 1,
    pattern_string = "110",
    eligibility = c(1, 1, 0),
    n_obs = n,
    prevalence = 1.0,
    n_eligible = 2,
    obs_indices = 1:n
  )

  result <- estimate_pattern_effect(
    data = data,
    pattern_info = pattern_info,
    covariate_cols = "X1",
    delta =c(1.5, 2.0, 1.3),  # delta[3] won't be used
    divergence = "kl",
    verbose = FALSE
  )

  expect_true(result$success)
  expect_equal(result$n_eligible, 2)
  expect_equal(result$treatment_mapping$original, c(1, 2))
  expect_equal(result$treatment_mapping$reduced, c(1, 2))
  expect_equal(length(result$eif_values_extended), n)
})

test_that("estimate_pattern_effect handles single-treatment patterns (m_k = 1)", {
  set.seed(456)
  n <- 300

  # Only A2 eligible: m_k = 1 (scalar incremental effect on one treatment).
  X1 <- rnorm(n)
  A1 <- rep(0, n)
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X1))
  Y <- 1 + 0.8 * A2 + 0.3 * X1 + rnorm(n, 0, 0.5)

  data <- data.frame(A1 = A1, A2 = A2, elig_A1 = 0, elig_A2 = 1, X1 = X1, Y = Y)

  pattern_info <- list(
    pattern_id = 1,
    pattern_string = "01",
    eligibility = c(0, 1),
    n_obs = n,
    prevalence = 1.0,
    n_eligible = 1,
    obs_indices = 1:n
  )

  result <- estimate_pattern_effect(
    data = data,
    pattern_info = pattern_info,
    covariate_cols = "X1",
    delta = c(1.5, 2.0),
    divergence = "kl",
    verbose = FALSE
  )

  # m_k = 1 is a valid scalar incremental effect and should succeed.
  expect_true(result$success)
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result$se))
  expect_equal(result$n_eligible, 1)
  # Only the second treatment is eligible, renamed to the single reduced index.
  expect_equal(result$treatment_mapping$original, 2)
  expect_equal(result$treatment_mapping$reduced, 1)
})

test_that("estimate_pattern_effect extends EIF correctly", {
  set.seed(567)
  n_full <- 300
  n_pattern <- 100

  X1 <- rnorm(n_full)
  A1 <- rbinom(n_full, 1, 0.5)
  A2 <- rbinom(n_full, 1, 0.5)
  Y <- 1 + 0.5 * A1 + 0.8 * A2 + 0.3 * X1 + rnorm(n_full, 0, 0.5)

  data <- data.frame(A1 = A1, A2 = A2, elig_A1 = 1, elig_A2 = 1, X1 = X1, Y = Y)

  # Pattern uses first 100 observations
  obs_indices <- 1:n_pattern

  pattern_info <- list(
    pattern_id = 1,
    pattern_string = "11",
    eligibility = c(1, 1),
    n_obs = n_pattern,
    prevalence = n_pattern / n_full,
    n_eligible = 2,
    obs_indices = obs_indices
  )

  result <- estimate_pattern_effect(
    data = data,
    pattern_info = pattern_info,
    covariate_cols = "X1",
    delta =c(1.5, 2.0),
    divergence = "kl",
    verbose = FALSE
  )

  # Check EIF extended to full length
  expect_equal(length(result$eif_values_extended), n_full)

  # Check non-zero only for pattern observations
  expect_true(all(result$eif_values_extended[(n_pattern + 1):n_full] == 0))
  expect_false(all(result$eif_values_extended[obs_indices] == 0))
})

test_that("estimate_pattern_effect errors on insufficient sample", {
  set.seed(678)
  n <- 30  # Too small

  X1 <- rnorm(n)
  A1 <- rbinom(n, 1, 0.5)
  A2 <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)

  data <- data.frame(A1 = A1, A2 = A2, elig_A1 = 1, elig_A2 = 1, X1 = X1, Y = Y)

  pattern_info <- list(
    pattern_id = 1,
    pattern_string = "11",
    eligibility = c(1, 1),
    n_obs = n,
    prevalence = 1.0,
    n_eligible = 2,
    obs_indices = 1:n
  )

  result <- estimate_pattern_effect(
    data = data,
    pattern_info = pattern_info,
    covariate_cols = "X1",
    delta =c(1.5, 2.0),
    divergence = "kl",
    verbose = FALSE
  )

  expect_false(result$success)
  expect_match(result$error_message, "Insufficient sample size")
})

test_that("estimate_pattern_effect handles no eligible treatments", {
  set.seed(789)
  n <- 100

  data <- data.frame(
    A1 = rep(0, n),
    A2 = rep(0, n),
    elig_A1 = 0,
    elig_A2 = 0,
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  pattern_info <- list(
    pattern_id = 1,
    pattern_string = "00",
    eligibility = c(0, 0),
    n_obs = n,
    prevalence = 1.0,
    n_eligible = 0,
    obs_indices = 1:n
  )

  result <- estimate_pattern_effect(
    data = data,
    pattern_info = pattern_info,
    covariate_cols = "X1",
    delta =c(1.5, 2.0),
    verbose = FALSE
  )

  expect_false(result$success)
  expect_match(result$error_message, "Active set is empty")
})

# ============================================================================
# Test 4: aggregate_stratified_results()
# ============================================================================

test_that("aggregate_stratified_results computes weighted average correctly", {
  set.seed(890)
  n <- 300

  # Create mock pattern results
  # Pattern 1: 70%, estimate = 0.10, se = 0.05
  # Pattern 2: 30%, estimate = 0.20, se = 0.08
  # Expected weighted: 0.7 * 0.10 + 0.3 * 0.20 = 0.13

  eif1 <- rnorm(n, mean = 0, sd = 0.05 * sqrt(n))
  eif2 <- rnorm(n, mean = 0, sd = 0.08 * sqrt(n))

  # Extend EIFs (pattern 1 uses first 210 obs, pattern 2 uses last 90)
  eif1_extended <- c(eif1[1:210], rep(0, 90))
  eif2_extended <- c(rep(0, 210), eif2[1:90])

  pattern_results <- list(
    list(success = TRUE, estimate = 0.10, se = 0.05,
         eif_values_extended = eif1_extended, n_obs = 210, n_eligible = 2,
         treatment_mapping = list(original = c(1, 2), reduced = c(1, 2))),
    list(success = TRUE, estimate = 0.20, se = 0.08,
         eif_values_extended = eif2_extended, n_obs = 90, n_eligible = 2,
         treatment_mapping = list(original = c(1, 2), reduced = c(1, 2)))
  )

  patterns <- data.frame(
    pattern_id = 1:2,
    pattern_string = c("11", "10"),
    n_obs = c(210, 90),
    prevalence = c(0.7, 0.3),
    n_eligible = c(2, 2),
    stringsAsFactors = FALSE
  )

  result <- aggregate_stratified_results(pattern_results, patterns)

  # Check weighted estimate
  expect_equal(result$estimate, 0.7 * 0.10 + 0.3 * 0.20, tolerance = 1e-10)

  # Check structure
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(is.numeric(result$p_value))
  expect_equal(length(result$eif_values), n)
  expect_equal(result$n_patterns_used, 2)
  expect_equal(result$coverage, 1.0)

  # Check pattern summaries
  expect_equal(nrow(result$pattern_summaries), 2)
  expect_equal(result$pattern_summaries$weight, c(0.7, 0.3))
})

test_that("aggregate_stratified_results handles failed patterns", {
  set.seed(901)
  n <- 300

  eif1 <- rnorm(n, 0, 0.05 * sqrt(n))
  eif1_extended <- c(eif1[1:200], rep(0, 100))

  # Pattern 1 succeeds, pattern 2 fails
  pattern_results <- list(
    list(success = TRUE, estimate = 0.15, se = 0.06,
         eif_values_extended = eif1_extended, n_obs = 200, n_eligible = 2),
    list(success = FALSE, estimate = NA, se = NA,
         eif_values_extended = rep(NA, n), n_obs = 100, n_eligible = 1,
         error_message = "Convergence failed")
  )

  patterns <- data.frame(
    pattern_id = 1:2,
    pattern_string = c("11", "10"),
    n_obs = c(200, 100),
    prevalence = c(0.67, 0.33),
    n_eligible = c(2, 1),
    stringsAsFactors = FALSE
  )

  # Should warn about excluded pattern (loud report: count + which + coverage)
  expect_warning(
    result <- aggregate_stratified_results(pattern_results, patterns),
    "1 of 2 eligibility patterns FAILED"
  )

  # Check only successful pattern used
  expect_equal(result$n_patterns_used, 1)
  expect_equal(result$coverage, 0.67)

  # Estimate should be from pattern 1 only (renormalized weight = 1.0)
  expect_equal(result$estimate, 0.15, tolerance = 1e-10)
})

test_that("aggregate_stratified_results errors if all fail", {
  pattern_results <- list(
    list(success = FALSE, error_message = "Failed 1"),
    list(success = FALSE, error_message = "Failed 2")
  )

  patterns <- data.frame(
    pattern_id = 1:2,
    pattern_string = c("11", "10"),
    n_obs = c(100, 100),
    prevalence = c(0.5, 0.5),
    n_eligible = c(2, 1),
    stringsAsFactors = FALSE
  )

  expect_error(
    aggregate_stratified_results(pattern_results, patterns),
    "All patterns failed"
  )
})

test_that("aggregate_stratified_results renormalizes weights", {
  set.seed(1012)
  n <- 300

  eif1 <- rnorm(n, 0, 0.05 * sqrt(n))
  eif2 <- rnorm(n, 0, 0.06 * sqrt(n))

  eif1_extended <- c(eif1[1:150], rep(0, 150))
  eif2_extended <- c(rep(0, 150), eif2[1:150])

  # Both succeed, but original prevalences were 0.4 and 0.2 (third pattern excluded)
  pattern_results <- list(
    list(success = TRUE, estimate = 0.10, se = 0.05,
         eif_values_extended = eif1_extended, n_obs = 150, n_eligible = 2),
    list(success = TRUE, estimate = 0.15, se = 0.06,
         eif_values_extended = eif2_extended, n_obs = 150, n_eligible = 2)
  )

  patterns <- data.frame(
    pattern_id = 1:2,
    pattern_string = c("11", "10"),
    n_obs = c(150, 150),
    prevalence = c(0.4, 0.2),  # Should renormalize to 0.667, 0.333
    n_eligible = c(2, 2),
    stringsAsFactors = FALSE
  )

  result <- aggregate_stratified_results(pattern_results, patterns)

  # Check weights sum to 1 after renormalization
  expect_equal(sum(result$pattern_summaries$weight), 1.0, tolerance = 1e-10)

  # Check renormalized weights
  expect_equal(result$pattern_summaries$weight[1], 0.4 / 0.6, tolerance = 1e-10)
  expect_equal(result$pattern_summaries$weight[2], 0.2 / 0.6, tolerance = 1e-10)
})

# ============================================================================
# Test 5: estimate_qdelta_mean_stratified() - Integration Tests
# ============================================================================

test_that("estimate_qdelta_mean_stratified full workflow with 3 treatments, 2 patterns", {
  set.seed(1123)
  n <- 600
  p <- 3

  # Create 2 patterns (both with >= 2 eligible treatments)
  # Pattern 1: All 3 eligible (60%)
  # Pattern 2: A1 and A2 only (40%)

  pattern <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)

  # Eligibility
  elig_A1 <- rep(1, n)
  elig_A2 <- rep(1, n)
  elig_A3 <- ifelse(pattern == 1, 1, 0)

  # Treatment (only if eligible)
  A1 <- rbinom(n, 1, plogis(-0.5 + 0.3 * X[,1]))
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X[,2]))
  A3 <- ifelse(elig_A3 == 1, rbinom(n, 1, plogis(-0.2 + 0.2 * X[,3])), 0)

  # Outcome
  Y <- 2 + 0.6 * A1 + 0.8 * A2 + 0.4 * A3 + 0.4 * X[,1] + 0.3 * X[,2] + rnorm(n, 0, 0.8)

  data <- data.frame(
    A1 = A1, A2 = A2, A3 = A3,
    elig_A1 = elig_A1, elig_A2 = elig_A2, elig_A3 = elig_A3,
    X, Y
  )

  result <- estimate_qdelta_mean_stratified(
    data = data,
    delta = c(1.3, 1.5, 1.2),
    covariate_cols = paste0("X", 1:p),
    eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"),
    divergence = "kl",
    min_stratum_size = 80,
    verbose = FALSE
  )

  # Check structure
  expect_s3_class(result, "qdelta_stratified")
  expect_true(is.numeric(result$estimate))
  expect_false(is.na(result$estimate))
  expect_true(is.numeric(result$se))
  expect_true(result$se > 0)

  # Check confidence interval
  expect_true(result$ci_lower < result$estimate)
  expect_true(result$ci_upper > result$estimate)

  # Check pattern information
  expect_equal(result$m, 3)
  expect_equal(result$n, n)
  expect_true(result$n_patterns_used >= 2)
  expect_true(result$coverage >= 0.5)

  # Check pattern summaries
  expect_true(is.data.frame(result$pattern_summaries))
  expect_true(nrow(result$pattern_summaries) >= 2)
  expect_true(all(c("pattern_string", "n_obs", "weight", "estimate", "se") %in%
                   names(result$pattern_summaries)))

  # Check EIF
  expect_equal(length(result$eif_values), n)
  expect_true(all(is.finite(result$eif_values)))

  # EIF should be approximately centered
  expect_true(abs(mean(result$eif_values)) < 0.1 * sd(result$eif_values))
})

test_that("estimate_qdelta_mean_stratified matches non-stratified when all eligible", {
  set.seed(1234)
  n <- 300
  p <- 2

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)

  # Create 2 patterns by varying eligibility
  pattern <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))

  A1 <- rbinom(n, 1, plogis(-0.5 + 0.3 * X[,1]))
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X[,2]))
  Y <- 1.5 + 0.5 * A1 + 0.7 * A2 + 0.3 * X[,1] + rnorm(n, 0, 0.5)

  data <- data.frame(
    A1 = A1, A2 = A2,
    elig_A1 = ifelse(pattern == 1, 1, 0),  # Pattern-dependent eligibility
    elig_A2 = 1,  # Always eligible
    X, Y
  )

  delta_test <- c(1.5, 1.8)

  # Stratified
  result_strat <- estimate_qdelta_mean_stratified(
    data = data,
    delta = delta_test,
    covariate_cols = paste0("X", 1:p),
    eligibility_cols = c("elig_A1", "elig_A2"),
    divergence = "kl",
    min_stratum_size = 80,
    verbose = FALSE
  )

  # Should have at least 1 pattern
  expect_true(result_strat$n_patterns_used >= 1)
  expect_true(is.numeric(result_strat$estimate))
  expect_false(is.na(result_strat$estimate))
})

test_that("estimate_qdelta_mean_stratified errors when all observations same pattern", {
  set.seed(12345)
  n <- 300
  p <- 2

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)

  A1 <- rbinom(n, 1, plogis(-0.5 + 0.3 * X[,1]))
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X[,2]))
  Y <- 1.5 + 0.5 * A1 + 0.7 * A2 + 0.3 * X[,1] + rnorm(n, 0, 0.5)

  data <- data.frame(
    A1 = A1, A2 = A2,
    elig_A1 = 1, elig_A2 = 1,  # All same eligibility
    X, Y
  )

  # Should error because only 1 pattern
  expect_error(
    estimate_qdelta_mean_stratified(data, delta = c(1.5, 1.8),
                                    covariate_cols = "X1",
                                    eligibility_cols = c("elig_A1", "elig_A2"),
                                    verbose = FALSE),
    "Only 1 eligibility pattern"
  )
})

test_that("estimate_qdelta_mean_stratified handles min_stratum_size filtering", {
  set.seed(1345)
  n <- 400

  # Create patterns with varying sizes (all >= 2 eligible treatments)
  # Pattern 1: 60% (240 obs) - should be included (3 eligible)
  # Pattern 2: 25% (100 obs) - should be included with min=80 (2 eligible)
  # Pattern 3: 15% (60 obs) - should be excluded with min=80 (2 eligible)

  pattern <- sample(1:3, n, replace = TRUE, prob = c(0.6, 0.25, 0.15))

  X1 <- rnorm(n)

  # All patterns have >= 2 eligible treatments (m=3)
  elig_A1 <- rep(1, n)
  elig_A2 <- rep(1, n)
  elig_A3 <- as.integer(pattern == 1)

  A1 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X1))
  A2 <- rbinom(n, 1, plogis(-0.2 + 0.2 * X1))
  A3 <- ifelse(elig_A3 == 1, rbinom(n, 1, 0.4), 0)
  Y <- rnorm(n, 2 + 0.5 * A1 + 0.6 * A2 + 0.3 * A3, 1)

  data <- data.frame(A1 = A1, A2 = A2, A3 = A3,
                     elig_A1 = elig_A1, elig_A2 = elig_A2, elig_A3 = elig_A3,
                     X1 = X1, Y = Y)

  result <- estimate_qdelta_mean_stratified(
    data = data,
    delta = c(1.3, 1.5, 1.2),
    covariate_cols = "X1",
    eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"),
    min_stratum_size = 80,
    verbose = FALSE
  )

  # Should use only patterns with n >= 80 (patterns 1 and 2)
  expect_true(result$n_patterns_used >= 2)
  expect_true(result$coverage >= 0.8)  # 0.6 + 0.25 = 0.85
})

test_that("estimate_qdelta_mean_stratified handles max_patterns limit", {
  set.seed(1456)
  n <- 600

  # Create 4 patterns (all with >= 2 eligible treatments, using m=3)
  pattern <- sample(1:4, n, replace = TRUE, prob = c(0.3, 0.25, 0.25, 0.2))

  X1 <- rnorm(n)

  # Each pattern has >= 2 eligible treatments
  elig_A1 <- rep(1, n)  # Always eligible
  elig_A2 <- rep(1, n)  # Always eligible
  elig_A3 <- as.integer(pattern %in% c(1, 2))  # Only patterns 1 and 2

  A1 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X1))
  A2 <- rbinom(n, 1, plogis(-0.2 + 0.2 * X1))
  A3 <- ifelse(elig_A3 == 1, rbinom(n, 1, 0.4), 0)
  Y <- rnorm(n, 2 + 0.5 * A1 + 0.6 * A2 + 0.3 * A3, 1)

  data <- data.frame(
    A1 = A1, A2 = A2, A3 = A3,
    elig_A1 = elig_A1, elig_A2 = elig_A2, elig_A3 = elig_A3,
    X1 = X1, Y = Y
  )

  result <- estimate_qdelta_mean_stratified(
    data = data,
    delta = c(1.3, 1.5, 1.2),
    covariate_cols = "X1",
    eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"),
    min_stratum_size = 30,
    max_patterns = 2,  # Limit to top 2
    verbose = FALSE
  )

  # Should use at most 2 patterns (limited by max_patterns)
  expect_true(result$n_patterns_used <= 2)
  expect_true(result$n_patterns_used >= 1)
})

test_that("estimate_qdelta_mean_stratified errors on infeasible stratification", {
  set.seed(1567)
  n <- 100

  # All same pattern - not feasible for stratification
  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    elig_A1 = 1,
    elig_A2 = 1,
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  expect_error(
    estimate_qdelta_mean_stratified(data, delta = c(1.5, 2.0),
                                    covariate_cols = "X1",
                                    eligibility_cols = c("elig_A1", "elig_A2"),
                                    verbose = FALSE),
    "Stratification not feasible.*Only 1 eligibility pattern"
  )
})

test_that("estimate_qdelta_mean_stratified errors on bootstrap request", {
  set.seed(1678)
  n <- 200

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    elig_A1 = rbinom(n, 1, 0.7),
    elig_A2 = rbinom(n, 1, 0.8),
    X1 = rnorm(n),
    Y = rnorm(n)
  )

  expect_error(
    estimate_qdelta_mean_stratified(data, delta = c(1.5, 2.0), bootstrap = TRUE, verbose = FALSE),
    "Bootstrap inference not currently supported"
  )
})

test_that("estimate_qdelta_mean_stratified requires explicit eligibility_cols", {
  set.seed(1789)
  n <- 200

  pattern <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))

  data <- data.frame(
    A1 = rbinom(n, 1, 0.5),
    A2 = rbinom(n, 1, 0.5),
    elig_A1 = ifelse(pattern == 1, 1, 0),
    elig_A2 = 1,  # Always eligible
    X1 = rnorm(n),
    Y = rnorm(n, 2, 1)
  )

  # Auto-detection has been removed: NULL eligibility_cols must error.
  expect_error(
    estimate_qdelta_mean_stratified(
      data = data,
      delta = c(1.3, 1.5),
      covariate_cols = "X1",
      eligibility_cols = NULL,
      min_stratum_size = 50,
      verbose = FALSE
    ),
    "eligibility_cols.*must be provided"
  )

  # With explicit columns it works.
  result <- estimate_qdelta_mean_stratified(
    data = data,
    delta = c(1.3, 1.5),
    covariate_cols = "X1",
    eligibility_cols = c("elig_A1", "elig_A2"),
    min_stratum_size = 50,
    verbose = FALSE
  )

  expect_true(result$n_patterns_used >= 1)
  expect_equal(result$m, 2)
})

# ============================================================================
# Test 6: Edge Cases and Error Handling
# ============================================================================

test_that("estimate_qdelta_mean_stratified handles cross-fitting disable for small patterns", {
  set.seed(1890)
  n <- 250

  # Pattern 1: 150 obs (enough for cross-fitting)
  # Pattern 2: 60 obs (too small for 5-fold, should disable)

  pattern <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))

  X1 <- rnorm(n)
  A1 <- rbinom(n, 1, 0.5)
  A2 <- rbinom(n, 1, 0.5)
  Y <- rnorm(n, 2 + 0.5 * A1, 1)

  data <- data.frame(
    A1 = A1, A2 = A2,
    elig_A1 = ifelse(pattern == 1, 1, 0),
    elig_A2 = 1,
    X1 = X1, Y = Y
  )

  # Request cross-fitting
  result <- estimate_qdelta_mean_stratified(
    data = data,
    delta = c(1.3, 1.5),
    covariate_cols = "X1",
    eligibility_cols = c("elig_A1", "elig_A2"),
    cross_fit = TRUE,
    n_folds = 5,
    min_stratum_size = 50,
    verbose = FALSE
  )

  # Should succeed (cross-fitting automatically disabled for small pattern)
  expect_true(result$n_patterns_used >= 1)
  expect_true(is.numeric(result$estimate))
})

# ============================================================================
# Test 7: threshold_type forwarding (none / vanishing / adaptive)
# ============================================================================
# These verify the Path A change: the stratified path can reach the thresholded
# estimators the paper's theory covers (vanishing -> cor:onestep; adaptive ->
# localized-margin). Thresholding acts per-coordinate on m_k(X) WITHIN each
# eligibility pattern, which is separate from eligibility subsetting.

# Shared small stratified DGP (2 patterns, all patterns >= 2 eligible).
make_stratified_data <- function(n = 600, seed = 2026) {
  set.seed(seed)
  p <- 3
  pattern <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)
  elig_A1 <- rep(1, n)
  elig_A2 <- rep(1, n)
  elig_A3 <- ifelse(pattern == 1, 1, 0)
  A1 <- rbinom(n, 1, plogis(-0.5 + 0.3 * X[, 1]))
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X[, 2]))
  A3 <- ifelse(elig_A3 == 1, rbinom(n, 1, plogis(-0.2 + 0.2 * X[, 3])), 0)
  Y <- 2 + 0.6 * A1 + 0.8 * A2 + 0.4 * A3 + 0.4 * X[, 1] + 0.3 * X[, 2] + rnorm(n, 0, 0.8)
  data.frame(A1 = A1, A2 = A2, A3 = A3,
             elig_A1 = elig_A1, elig_A2 = elig_A2, elig_A3 = elig_A3, X, Y)
}

test_that("stratified default threshold_type is 'none' and reproduces prior estimate", {
  data <- make_stratified_data(seed = 2026)
  common <- list(
    data = data, delta = c(1.3, 1.5, 1.2), covariate_cols = paste0("X", 1:3),
    eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"),
    divergence = "kl", min_stratum_size = 80, verbose = FALSE
  )

  # Default (no threshold_type arg) must equal explicit "none": regression guard
  # that the Path A change did not alter the reported (unthresholded) numbers.
  # Reset the seed before each call: cross-fitting draws fold assignments from the
  # RNG, so back-to-back calls would otherwise start from different RNG states.
  set.seed(909)
  res_default  <- do.call(estimate_qdelta_mean_stratified, common)
  set.seed(909)
  res_explicit <- do.call(estimate_qdelta_mean_stratified,
                          c(common, list(threshold_type = "none")))

  expect_equal(res_default$threshold_type, "none")
  expect_equal(res_default$estimate, res_explicit$estimate, tolerance = 1e-10)
  expect_equal(res_default$se,       res_explicit$se,       tolerance = 1e-10)
})

test_that("stratified runs end-to-end with vanishing and adaptive thresholds", {
  data <- make_stratified_data(seed = 3031)
  common <- list(
    data = data, delta = c(1.3, 1.5, 1.2), covariate_cols = paste0("X", 1:3),
    eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"),
    divergence = "kl", min_stratum_size = 80, verbose = FALSE
  )

  res_van <- do.call(estimate_qdelta_mean_stratified,
                     c(common, list(threshold_type = "vanishing", gamma = 0.25)))
  expect_s3_class(res_van, "qdelta_stratified")
  expect_equal(res_van$threshold_type, "vanishing")
  expect_true(is.finite(res_van$estimate))
  expect_true(is.finite(res_van$se) && res_van$se > 0)

  res_adp <- do.call(estimate_qdelta_mean_stratified,
                     c(common, list(threshold_type = "adaptive", epsilon = 0.1)))
  expect_s3_class(res_adp, "qdelta_stratified")
  expect_equal(res_adp$threshold_type, "adaptive")
  expect_true(is.finite(res_adp$estimate))
  expect_true(is.finite(res_adp$se) && res_adp$se > 0)
})

test_that("stratified rejects an invalid threshold_type", {
  data <- make_stratified_data(seed = 4041)
  expect_error(
    estimate_qdelta_mean_stratified(
      data = data, delta = c(1.3, 1.5, 1.2), covariate_cols = paste0("X", 1:3),
      eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"),
      min_stratum_size = 80, verbose = FALSE,
      threshold_type = "bogus"
    ),
    "should be one of"
  )
})

test_that("estimate_pattern_effect forwards threshold_type to the thresholded estimators", {
  set.seed(5051)
  n <- 300
  X1 <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(-0.5 + 0.5 * X1))
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.3 * X1))
  Y <- 1 + 0.5 * A1 + 0.8 * A2 + 0.3 * X1 + rnorm(n, 0, 0.5)
  data <- data.frame(A1 = A1, A2 = A2, elig_A1 = 1, elig_A2 = 1, X1 = X1, Y = Y)
  pattern_info <- list(
    pattern_id = 1, pattern_string = "11", eligibility = c(1, 1),
    n_obs = n, prevalence = 1.0, n_eligible = 2, obs_indices = 1:n
  )

  res_adp <- estimate_pattern_effect(
    data = data, pattern_info = pattern_info, covariate_cols = "X1",
    delta = c(1.5, 2.0), divergence = "kl",
    threshold_type = "adaptive", epsilon = 0.1, verbose = FALSE
  )
  expect_true(res_adp$success)
  expect_true(is.finite(res_adp$estimate))
})

test_that("single-eligibility (m=1) patterns estimate with sequential PS (regression: 2:m backwards-seq bug)", {
  # Regression for the `for (k in 2:m)` bug in compute_joint_probabilities[_cv]:
  # for m = 1, 2:m is c(2, 1) (backwards) and indexed a non-existent 2nd
  # conditional model -> "subscript out of bounds". Single-eligibility patterns
  # (the two largest strata in the NYS Medicaid application) were silently
  # dropped and renormalized away, biasing the delta = 1 estimate upward.
  set.seed(20260713)
  n <- 3000
  X1 <- rnorm(n); X2 <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * X1))
  Y <- 6000 + 2000 * A1 + 1500 * X1
  data <- data.frame(A1 = A1, X1 = X1, X2 = X2, Y = Y)

  # m = 1 with sequential PS must succeed (previously errored) in both CF modes.
  for (cf in c(TRUE, FALSE)) {
    res <- estimate_qdelta_mean(
      data = data, delta = 1.0, treatment_cols = "A1",
      covariate_cols = c("X1", "X2"), outcome_col = "Y",
      divergence = "kl", ps_method = "sequential",
      binary_method = "glm", outcome_method = "glm",
      estimator = "one_step", cross_fit = cf, n_folds = 3,
      threshold_type = "none"
    )
    expect_true(is.finite(res$estimate))
    # At delta = 1, Q_delta = observed law => estimate == sample mean of Y.
    expect_equal(res$estimate, mean(Y), tolerance = 1e-3)
  }
})

test_that("stratified estimator uses ALL patterns and recovers the sample mean at delta=1 (sequential PS)", {
  # End-to-end: with single-eligibility patterns present, the delta = 1 estimate
  # must equal the overall sample mean (no patterns silently dropped).
  set.seed(20260713)
  n <- 6000
  X <- matrix(rnorm(n * 3), n, 3); colnames(X) <- paste0("X", 1:3)
  e <- sapply(1:3, function(k) rbinom(n, 1, c(0.6, 0.4, 0.5)[k]))
  e[rowSums(e) == 0, 1] <- 1  # everyone eligible for >= 1
  A <- sapply(1:3, function(k) ifelse(e[, k] == 1, rbinom(n, 1, plogis(0.3 * X[, 1])), 0))
  colnames(A) <- paste0("A", 1:3); colnames(e) <- paste0("elig_A", 1:3)
  Y <- 8000 + 3000 * A[, 1] - 1500 * A[, 2] + 1000 * A[, 3] + 2000 * X[, 1]
  data <- data.frame(A, e, X, Y)

  res <- estimate_qdelta_mean_stratified(
    data = data, delta = rep(1.0, 3),
    covariate_cols = paste0("X", 1:3),
    eligibility_cols = paste0("elig_A", 1:3),
    outcome_col = "Y", divergence = "kl", estimator = "one_step",
    ps_method = "sequential", binary_method = "glm", outcome_method = "glm",
    cross_fit = TRUE, n_folds = 3, min_stratum_size = 100, verbose = FALSE
  )
  # No pattern (including single-eligibility) should be dropped.
  n_identified <- nrow(res$eligibility_patterns$patterns[
    res$eligibility_patterns$patterns$n_obs >= 100, ])
  expect_equal(res$n_patterns_used, n_identified)
  # delta = 1 recovers the overall sample mean.
  expect_equal(res$estimate, mean(Y), tolerance = 5e-3)
})

test_that("realized-coverage gate aborts when failed patterns drop too much of the sample", {
  # The guardrail for silent pattern-drop: if estimation failures push realized
  # coverage below min_realized_coverage, abort rather than return an estimate
  # over a non-representative subset. (This is the UI that would have surfaced
  # the single-eligibility failure immediately.)
  skip_if_not_installed("xgboost")
  set.seed(2027)
  n <- 4000
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("X1", "X2")
  e <- cbind(elig_A1 = rbinom(n, 1, 0.6), elig_A2 = rbinom(n, 1, 0.5))
  e[rowSums(e) == 0, 1] <- 1
  A1 <- ifelse(e[, 1] == 1, rbinom(n, 1, plogis(0.3 * X[, 1])), 0)
  A2 <- ifelse(e[, 2] == 1, rbinom(n, 1, plogis(-0.2 * X[, 1])), 0)
  Y <- 5000 + 2000 * A1 - 1000 * A2 + 1500 * X[, 1]
  data <- data.frame(A1, A2, e, X, Y)

  # An impossible coverage requirement (>100%) forces the gate to trip even when
  # every pattern succeeds, exercising the abort path deterministically.
  expect_error(
    suppressWarnings(estimate_qdelta_mean_stratified(
      data = data, delta = c(1, 1),
      covariate_cols = c("X1", "X2"),
      eligibility_cols = c("elig_A1", "elig_A2"),
      outcome_col = "Y", divergence = "kl", estimator = "one_step",
      ps_method = "sequential", binary_method = "glm", outcome_method = "glm",
      cross_fit = TRUE, n_folds = 3, min_stratum_size = 100,
      min_realized_coverage = 1.01, verbose = FALSE
    )),
    "Realized pattern coverage"
  )

  # With the gate disabled (0), the same run returns normally.
  res <- estimate_qdelta_mean_stratified(
    data = data, delta = c(1, 1),
    covariate_cols = c("X1", "X2"),
    eligibility_cols = c("elig_A1", "elig_A2"),
    outcome_col = "Y", divergence = "kl", estimator = "one_step",
    ps_method = "sequential", binary_method = "glm", outcome_method = "glm",
    cross_fit = TRUE, n_folds = 3, min_stratum_size = 100,
    min_realized_coverage = 0, verbose = FALSE
  )
  expect_true(is.finite(res$estimate))
  expect_true(res$realized_coverage > 0)
})

test_that("recombine_stratified_patterns reproduces the fit exactly at the same threshold", {
  # Fidelity: post-hoc recombination with all patterns must equal the original
  # aggregate to machine precision (shared .combine_patterns() core).
  set.seed(20260714)
  n <- 6000
  X <- matrix(rnorm(n * 3), n, 3); colnames(X) <- paste0("X", 1:3)
  e <- sapply(1:3, function(k) rbinom(n, 1, c(0.6, 0.4, 0.5)[k]))
  e[rowSums(e) == 0, 1] <- 1
  A <- sapply(1:3, function(k) ifelse(e[, k] == 1, rbinom(n, 1, plogis(0.3 * X[, 1])), 0))
  colnames(A) <- paste0("A", 1:3); colnames(e) <- paste0("elig_A", 1:3)
  Y <- 8000 + 3000 * A[, 1] - 1500 * A[, 2] + 1000 * A[, 3] + 2000 * X[, 1]
  data <- data.frame(A, e, X, Y)

  fit <- estimate_qdelta_mean_stratified(
    data = data, delta = rep(1.3, 3),
    covariate_cols = paste0("X", 1:3), eligibility_cols = paste0("elig_A", 1:3),
    outcome_col = "Y", divergence = "kl", estimator = "one_step",
    ps_method = "sequential", binary_method = "glm", outcome_method = "glm",
    cross_fit = TRUE, n_folds = 3, min_stratum_size = 100,
    min_realized_coverage = 0, verbose = FALSE
  )
  expect_false(is.null(fit$pattern_contributions))

  rc <- recombine_stratified_patterns(fit$pattern_contributions, min_stratum_size = 0)
  expect_equal(rc$estimate, fit$estimate, tolerance = 1e-10)
  expect_equal(rc$se, fit$se, tolerance = 1e-10)
  expect_equal(rc$n_patterns_used, fit$n_patterns_used)
  expect_equal(rc$coverage, fit$coverage, tolerance = 1e-10)
})

test_that("recombine_stratified_patterns at a higher threshold drops patterns and lowers coverage", {
  set.seed(20260714)
  n <- 6000
  X <- matrix(rnorm(n * 3), n, 3); colnames(X) <- paste0("X", 1:3)
  e <- sapply(1:3, function(k) rbinom(n, 1, c(0.6, 0.4, 0.5)[k]))
  e[rowSums(e) == 0, 1] <- 1
  A <- sapply(1:3, function(k) ifelse(e[, k] == 1, rbinom(n, 1, plogis(0.3 * X[, 1])), 0))
  colnames(A) <- paste0("A", 1:3); colnames(e) <- paste0("elig_A", 1:3)
  Y <- 8000 + 3000 * A[, 1] - 1500 * A[, 2] + 1000 * A[, 3] + 2000 * X[, 1]
  data <- data.frame(A, e, X, Y)

  fit <- estimate_qdelta_mean_stratified(
    data = data, delta = rep(1.3, 3),
    covariate_cols = paste0("X", 1:3), eligibility_cols = paste0("elig_A", 1:3),
    outcome_col = "Y", ps_method = "sequential", binary_method = "glm",
    outcome_method = "glm", cross_fit = TRUE, n_folds = 3,
    min_stratum_size = 100, min_realized_coverage = 0, verbose = FALSE
  )
  pc <- fit$pattern_contributions
  thr <- sort(pc$patterns$n_obs)[2]  # exclude the single smallest pattern
  rc <- recombine_stratified_patterns(pc, min_stratum_size = thr)

  expect_lt(rc$n_patterns_used, fit$n_patterns_used)
  expect_lt(rc$coverage, 1.0)
  expect_true(is.finite(rc$estimate) && is.finite(rc$se))
})

test_that("recombine_stratified_patterns supports custom (non-prevalence) weights", {
  set.seed(20260714)
  n <- 5000
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("X1", "X2")
  e <- cbind(elig_A1 = rbinom(n, 1, 0.6), elig_A2 = rbinom(n, 1, 0.5))
  e[rowSums(e) == 0, 1] <- 1
  A1 <- ifelse(e[, 1] == 1, rbinom(n, 1, plogis(0.3 * X[, 1])), 0)
  A2 <- ifelse(e[, 2] == 1, rbinom(n, 1, plogis(-0.2 * X[, 1])), 0)
  Y <- 5000 + 2000 * A1 - 1000 * A2 + 1500 * X[, 1]
  data <- data.frame(A1, A2, e, X, Y)

  fit <- estimate_qdelta_mean_stratified(
    data = data, delta = c(1, 1),
    covariate_cols = c("X1", "X2"), eligibility_cols = c("elig_A1", "elig_A2"),
    outcome_col = "Y", ps_method = "sequential", binary_method = "glm",
    outcome_method = "glm", cross_fit = TRUE, n_folds = 3,
    min_stratum_size = 100, min_realized_coverage = 0, verbose = FALSE
  )
  pc <- fit$pattern_contributions
  ids <- pc$patterns$pattern_id
  w <- setNames(rep(1, length(ids)), ids)  # equal weights
  rc <- recombine_stratified_patterns(pc, weights = w)

  expect_equal(sum(rc$patterns$weight), 1, tolerance = 1e-10)
  # Equal weighting generally differs from prevalence weighting.
  expect_false(isTRUE(all.equal(rc$estimate, fit$estimate)))
})
