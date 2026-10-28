#' @file test-nuisance-cache.R
#' @title Tests for nuisance-model-fit caching in stratified estimation
#' @description
#' The stratified pipeline refits per-pattern propensity and outcome models for
#' every (cost type x delta) cell even though neither nuisance depends on delta
#' and the propensity does not depend on the cost type. `new_nuisance_cache()`
#' removes that redundancy. These tests pin the properties that make the
#' optimisation safe:
#'   (i)   cached and freshly-fitted runs are BIT-identical,
#'   (ii)  the cache is really used (hit counts across >= 3 deltas),
#'   (iii) a fingerprint mismatch is a hard error, never a silent misalignment,
#'   (iv)  fits are RNG-neutral to the caller's stream,
#'   (v)   passing a cache with a thresholded estimator errors,
#'   (vi)  in-place nuisance mutation (TMLE targeting) cannot corrupt the cache,
#'   (vii) fit counts scale as O(P) + O(P x cost types), not O(P x CT x deltas).

context("Nuisance-fit cache")

# ============================================================================
# Cache primitives
# ============================================================================

test_that("new_nuisance_cache builds a classed environment with empty stores", {
  cache <- new_nuisance_cache(100L, tags = c("a", "b"))

  expect_s3_class(cache, "nuisance_cache")
  expect_true(is.environment(cache))
  expect_equal(cache$seed_base, 100L)
  expect_equal(cache$tags, c("a", "b"))
  expect_true(cache$tags_fixed)
  expect_equal(length(cache$ps), 0)
  expect_equal(length(cache$outcome), 0)

  s <- nuisance_cache_stats(cache)
  expect_equal(s$n_ps_entries, 0L)
  expect_equal(s$n_outcome_entries, 0L)
  expect_equal(s$size_mb, 0)

  expect_error(new_nuisance_cache("x"), "seed_base")
  expect_error(new_nuisance_cache(1L, tags = c("a", "a")), "unique")
})

test_that("tag indices are position-based when tags are registered up front", {
  cache <- new_nuisance_cache(0L, tags = c("total", "inpatient", "mh"))

  # Deterministic regardless of which tag is seen first: a resumed sweep that
  # skips "total" must still seed "mh" identically.
  expect_equal(.cache_tag_index(cache, "mh"), 3L)
  expect_equal(.cache_tag_index(cache, "total"), 1L)
  expect_error(.cache_tag_index(cache, "typo"), "not one of the registered tags")

  # Without registration, tags are added on first use.
  dyn <- new_nuisance_cache(0L)
  expect_equal(.cache_tag_index(dyn, "b"), 1L)
  expect_equal(.cache_tag_index(dyn, "a"), 2L)
  expect_equal(.cache_tag_index(dyn, "b"), 1L)
})

test_that("per-fit seeds follow the documented namespaced formula", {
  cache <- new_nuisance_cache(20260706L, tags = c("total", "inpatient"))

  s1 <- .cache_seeds(cache, "total", 3L)
  expect_equal(s1$ps, 20260706L + 10000L + 3L)
  expect_equal(s1$outcome, 20260706L + 20000L + 1000L * 1L + 3L)

  s2 <- .cache_seeds(cache, "inpatient", 3L)
  # Propensity seed does NOT depend on the cost type; outcome seed does.
  expect_equal(s2$ps, s1$ps)
  expect_equal(s2$outcome, 20260706L + 20000L + 1000L * 2L + 3L)

  expect_error(.cache_seeds(cache, "total", NULL), "pattern_pos")
  expect_error(.cache_seeds(cache, "total", 1000L), "collide")
})

test_that(".clear_outcomes drops only outcome entries", {
  cache <- new_nuisance_cache(1L, tags = c("a", "b"))
  fp <- .nuisance_fingerprint(1:4, 2, list(method = "glm"))
  mat <- matrix(0.5, nrow = 4, ncol = 2)

  .cache_put(cache, "ps", NULL, "11", mat, fp)
  .cache_put(cache, "outcome", "a", "11", mat, fp)
  .cache_put(cache, "outcome", "b", "11", mat, fp)
  expect_equal(nuisance_cache_stats(cache)$n_outcome_entries, 2L)

  n_dropped <- .clear_outcomes(cache)
  expect_equal(n_dropped, 2L)
  expect_equal(nuisance_cache_stats(cache)$n_outcome_entries, 0L)
  # Propensity fits survive: they are valid for every cost type.
  expect_equal(nuisance_cache_stats(cache)$n_ps_entries, 1L)
  expect_true(!is.null(.cache_get(cache, "ps", NULL, "11", fp)))
})

test_that(".cache_put rejects shape/fingerprint disagreement", {
  cache <- new_nuisance_cache(1L)
  fp <- .nuisance_fingerprint(1:4, 2, list(method = "glm"))
  expect_error(.cache_put(cache, "ps", NULL, "11", matrix(0, 3, 2), fp), "disagrees")
  expect_error(.cache_put(cache, "ps", NULL, "11", 1:8, fp), "must be a matrix")
  expect_error(.cache_put(cache, "bogus", NULL, "11", matrix(0, 4, 2), fp), "'ps' or 'outcome'")
})

# ============================================================================
# (iii) Fingerprint mismatch is a hard error (unit level)
# ============================================================================

test_that("mismatched fingerprint on lookup is a hard, classed error", {
  cache <- new_nuisance_cache(1L)
  fp <- .nuisance_fingerprint(1:10, 4, list(method = "glm"))
  .cache_put(cache, "ps", NULL, "1100", matrix(0.25, nrow = 10, ncol = 4), fp)

  # Same key, different rows (one more observation) -> must NOT return a matrix.
  fp_rows <- .nuisance_fingerprint(1:11, 4, list(method = "glm"))
  expect_error(.cache_get(cache, "ps", NULL, "1100", fp_rows),
               class = "nuisance_cache_mismatch")

  # Same number of rows, different row identities (idx_sum differs).
  fp_shift <- .nuisance_fingerprint(2:11, 4, list(method = "glm"))
  expect_error(.cache_get(cache, "ps", NULL, "1100", fp_shift),
               class = "nuisance_cache_mismatch")

  # Same rows, different nuisance configuration.
  fp_config <- .nuisance_fingerprint(1:10, 4, list(method = "xgboost"))
  expect_error(.cache_get(cache, "ps", NULL, "1100", fp_config),
               class = "nuisance_cache_mismatch")

  # Exact match still returns the stored matrix.
  expect_equal(dim(.cache_get(cache, "ps", NULL, "1100", fp)), c(10L, 4L))
})

# ============================================================================
# (v) threshold_type guard
# ============================================================================

test_that("passing a cache with a thresholded estimator errors", {
  d <- make_cache_test_data(n = 200, n_cost_types = 1)
  cache <- new_nuisance_cache(1L, tags = d$cost_types)

  for (tt in c("vanishing", "adaptive")) {
    expect_error(
      run_cache_test_estimate(d, rep(1.1, 6), cache = cache, threshold_type = tt),
      "only supported with threshold_type = 'none'"
    )
  }

  # Same guard at the pattern level (in case a caller bypasses the wrapper).
  pat <- identify_eligibility_patterns(cache_test_dataset(d), d$eligibility_cols)
  row1 <- pat$patterns[1, ]
  pinfo <- list(pattern_id = row1$pattern_id, pattern_string = row1$pattern_string,
                eligibility = row1$eligibility[[1]], n_obs = row1$n_obs,
                prevalence = row1$prevalence, n_eligible = row1$n_eligible,
                obs_indices = row1$obs_indices[[1]])
  expect_error(
    estimate_pattern_effect(
      data = cache_test_dataset(d), pattern_info = pinfo, delta = rep(1.1, 6),
      covariate_cols = d$covariate_cols, outcome_col = "Y",
      threshold_type = "vanishing", verbose = FALSE,
      nuisance_cache = cache, nuisance_tag = d$cost_types[1], pattern_pos = 1L),
    "only supported with threshold_type = 'none'"
  )

  # pattern_pos is mandatory with a cache (it namespaces the seeds).
  expect_error(
    estimate_pattern_effect(
      data = cache_test_dataset(d), pattern_info = pinfo, delta = rep(1.1, 6),
      covariate_cols = d$covariate_cols, outcome_col = "Y", verbose = FALSE,
      nuisance_cache = cache, nuisance_tag = d$cost_types[1]),
    "pattern_pos"
  )
})

# ============================================================================
# (i) Bit-identical results: cache hit == fresh fit
# ============================================================================

test_that("cached and freshly-fitted runs are bit-identical", {
  d <- make_cache_test_data()
  delta <- rep(1.1, 6)

  # Run A: fresh cache -> every nuisance is fitted.
  cache_a <- new_nuisance_cache(777L, tags = d$cost_types)
  set.seed(11)
  res_a <- run_cache_test_estimate(d, delta, "ct1", cache_a)
  stats_a <- nuisance_cache_stats(cache_a)
  expect_equal(stats_a$n_ps_fit, res_a$n_patterns_used)
  expect_equal(stats_a$n_ps_hit, 0L)

  # Run B: same cache, same delta -> every nuisance is a hit, nothing refitted.
  res_b <- run_cache_test_estimate(d, delta, "ct1", cache_a)
  stats_b <- nuisance_cache_stats(cache_a)
  expect_equal(stats_b$n_ps_fit, stats_a$n_ps_fit)
  expect_equal(stats_b$n_outcome_fit, stats_a$n_outcome_fit)
  expect_true(stats_b$n_ps_hit > 0)

  expect_identical(res_a$estimate, res_b$estimate)
  expect_identical(res_a$se, res_b$se)
  expect_identical(res_a$eif_values, res_b$eif_values)

  # Run C: a DIFFERENT fresh cache under a different ambient RNG state. Because
  # each fit is seeded from (seed_base, kind, tag, pattern_pos), the fits -- and
  # therefore the estimates -- must be identical to run A. This is what makes a
  # crash-resumed sweep reproducible.
  cache_c <- new_nuisance_cache(777L, tags = d$cost_types)
  set.seed(99999)
  invisible(runif(37))  # arbitrary RNG consumption before the call
  res_c <- run_cache_test_estimate(d, delta, "ct1", cache_c)

  expect_identical(res_a$estimate, res_c$estimate)
  expect_identical(res_a$eif_values, res_c$eif_values)

  # A different seed_base must give genuinely different splits (guards against
  # the seeds being ignored, which would make the test above vacuous).
  cache_d <- new_nuisance_cache(778L, tags = d$cost_types)
  res_d <- run_cache_test_estimate(d, delta, "ct1", cache_d)
  expect_false(isTRUE(all.equal(res_a$estimate, res_d$estimate, tolerance = 1e-12)))
})

# ============================================================================
# (ii)/(vii) Hit counts and fit-count scaling
# ============================================================================

test_that("nuisance fits scale as O(patterns) + O(patterns x cost types)", {
  d <- make_cache_test_data(n_cost_types = 2)
  deltas <- list(rep(1.05, 6), rep(1.10, 6), rep(1.15, 6))

  cache <- new_nuisance_cache(555L, tags = d$cost_types)
  n_patterns <- NA_integer_

  for (ct in d$cost_types) {
    .clear_outcomes(cache)
    for (dl in deltas) {
      res <- run_cache_test_estimate(d, dl, ct, cache)
      if (is.na(n_patterns)) n_patterns <- res$n_patterns_used
    }
  }

  s <- nuisance_cache_stats(cache)
  n_cells <- length(d$cost_types) * length(deltas)

  # Propensity: fitted once per pattern for the whole sweep.
  expect_equal(s$n_ps_fit, n_patterns)
  expect_equal(s$n_ps_hit, n_patterns * (n_cells - 1))

  # Outcome: fitted once per pattern per cost type.
  expect_equal(s$n_outcome_fit, n_patterns * length(d$cost_types))
  expect_equal(s$n_outcome_hit,
               n_patterns * length(d$cost_types) * (length(deltas) - 1))

  # Total fits vs. the un-cached pipeline's P x CT x deltas.
  total_fits <- s$n_ps_fit + s$n_outcome_fit
  uncached_fits <- 2 * n_patterns * n_cells
  expect_lt(total_fits, uncached_fits)
  expect_equal(total_fits, n_patterns * (1 + length(d$cost_types)))

  # Miss on the first delta, hit on subsequent ones, for the same pattern within
  # the same cost type (>= 3 deltas exercised above).
  expect_equal(s$n_ps_miss, n_patterns)
  expect_equal(s$n_outcome_miss, n_patterns * length(d$cost_types))

  # Timing/diagnostics are populated.
  expect_true(s$ps_fit_seconds > 0)
  expect_true(s$outcome_fit_seconds > 0)
  expect_true(s$size_mb > 0)
  expect_output(print(cache), "nuisance_cache")
})

test_that("cached estimates agree with the un-cached pipeline within split noise", {
  d <- make_cache_test_data()
  delta <- rep(1.1, 6)

  cache <- new_nuisance_cache(31337L, tags = d$cost_types)
  res_cached <- run_cache_test_estimate(d, delta, "ct1", cache)

  set.seed(2468)
  res_plain <- run_cache_test_estimate(d, delta, "ct1", cache = NULL)

  # Different cross-fitting splits => not identical, but the same estimand.
  rel_diff <- abs(res_cached$estimate - res_plain$estimate) / abs(res_plain$estimate)
  expect_lt(rel_diff, 0.05)
  expect_lt(abs(res_cached$se - res_plain$se) / res_plain$se, 0.25)
})

# ============================================================================
# (iii) Fingerprint mismatch end-to-end (data changed under a warm cache)
# ============================================================================

test_that("reusing a warm cache on different rows aborts the run", {
  d <- make_cache_test_data()
  cache <- new_nuisance_cache(4242L, tags = d$cost_types)
  invisible(run_cache_test_estimate(d, rep(1.1, 6), "ct1", cache))

  # Same eligibility patterns, different rows (dropped observations): the cached
  # matrices no longer align, and that must abort rather than misalign silently.
  d_sub <- d
  keep <- setdiff(seq_len(nrow(d$base)), c(3, 17, 200))
  d_sub$base <- d$base[keep, , drop = FALSE]
  d_sub$outcomes <- lapply(d$outcomes, function(y) y[keep])
  names(d_sub$outcomes) <- names(d$outcomes)

  expect_error(run_cache_test_estimate(d_sub, rep(1.1, 6), "ct1", cache),
               class = "nuisance_cache_mismatch")

  # Corrupting a stored fingerprint by hand is caught the same way.
  cache2 <- new_nuisance_cache(4242L, tags = d$cost_types)
  invisible(run_cache_test_estimate(d, rep(1.1, 6), "ct1", cache2))
  key <- names(cache2$ps)[1]
  cache2$ps[[key]]$fingerprint$idx_sum <- cache2$ps[[key]]$fingerprint$idx_sum + 1
  expect_error(run_cache_test_estimate(d, rep(1.1, 6), "ct1", cache2),
               class = "nuisance_cache_mismatch")
})

# ============================================================================
# (iv) RNG neutrality
# ============================================================================

test_that("cached runs are RNG-neutral to the caller's stream", {
  d <- make_cache_test_data()
  delta <- rep(1.1, 6)

  # Reference: what the RNG stream yields with no estimation at all.
  set.seed(8080)
  reference_draws <- runif(5)

  # Fresh cache (all fits) -- fits happen inside withr::with_seed(), which
  # restores .Random.seed, so the ambient stream must be untouched.
  cache <- new_nuisance_cache(606L, tags = d$cost_types)
  set.seed(8080)
  seed_before <- .Random.seed
  invisible(run_cache_test_estimate(d, delta, "ct1", cache))
  expect_identical(.Random.seed, seed_before)
  expect_identical(runif(5), reference_draws)

  # Warm cache (all hits) -- same stream position, which is the property that
  # makes a resumed sweep's delta-scenario-level RNG independent of cache state.
  set.seed(8080)
  seed_before2 <- .Random.seed
  invisible(run_cache_test_estimate(d, delta, "ct1", cache))
  expect_identical(.Random.seed, seed_before2)
  expect_identical(runif(5), reference_draws)
})

# ============================================================================
# (vi) In-place nuisance mutation cannot corrupt the cache (TMLE)
# ============================================================================

test_that("TMLE targeting and feasibility subsetting leave cached matrices intact", {
  d <- make_cache_test_data()
  cache <- new_nuisance_cache(9090L, tags = d$cost_types)

  res1 <- run_cache_test_estimate(d, rep(1.1, 6), "ct1", cache, estimator = "tmle")

  # Deep copies (value comparison, not identity) of every cached matrix.
  snap_ps <- lapply(cache$ps, function(e) e$value + 0)
  snap_out <- lapply(cache$outcome[["ct1"]], function(e) e$value + 0)
  expect_true(length(snap_ps) > 0 && length(snap_out) > 0)

  # A second TMLE run at a different delta re-uses the same cached objects; the
  # TMLE targeting step assigns into a local copy, so the cache must be unchanged.
  res2 <- run_cache_test_estimate(d, rep(1.2, 6), "ct1", cache, estimator = "tmle")

  for (k in names(snap_ps)) {
    expect_identical(cache$ps[[k]]$value, snap_ps[[k]])
  }
  for (k in names(snap_out)) {
    expect_identical(cache$outcome[["ct1"]][[k]]$value, snap_out[[k]])
  }

  # And the second run's numbers match a fresh-cache TMLE run at the same delta:
  # if run 1 had corrupted the cached inputs, these would differ.
  cache_fresh <- new_nuisance_cache(9090L, tags = d$cost_types)
  res2_fresh <- run_cache_test_estimate(d, rep(1.2, 6), "ct1", cache_fresh,
                                        estimator = "tmle")
  expect_identical(res2$estimate, res2_fresh$estimate)
  expect_identical(res2$eif_values, res2_fresh$eif_values)
  expect_false(isTRUE(all.equal(res1$estimate, res2$estimate, tolerance = 1e-12)))
})

# ============================================================================
# Copy-on-modify: in-place nuisance edits must not touch cached matrices
# ============================================================================

test_that("feasibility subsetting does not mutate pre-supplied nuisance matrices", {
  # The other in-place assignment on nuisance matrices is the feasibility filter
  # (ps_est$propensity_scores <- ps_est$propensity_scores[feasible_obs, ]).
  # It is unreachable through estimate_pattern_effect() (which leaves
  # prob_infeasible at 0, so infeasibility errors instead of filtering), so it is
  # exercised here directly against pre-supplied (i.e. cache-style) matrices.
  n <- 100
  m <- 2
  # Column order for m = 2 is (A1,A2) = (0,0), (1,0), (0,1), (1,1).
  ps <- matrix(0.25, nrow = n, ncol = 4)
  # Five rows with support only on {(1,0), (0,1)}: the two marginals must sum to
  # 1 under any Q with that support, so delta = 3 (targets 0.75 + 0.75) is
  # infeasible for exactly those rows.
  bad <- 1:5
  ps[bad, ] <- matrix(rep(c(0, 0.5, 0.5, 0), each = length(bad)),
                      nrow = length(bad))
  outcomes <- matrix(seq_len(n * 4) / (n * 4), nrow = n, ncol = 4)

  set.seed(5)
  dat <- data.frame(A1 = rbinom(n, 1, 0.5), A2 = rbinom(n, 1, 0.5),
                    X1 = rnorm(n), Y = rnorm(n))

  ps_snapshot <- ps + 0
  outcomes_snapshot <- outcomes + 0

  res <- suppressWarnings(estimate_qdelta_mean(
    data = dat, delta = rep(3, m), treatment_cols = c("A1", "A2"),
    covariate_cols = "X1", outcome_col = "Y", divergence = "kl",
    oracle_propensities = ps, oracle_outcomes = outcomes,
    estimator = "one_step", cross_fit = FALSE, prob_infeasible = 0.1
  ))

  # The filter really engaged ...
  expect_equal(res$n_infeasible_excluded, length(bad))
  expect_equal(res$n, n - length(bad))
  # ... and the caller's matrices (which a cache would own) are untouched.
  expect_identical(ps, ps_snapshot)
  expect_identical(outcomes, outcomes_snapshot)
  expect_equal(dim(ps), c(n, 4L))
})

# ============================================================================
# Backward compatibility
# ============================================================================

test_that("stratified estimation without a cache is unaffected", {
  d <- make_cache_test_data(n = 400, n_cost_types = 1)
  delta <- rep(1.1, 6)

  set.seed(1234)
  res1 <- run_cache_test_estimate(d, delta, cache = NULL)
  set.seed(1234)
  res2 <- run_cache_test_estimate(d, delta, cache = NULL)

  # Same ambient seed => same result, exactly as before this feature existed.
  expect_identical(res1$estimate, res2$estimate)
  expect_true(is.finite(res1$estimate))
  # Timing block reports NA when nuisance fitting is not separable (no cache).
  expect_true(is.na(res1$nuisance_timing$nuisance_fit_seconds))
  expect_false(res1$nuisance_timing$cache_used)
})
