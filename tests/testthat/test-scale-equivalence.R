# Equivalence tests for the scaling fixes (perf/scale-in-m): the optimised helpers must reproduce
# the previous implementations, which are kept below as slow references.

reference_binary_vectors <- function(m) {
  vectors <- matrix(0, nrow = 2^m, ncol = m)
  for (i in 0:(2^m - 1)) vectors[i + 1, ] <- as.integer(intToBits(i))[1:m]
  vectors
}

reference_marginal_polytope <- function(support_vector, m) {
  vectors <- reference_binary_vectors(m)
  vertices <- list()
  for (i in which(support_vector)) {
    marginals <- vectors[i, ]
    dup <- FALSE
    for (v in vertices) if (all(abs(v - marginals) < 1e-10)) { dup <- TRUE; break }
    if (!dup) vertices[[length(vertices) + 1]] <- marginals
  }
  do.call(rbind, vertices)
}

test_that("generate_binary_vectors matches the loop implementation exactly (m = 1..10)", {
  for (m in 1:10) {
    expect_identical(generate_binary_vectors(m), reference_binary_vectors(m))
  }
})

test_that("generate_binary_vectors is memoised and callers cannot corrupt the cache", {
  a <- generate_binary_vectors(4)
  a[1, 1] <- 99
  b <- generate_binary_vectors(4)
  expect_identical(b, reference_binary_vectors(4))
  expect_false(identical(a, b))
})

test_that("generate_binary_vectors rejects invalid and unsupported m", {
  expect_error(generate_binary_vectors(0), "positive integer")
  expect_error(generate_binary_vectors(2.5), "positive integer")
  expect_error(generate_binary_vectors(21), "not supported")
  expect_error(generate_binary_vectors("3"), "positive integer")
})

test_that("compute_marginal_polytope matches the de-duplicating reference", {
  set.seed(20261027)
  for (m in 2:6) {
    full <- rep(TRUE, 2^m)
    expect_identical(compute_marginal_polytope(full, m), reference_marginal_polytope(full, m))
    for (rep in 1:5) {
      support <- runif(2^m) < 0.6
      if (!any(support)) support[1] <- TRUE
      expect_identical(compute_marginal_polytope(support, m),
                       reference_marginal_polytope(support, m))
    }
  }
})

test_that("full-support shortcut agrees with the LP on feasibility and keeps the vertices", {
  set.seed(20261027)
  for (m in 2:5) {
    cube <- reference_marginal_polytope(rep(TRUE, 2^m), m)
    for (rep in 1:5) {
      p0 <- runif(m, 0.05, 0.95)
      delta <- runif(m, 0.3, 3)
      res <- check_delta_feasibility(delta, p0, rep(TRUE, 2^m))
      target <- delta * p0 / (1 + (delta - 1) * p0)
      lp <- check_point_in_marginal_polytope(target, cube)
      expect_equal(res$feasible, lp$feasible)
      expect_identical(res$polytope_vertices, cube)
      expect_equal(res$target_marginals, target)
      expect_true(res$distance_to_boundary > 0)
    }
  }
})

test_that("analyze_delta_feasibility_all_obs gives the same counts as per-observation LP checks", {
  set.seed(20261027)
  m <- 3
  n <- 40
  ps <- matrix(runif(n * 2^m), n, 2^m)
  ps[sample(length(ps), 60)] <- 0                      # partial support in some rows
  ps[1, ] <- 1                                         # make sure a row has full support
  ps <- ps / rowSums(ps)
  delta <- c(1.5, 2, 0.7)
  vectors <- reference_binary_vectors(m)

  expected <- vapply(seq_len(n), function(i) {
    p0 <- compute_marginals(ps[i, ], vectors)
    support <- ps[i, ] > 1e-8
    target <- delta * p0 / (1 + (delta - 1) * p0)
    if (any(target < 0 | target > 1)) return(FALSE)
    if (any(target == 0 | target == 1)) {
      for (k in which(target == 0 | target == 1)) {
        bad <- if (target[k] == 0) support & vectors[, k] == 1 else support & vectors[, k] == 0
        if (any(bad)) return(FALSE)
      }
    }
    check_point_in_marginal_polytope(target, reference_marginal_polytope(support, m))$feasible
  }, logical(1))

  res <- analyze_delta_feasibility_all_obs(ps, delta, early_stop = FALSE)
  expect_equal(res$feasible_per_obs, expected)
  expect_equal(res$n_feasible, sum(expected))
})

# ---- sequential propensity assembly -----------------------------------------------------------

reference_assemble <- function(pred_cache, m) {
  n <- nrow(pred_cache[[1]])
  vectors <- reference_binary_vectors(m)
  ps <- matrix(1, n, 2^m)
  for (j in seq_len(2^m)) {
    a <- vectors[j, ]
    for (k in seq_len(m)) {
      if (k == 1) {
        p_k <- pred_cache[[1]][, 1]
      } else {
        cond_vectors <- reference_binary_vectors(k - 1)
        idx <- which(apply(cond_vectors, 1, function(v) all(v == a[1:(k - 1)])))
        p_k <- pred_cache[[k]][, idx]
      }
      ps[, j] <- ps[, j] * if (a[k] == 1) p_k else 1 - p_k
    }
  }
  ps
}

random_pred_cache <- function(n, m) {
  lapply(seq_len(m), function(k) matrix(runif(n * 2^(k - 1), 0.05, 0.95), n, 2^(k - 1)))
}

test_that("assemble_sequential_joint reproduces the cell-by-cell chain-rule loop (m = 1..7)", {
  set.seed(20261027)
  for (m in 1:7) {
    pc <- random_pred_cache(25, m)
    expect_equal(assemble_sequential_joint(pc, m), reference_assemble(pc, m), tolerance = 1e-14)
  }
})

test_that("assemble_sequential_joint rows sum to one", {
  set.seed(20261027)
  ps <- assemble_sequential_joint(random_pred_cache(30, 6), 6)
  expect_equal(rowSums(ps), rep(1, 30), tolerance = 1e-12)
})

test_that("a non-default treatment_order is mapped back to treatment_cols order", {
  # Conditionals that ignore history => the joint is a product of marginals in ANY order, so the
  # answer is known exactly. The previous code evaluated the order-space joint at the wrong cells.
  set.seed(20261027)
  m <- 4; n <- 10
  marg <- matrix(runif(n * m, 0.1, 0.9), n, m)        # marginal P(A_c = 1), columns = treatment_cols
  order <- c(3, 1, 4, 2)                              # fitted order: treatment_cols[order]
  pc <- lapply(seq_len(m), function(k) {
    matrix(marg[, order[k]], n, 2^(k - 1))
  })
  ps <- assemble_sequential_joint(pc, m, perm = order)
  vectors <- reference_binary_vectors(m)
  truth <- sapply(seq_len(2^m), function(j) {
    apply(marg, 1, function(p) prod(ifelse(vectors[j, ] == 1, p, 1 - p)))
  })
  expect_equal(ps, truth, tolerance = 1e-12)
})

test_that("compute_joint_probabilities (glm conditionals) matches the reference assembly", {
  set.seed(20261027)
  n <- 300; m <- 4
  X <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
  A <- sapply(seq_len(m), function(k) rbinom(n, 1, plogis(0.3 * X[, 1] - 0.2 * k)))
  colnames(A) <- paste0("A", seq_len(m))
  d <- data.frame(X, A)
  models <- lapply(seq_len(m), function(k) {
    preds <- c("x1", "x2", if (k > 1) paste0("A", seq_len(k - 1)))
    glm(reformulate(preds, paste0("A", k)), data = d, family = binomial())
  })
  new <- compute_joint_probabilities(d, models, paste0("A", seq_len(m)),
                                     generate_binary_vectors(m), c("x1", "x2"))
  pc <- vector("list", m)
  pc[[1]] <- matrix(predict(models[[1]], d, type = "response"), n, 1)
  for (k in 2:m) {
    cv <- reference_binary_vectors(k - 1)
    pc[[k]] <- sapply(seq_len(2^(k - 1)), function(ci) {
      nd <- d
      for (jp in seq_len(k - 1)) nd[[paste0("A", jp)]] <- cv[ci, jp]
      predict(models[[k]], nd, type = "response")
    })
  }
  expect_equal(new, reference_assemble(pc, m), tolerance = 1e-12)
})

# ---- EIF chunk sizing --------------------------------------------------------------------------

test_that("eif_default_chunk_size respects the per-worker memory budget at large m", {
  for (m in c(6, 11, 12, 14, 16)) {
    for (K in c(1, 9)) {
      cs <- eif_default_chunk_size(m, K)
      expect_gte(cs, EIF_MIN_CHUNK)
      expect_lte(cs, EIF_MAX_CHUNK)
      if (cs > EIF_MIN_CHUNK) {
        expect_lte(cs * EIF_LIVE_MATRICES * 2^m * 8 * K, EIF_CHUNK_BUDGET_BYTES * 1.0001)
      }
    }
  }
  expect_equal(eif_default_chunk_size(6, 1), EIF_MAX_CHUNK)   # small m: the upper bound binds
  expect_lt(eif_default_chunk_size(14, 9), eif_default_chunk_size(11, 9))
})

# ---- Q_delta solver failure policy and tolerance -----------------------------------------------

test_that("resolve_qdelta_tol defaults to QDELTA_TOLERANCE and validates overrides", {
  expect_equal(resolve_qdelta_tol(NULL), QDELTA_TOLERANCE)
  expect_equal(resolve_qdelta_tol(1e-6), 1e-6)
  expect_error(resolve_qdelta_tol(-1), "single positive")
  expect_error(resolve_qdelta_tol(c(1e-6, 1e-8)), "single positive")
  expect_error(resolve_qdelta_tol("a"), "single positive")
})

test_that("resolve_max_failure_rate accepts [0, 1] and rejects everything else", {
  expect_equal(resolve_max_failure_rate(0), 0)
  expect_equal(resolve_max_failure_rate(0.02), 0.02)
  expect_error(resolve_max_failure_rate(-0.1), "single number in")
  expect_error(resolve_max_failure_rate(1.5), "single number in")
  expect_error(resolve_max_failure_rate(c(0, 0.1)), "single number in")
  expect_error(resolve_max_failure_rate(NA_real_), "single number in")
})

test_that("qdelta_fallback_or_stop stops at the first failure by default", {
  ps <- c(0.4, 0.3, 0.2, 0.1)
  e <- simpleError("did not converge")
  expect_error(qdelta_fallback_or_stop(7, e, ps, c(1.5, 2), 1L, 100L, 0),
               "observation 7.*did not converge")
  expect_error(qdelta_fallback_or_stop(7, e, ps, c(1.5, 2), 1L, 100L, 0), "max_qdelta_failure_rate")
  cond <- tryCatch(qdelta_fallback_or_stop(7, e, ps, c(1.5, 2), 1L, 100L, 0), error = identity)
  expect_s3_class(cond, "qdelta_solver_failure")
})

test_that("qdelta_fallback_or_stop tolerates failures up to rate * n, then stops", {
  ps <- c(0.4, 0.3, 0.2, 0.1)
  e <- simpleError("did not converge")
  # 2% of 100 = 2 failures allowed: the 1st and 2nd fall back, the 3rd stops
  expect_warning(r1 <- qdelta_fallback_or_stop(7, e, ps, c(1.5, 2), 1L, 100L, 0.02), "fallback")
  expect_equal(r1$q_delta, ps)
  expect_equal(r1$lambda, rep(0, 3))
  expect_false(r1$converged)
  expect_warning(qdelta_fallback_or_stop(8, e, ps, c(1.5, 2), 2L, 100L, 0.02), "fallback")
  expect_error(qdelta_fallback_or_stop(9, e, ps, c(1.5, 2), 3L, 100L, 0.02),
               "exceeded max_qdelta_failure_rate: 3 of 100")
  # rate 1 always falls back
  expect_warning(qdelta_fallback_or_stop(9, e, ps, c(1.5, 2), 100L, 100L, 1), "fallback")
})

test_that("estimate_qdelta_mean stops on a solver failure unless max_qdelta_failure_rate allows it", {
  set.seed(20261027)
  n <- 120
  d <- data.frame(A1 = rbinom(n, 1, 0.5), A2 = rbinom(n, 1, 0.4), X1 = rnorm(n), Y = rnorm(n))
  original <- compute_qdelta_per_obs
  calls <- 0L
  fail_obs_3 <- function(ps, delta, divergence, ...) {
    calls <<- calls + 1L
    if (calls == 3L) stop("synthetic solver failure")
    original(ps, delta, divergence, ...)
  }
  call_est <- function(...) {
    estimate_qdelta_mean(d, delta = c(1.5, 1.5), treatment_cols = c("A1", "A2"),
                         covariate_cols = "X1", outcome_col = "Y", cross_fit = FALSE, ...)
  }
  testthat::local_mocked_bindings(compute_qdelta_per_obs = fail_obs_3)
  expect_error(suppressMessages(call_est()), "Q_delta solver failed for observation 3")
  calls <- 0L
  expect_warning(
    res <- suppressMessages(call_est(max_qdelta_failure_rate = 0.05)),   # 1 of 120 = 0.8%
    "fallback to P0 used for 1 of"
  )
  expect_true(is.finite(res$estimate))
  calls <- 0L
  expect_error(suppressWarnings(suppressMessages(call_est(max_qdelta_failure_rate = 0.001))),
               "exceeded max_qdelta_failure_rate")
})

test_that("a Q_delta solver failure aborts the stratified estimate instead of dropping the pattern", {
  set.seed(20261028)
  n <- 400
  d <- data.frame(E1 = 1, E2 = rbinom(n, 1, 0.7), A1 = rbinom(n, 1, 0.5), A2 = rbinom(n, 1, 0.4),
                  X1 = rnorm(n), Y = rnorm(n))
  d$A2[d$E2 == 0] <- 0
  original <- compute_qdelta_per_obs
  calls <- 0L
  failing <- function(ps, delta, divergence, ...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("synthetic solver failure")
    original(ps, delta, divergence, ...)
  }
  testthat::local_mocked_bindings(compute_qdelta_per_obs = failing)
  run <- function(...) {
    estimate_qdelta_mean_stratified(
      d, delta = c(1.5, 1.5), covariate_cols = "X1", eligibility_cols = c("E1", "E2"),
      outcome_col = "Y", min_stratum_size = 50, cross_fit = FALSE, verbose = FALSE, ...)
  }
  skip_if(!is.function(estimate_qdelta_mean_stratified))
  expect_error(suppressWarnings(suppressMessages(run())), "Q_delta solver failed for observation")
})

test_that("pattern checkpoints are saved, reloaded on an identical rerun, and invalidated by a changed delta", {
  set.seed(20261029)
  n <- 500
  d <- data.frame(elig_A1 = 1L, elig_A2 = rbinom(n, 1, 0.6), X1 = rnorm(n), X2 = rnorm(n), Y = rnorm(n))
  d$A1 <- rbinom(n, 1, plogis(0.3 * d$X1))
  d$A2 <- ifelse(d$elig_A2 == 1, rbinom(n, 1, plogis(0.2 * d$X2)), 0L)
  ck <- withr::local_tempdir()
  run <- function(delta) {
    estimate_qdelta_mean_stratified(
      d, delta = delta, covariate_cols = c("X1", "X2"), eligibility_cols = c("elig_A1", "elig_A2"),
      outcome_col = "Y", min_stratum_size = 50, cross_fit = FALSE, verbose = FALSE,
      pattern_checkpoint_dir = ck, pattern_checkpoint_tag = "total__d1.5")
  }
  expect_message(first <- suppressWarnings(run(c(1.5, 1.5))), "saved")
  files <- list.files(ck, pattern = "\\.rds$")
  expect_gte(length(files), 2L)
  expect_false(any(grepl("\\.tmp$", list.files(ck))))
  expect_message(second <- suppressWarnings(run(c(1.5, 1.5))), "reloaded from checkpoint")
  expect_equal(second$estimate, first$estimate)
  expect_equal(second$se, first$se)
  expect_message(third <- suppressWarnings(run(c(1.2, 1.2))), "stale or unusable")
  expect_false(isTRUE(all.equal(third$estimate, first$estimate)))
})

test_that("a permuted order with history-DEPENDENT conditionals matches a brute-force joint", {
  # The previous test uses history-ignoring conditionals, which cannot expose a mismapping of the
  # history bits. Here conditional k depends on the earlier treatments in fitted order, bit p-1 of the
  # history index = value of fitted treatment p; the truth is built cell by cell in treatment_cols order.
  set.seed(20261030)
  m <- 4; n <- 6
  order <- c(3, 1, 4, 2)
  pc <- lapply(seq_len(m), function(k) matrix(runif(n * 2^(k - 1), 0.1, 0.9), n, 2^(k - 1)))
  ps <- assemble_sequential_joint(pc, m, perm = order)
  vectors <- reference_binary_vectors(m)            # rows in treatment_cols order, column 1 fastest
  truth <- matrix(1, n, 2^m)
  for (j in seq_len(2^m)) {
    a <- vectors[j, ]
    for (k in seq_len(m)) {
      hist_idx <- if (k == 1) 0 else sum(a[order[seq_len(k - 1)]] * 2^(seq_len(k - 1) - 1))
      p1 <- pc[[k]][, hist_idx + 1]
      truth[, j] <- truth[, j] * if (a[order[k]] == 1) p1 else 1 - p1
    }
  }
  expect_equal(ps, truth, tolerance = 1e-12)
  expect_equal(rowSums(ps), rep(1, n), tolerance = 1e-12)
})

test_that("failure-rate edge cases: float boundary, reverse_kl, and the count is returned", {
  ps <- c(0.4, 0.3, 0.2, 0.1)
  e <- simpleError("did not converge")
  # 0.29 * 100 is 28.999999999999996 in floating point: exactly 29 failures must still be tolerated
  expect_warning(qdelta_fallback_or_stop(1, e, ps, c(1.5, 2), 29L, 100L, 0.29), "fallback")
  expect_error(qdelta_fallback_or_stop(1, e, ps, c(1.5, 2), 30L, 100L, 0.29), "exceeded")
  expect_error(suppressWarnings(qdelta_fallback_or_stop(1, e, ps, c(1.5, 2), 1L, 100L, 0.5, "reverse_kl")),
               "reverse_kl")

  set.seed(20261031)
  n <- 120
  d <- data.frame(A1 = rbinom(n, 1, 0.5), A2 = rbinom(n, 1, 0.4), X1 = rnorm(n), Y = rnorm(n))
  original <- compute_qdelta_per_obs
  calls <- 0L
  fail_obs_3 <- function(ps, delta, divergence, ...) {
    calls <<- calls + 1L
    if (calls == 3L) stop("synthetic solver failure")
    original(ps, delta, divergence, ...)
  }
  testthat::local_mocked_bindings(compute_qdelta_per_obs = fail_obs_3)
  expect_warning(
    res <- suppressMessages(estimate_qdelta_mean(d, delta = c(1.5, 1.5), treatment_cols = c("A1", "A2"),
                                                 covariate_cols = "X1", outcome_col = "Y", cross_fit = FALSE,
                                                 max_qdelta_failure_rate = 0.05)),
    "fallback to P0 used for 1 of")
  expect_identical(res$n_qdelta_fallback, 1L)
})

test_that("the vanishing-threshold estimator honours max_qdelta_failure_rate", {
  set.seed(20261101)
  n <- 150
  d <- data.frame(A1 = rbinom(n, 1, 0.5), A2 = rbinom(n, 1, 0.4), X1 = rnorm(n), Y = rnorm(n))
  original <- compute_qdelta_per_obs
  calls <- 0L
  failing <- function(ps, delta, divergence, ...) {
    calls <<- calls + 1L
    if (calls == 5L) stop("synthetic solver failure")
    original(ps, delta, divergence, ...)
  }
  testthat::local_mocked_bindings(compute_qdelta_per_obs = failing)
  run <- function(...) {
    estimate_qdelta_vanishing_threshold(d, delta = c(1.5, 1.5), treatment_cols = c("A1", "A2"),
                                        covariate_cols = "X1", outcome_col = "Y", cross_fit = FALSE, ...)
  }
  expect_error(suppressWarnings(suppressMessages(run())), "Q_delta solver failed for observation")
  calls <- 0L
  expect_warning(res <- suppressMessages(run(max_qdelta_failure_rate = 0.05)), "fallback to P0")
  expect_true(is.finite(res$estimate))
})

test_that("treatment_order with a repeated column is rejected", {
  expect_error(suppressWarnings(
    estimate_propensity_scores(data.frame(A1 = c(0, 1, 0, 1), A2 = c(1, 1, 0, 0), X1 = 1:4, Y = 1:4),
                               c("A1", "A2"), "X1", method = "sequential", binary_method = "glm",
                               cross_fit = FALSE, treatment_order = c("A1", "A1"))),
    "exactly once")
})

test_that("a pattern checkpoint is not reloaded when the outcome or covariates change but the row counts do not", {
  set.seed(20261102)
  n <- 500
  d <- data.frame(elig_A1 = 1L, elig_A2 = rbinom(n, 1, 0.6), X1 = rnorm(n), X2 = rnorm(n), Y = rnorm(n))
  d$A1 <- rbinom(n, 1, plogis(0.3 * d$X1))
  d$A2 <- ifelse(d$elig_A2 == 1, rbinom(n, 1, plogis(0.2 * d$X2)), 0L)
  ck <- withr::local_tempdir()
  run <- function(dat) {
    estimate_qdelta_mean_stratified(
      dat, delta = c(1.5, 1.5), covariate_cols = c("X1", "X2"), eligibility_cols = c("elig_A1", "elig_A2"),
      outcome_col = "Y", min_stratum_size = 50, cross_fit = FALSE, verbose = FALSE,
      pattern_checkpoint_dir = ck, pattern_checkpoint_tag = "total__d1.5")
  }
  suppressMessages(suppressWarnings(run(d)))
  d2 <- d; d2$Y <- d2$Y + 1                      # same n, same pattern sizes, different outcome
  expect_message(suppressWarnings(run(d2)), "stale or unusable")
  d3 <- d; d3$X2 <- rev(d3$X2)                   # a covariate changed
  expect_message(suppressWarnings(run(d3)), "stale or unusable")
  expect_message(suppressWarnings(run(d)), "stale or unusable")   # d2 and d3 overwrote the files; d is recomputed
})

test_that("stratified_contrast_vs_observed: observed-mean contrast needs no delta = 1 run and checks alignment", {
  set.seed(20261103)
  n <- 900
  d <- data.frame(elig_A1 = 1L, elig_A2 = rbinom(n, 1, 0.6), X1 = rnorm(n), X2 = rnorm(n))
  d$A1 <- rbinom(n, 1, plogis(0.3 * d$X1))
  d$A2 <- ifelse(d$elig_A2 == 1, rbinom(n, 1, plogis(0.2 * d$X2)), 0L)
  d$Y <- 1 + 0.8 * d$A1 + 0.5 * d$A2 + d$X1 + rnorm(n)
  fit <- function(delta) suppressWarnings(suppressMessages(estimate_qdelta_mean_stratified(
    d, delta = delta, covariate_cols = c("X1", "X2"), eligibility_cols = c("elig_A1", "elig_A2"),
    outcome_col = "Y", min_stratum_size = 50, cross_fit = FALSE, verbose = FALSE)))
  f2 <- fit(c(2, 2)); f1 <- fit(c(1, 1))
  ct <- stratified_contrast_vs_observed(f2$pattern_contributions, d$Y)
  # level identities: the observed side is the mean of Y over the included patients, the other side is the fit
  expect_equal(ct$estimate_a, mean(d$Y), tolerance = 1e-10)
  expect_equal(ct$estimate_b, f2$estimate, tolerance = 1e-10)
  expect_equal(ct$estimate, f2$estimate - mean(d$Y), tolerance = 1e-10)
  # a real delta = 1 run gives (up to nuisance noise) the same contrast
  ct1 <- stratified_contrast(f1$pattern_contributions, f2$pattern_contributions)
  expect_equal(ct$estimate, ct1$estimate, tolerance = 0.1)
  # the no-shift contrast is about zero
  zero <- stratified_contrast_vs_observed(f1$pattern_contributions, d$Y)
  expect_lt(abs(zero$estimate), 1e-6)    # at delta = 1 the one-step estimate IS the observed mean
  # the contrast SE is smaller than the SE of the level (it removes the outcome's own variance)
  expect_lt(ct$se, f2$se)
  # a misaligned or incomplete outcome vector is refused
  expect_error(stratified_contrast_vs_observed(f2$pattern_contributions, rev(d$Y)[-1]), "one value per row")
  # (a permutation of y with the right pattern sizes cannot be detected here; the applied script checks
  # row alignment against the eligibility patterns before calling this)
})

test_that("skip_null_shift: unshifted patterns contribute the observed mean exactly and cost no estimation", {
  set.seed(20261106)
  n <- 800
  d <- data.frame(elig_A1 = 1L, elig_A2 = rbinom(n, 1, 0.5), X1 = rnorm(n))
  d$A1 <- rbinom(n, 1, plogis(0.3 * d$X1)); d$A2 <- ifelse(d$elig_A2 == 1, rbinom(n, 1, 0.4), 0L)
  d$Y <- 1 + d$A1 + d$A2 + d$X1 + rnorm(n)
  run <- function(delta, skip) suppressWarnings(suppressMessages(estimate_qdelta_mean_stratified(
    d, delta = delta, covariate_cols = "X1", eligibility_cols = c("elig_A1", "elig_A2"), outcome_col = "Y",
    min_stratum_size = 50, cross_fit = FALSE, verbose = FALSE, skip_null_shift = skip)))
  # shift only A2: the pattern with only A1 eligible is a null shift
  a <- run(c(1, 2), TRUE); b <- run(c(1, 2), FALSE)
  pc <- a$pattern_contributions$patterns
  only_a1 <- which(pc$n_eligible == 1)
  expect_length(only_a1, 1L)
  ids <- a$pattern_contributions$obs_pattern_id
  y_pat <- d$Y[!is.na(ids) & ids == pc$pattern_id[only_a1]]
  expect_equal(pc$estimate[only_a1], mean(y_pat), tolerance = 1e-12)
  # the shifted pattern is unchanged, and the overall estimate agrees with the full machinery up to nuisance noise
  expect_equal(a$pattern_contributions$patterns$estimate[-only_a1], b$pattern_contributions$patterns$estimate[-only_a1],
               tolerance = 1e-10)
  expect_equal(a$estimate, b$estimate, tolerance = 0.1)
  # an all-ones delta: every pattern is skipped and the contrast against the observed mean is exactly 0
  z <- run(c(1, 1), TRUE)
  zc <- stratified_contrast_vs_observed(z$pattern_contributions, d$Y)
  expect_equal(zc$estimate, 0, tolerance = 1e-10)
  expect_equal(zc$se, 0, tolerance = 1e-8)
  # the checkpoint signature separates skipped from full runs
  ck <- withr::local_tempdir()
  go <- function(skip) suppressWarnings(estimate_qdelta_mean_stratified(
    d, delta = c(1, 2), covariate_cols = "X1", eligibility_cols = c("elig_A1", "elig_A2"), outcome_col = "Y",
    min_stratum_size = 50, cross_fit = FALSE, verbose = FALSE, skip_null_shift = skip,
    pattern_checkpoint_dir = ck, pattern_checkpoint_tag = "t"))
  suppressMessages(go(TRUE))
  expect_message(go(FALSE), "stale or unusable")
})

test_that("qdelta_memo: a second outcome reuses the Q_delta solves and gives identical results", {
  set.seed(20261107)
  n <- 700
  d <- data.frame(elig_A1 = 1L, elig_A2 = rbinom(n, 1, 0.5), X1 = rnorm(n), X2 = rnorm(n))
  d$A1 <- rbinom(n, 1, plogis(0.3 * d$X1)); d$A2 <- ifelse(d$elig_A2 == 1, rbinom(n, 1, 0.4), 0L)
  d$Y1 <- 1 + d$A1 + d$X1 + rnorm(n); d$Y2 <- 5 - d$A2 + d$X2 + rnorm(n)
  calls <- 0L
  original <- compute_qdelta_per_obs
  counting <- function(...) { calls <<- calls + 1L; original(...) }
  testthat::local_mocked_bindings(compute_qdelta_per_obs = counting)
  cache <- new_nuisance_cache(20260706L, tags = c("y1", "y2"))
  run <- function(y, tag, memo) {
    dd <- d; dd$Y <- d[[y]]
    suppressWarnings(suppressMessages(estimate_qdelta_mean_stratified(
      dd, delta = c(1.5, 1.5), covariate_cols = c("X1", "X2"), eligibility_cols = c("elig_A1", "elig_A2"),
      outcome_col = "Y", min_stratum_size = 50, cross_fit = TRUE, n_folds = 3, verbose = FALSE,
      nuisance_cache = cache, nuisance_tag = tag, qdelta_memo = memo)))
  }
  memo <- new.env()
  a1 <- run("Y1", "y1", memo); n1 <- calls
  a2 <- run("Y2", "y2", memo); n2 <- calls - n1
  expect_gt(n1, 0L)
  expect_identical(n2, 0L)                         # nothing re-solved for the second outcome
  # identical to a run that solves everything itself (no memo, fresh cache with the same seeds)
  cache <- new_nuisance_cache(20260706L, tags = c("y1", "y2"))
  calls <- 0L
  b1 <- run("Y1", "y1", NULL); b2 <- run("Y2", "y2", NULL)
  expect_identical(a1$estimate, b1$estimate); expect_identical(a1$se, b1$se)
  expect_identical(a2$estimate, b2$estimate); expect_identical(a2$se, b2$se)
  expect_identical(a2$eif_values, b2$eif_values)
  # a different delta does not hit the memo
  before <- calls
  dd <- d; dd$Y <- d$Y1
  suppressWarnings(suppressMessages(estimate_qdelta_mean_stratified(
    dd, delta = c(2, 2), covariate_cols = c("X1", "X2"), eligibility_cols = c("elig_A1", "elig_A2"),
    outcome_col = "Y", min_stratum_size = 50, cross_fit = TRUE, n_folds = 3, verbose = FALSE,
    nuisance_cache = cache, nuisance_tag = "y1", qdelta_memo = memo)))
  expect_gt(calls, before)
})

test_that("the pattern-major grid reproduces the per-run stratified estimator exactly, and resumes from disk", {
  set.seed(20261108)
  n <- 900
  d <- data.frame(elig_A1 = 1L, elig_A2 = rbinom(n, 1, 0.6), elig_A3 = rbinom(n, 1, 0.5),
                  X1 = rnorm(n), X2 = rnorm(n))
  d$A1 <- rbinom(n, 1, plogis(0.3 * d$X1))
  d$A2 <- ifelse(d$elig_A2 == 1, rbinom(n, 1, plogis(0.2 * d$X2)), 0L)
  d$A3 <- ifelse(d$elig_A3 == 1, rbinom(n, 1, 0.4), 0L)
  d$Y_a <- 1 + d$A1 + d$X1 + rnorm(n); d$Y_b <- 5 - d$A2 + d$X2 + 0.5 * d$A3 + rnorm(n)
  oc <- c(a = "Y_a", b = "Y_b")
  deltas <- list(up = c(2, 2, 2), dom = c(1, 2, 1))
  common <- list(covariate_cols = c("X1", "X2"), eligibility_cols = c("elig_A1", "elig_A2", "elig_A3"),
                 min_stratum_size = 50, cross_fit = TRUE, n_folds = 3, ps_method = "sequential",
                 binary_method = "glm", outcome_method = "glm")
  ck <- withr::local_tempdir()
  g <- suppressWarnings(suppressMessages(do.call(estimate_qdelta_stratified_grid, c(
    list(data = d, outcome_cols = oc, deltas = deltas, seed_base = 20260706L, checkpoint_dir = ck, verbose = FALSE,
         skip_null_shift = TRUE), common))))
  expect_length(g$runs, 4L)
  # the reference: the existing estimator, one run at a time, with a shared cache and the same seeds
  cache <- new_nuisance_cache(20260706L, tags = names(oc))
  for (nm in names(oc)) for (sc in names(deltas)) {
    dd <- d; dd$Y <- d[[oc[[nm]]]]
    ref <- suppressWarnings(suppressMessages(do.call(estimate_qdelta_mean_stratified, c(
      list(data = dd, delta = deltas[[sc]], outcome_col = "Y", verbose = FALSE, nuisance_cache = cache,
           nuisance_tag = nm, skip_null_shift = TRUE), common))))
    got <- g$runs[[paste(nm, sc, sep = "__")]]
    expect_true(got$success)
    expect_identical(got$estimate, ref$estimate)
    expect_identical(got$se, ref$se)
    expect_identical(got$eif_values, ref$eif_values)
    expect_identical(got$pattern_contributions$obs_eif, ref$pattern_contributions$obs_eif)
  }
  # resume: nothing is recomputed, results identical; a new scenario computes only itself
  files <- list.files(ck, pattern = "^pattern_.*rds$")
  expect_gte(length(files), 2L)
  calls <- 0L
  original <- compute_qdelta_per_obs
  testthat::local_mocked_bindings(compute_qdelta_per_obs = function(...) { calls <<- calls + 1L; original(...) })
  g2 <- suppressWarnings(suppressMessages(do.call(estimate_qdelta_stratified_grid, c(
    list(data = d, outcome_cols = oc, deltas = deltas, seed_base = 20260706L, checkpoint_dir = ck, verbose = FALSE,
         skip_null_shift = TRUE), common))))
  expect_identical(calls, 0L)
  expect_identical(g2$runs$a__up$estimate, g$runs$a__up$estimate)
  g3 <- suppressWarnings(suppressMessages(do.call(estimate_qdelta_stratified_grid, c(
    list(data = d, outcome_cols = oc, deltas = c(deltas, list(low = c(0.5, 0.5, 0.5))), seed_base = 20260706L,
         checkpoint_dir = ck, verbose = FALSE, skip_null_shift = TRUE), common))))
  expect_gt(calls, 0L)
  expect_identical(g3$runs$a__up$estimate, g$runs$a__up$estimate)
  expect_true(g3$runs$b__low$success)
  # a changed outcome invalidates only the units that use it
  d2 <- d; d2$Y_a <- d2$Y_a + 1
  calls <- 0L
  g4 <- suppressWarnings(suppressMessages(do.call(estimate_qdelta_stratified_grid, c(
    list(data = d2, outcome_cols = oc, deltas = deltas, seed_base = 20260706L, checkpoint_dir = ck, verbose = FALSE,
         skip_null_shift = TRUE), common))))
  expect_gt(calls, 0L)
  expect_identical(g4$runs$b__up$estimate, g$runs$b__up$estimate)
  expect_false(identical(g4$runs$a__up$estimate, g$runs$a__up$estimate))
})
