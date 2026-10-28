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
  expect_error(generate_binary_vectors(25), "not supported")
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
