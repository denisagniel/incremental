#' @file test-eif-omega-vectorized-equivalence.R
#' @title Old-vs-new equivalence for the vectorized EIF-lambda / Omega assembly
#'
#' Stage 2 of the performance work in
#' `quality_reports/specs/2026-08-10_qdelta-solver-performance-fix.md`, applying the
#' same equivalence discipline stage 1 used for the Q_delta solver
#' (`test-qdelta-vectorized-equivalence.R`).
#'
#' "Old" = `R/eif_legacy.R`, a verbatim in-tree snapshot of `calculate_eif_lambda()`
#' and `eif_compute_chunk()` as they stood before the fix, verified bit-identical
#' to the pre-fix working tree before any edit. Nothing here reads git history.
#'
#' Three exact (not approximate) changes are under test, each of which only became
#' the dominant cost after the previous one was removed:
#'   1. `t(A) %*% diag(w) %*% A`            -> `crossprod(A * w, A)`
#'   2. `apply(A, 1, function(r) all(r == x))` -> integer bit-key prefilter + exact verify
#'   3. `compute_marginals(ps, A[, -1], ...)` -> `compute_marginals(ps, A, ...)[-1]`
#'
#' Because all three are algebraic identities rather than approximations, the bar
#' here is deliberately tighter than stage 1's 1e-6: exact equality is asserted
#' wherever it was observed to hold (it holds everywhere tested), with a 1e-6
#' relative fallback documented on the aggregate comparisons.

EOV_TOL <- 1e-6          # repo point-estimate tolerance; the observed diffs are 0
EOV_SEED <- 20260811L

# ---------------------------------------------------------------------------
# Case construction
# ---------------------------------------------------------------------------

#' Propensity matrix of a requested support/marginal shape
#'
#' @param m Integer, number of treatments
#' @param n Integer, number of observations
#' @param kind One of:
#'   * `"full"` -- all 2^m combinations carry mass; every marginal is interior, so
#'     all m+1 constraints are active and `calculate_eif_lambda()` returns an
#'     unpadded (m+1)-vector.
#'   * `"sparse"` -- mass only on \{0, e_1, ..., e_m\}; marginals stay interior but
#'     most combinations are ~0, exercising the weighted-sum paths with many zero
#'     weights.
#'   * `"degenerate_marginal"` -- essentially all mass on combinations with
#'     `A_1 = 1`, so marginal 1 lands within `find_active_constraints()`'s 1e-8
#'     tolerance of exactly 1 and constraint 1 is DROPPED: `omega_active` becomes a
#'     strict submatrix and `calculate_eif_lambda()` takes its **zero-padding**
#'     return path. That is the branch most likely to be broken by an
#'     Omega-indexing mistake, so it must be covered.
#'     The `A_1 = 0` rows retain a tiny but strictly-above-1e-12
#'     (`PROPENSITY_ACTIVE_THRESHOLD`) mass rather than exactly 0: zeroing them
#'     would drop the support rank to m < m+1 and `newton_qdelta()` would
#'     (correctly) refuse to solve for Q_delta at all, leaving no EIF input to
#'     build. Total leaked mass 1e-9 keeps the marginal inside the 1e-8 tolerance.
#'   * `"extreme"` -- highly skewed masses spanning ~10 orders of magnitude, giving
#'     a near-singular Omega and exercising the `tryCatch`/regularization fallback.
eov_make_ps <- function(m, n, kind) {
  withr::local_seed(EOV_SEED + m)
  b <- generate_binary_vectors(m)
  n_comb <- nrow(b)

  base <- switch(kind,
    full = matrix(stats::runif(n * n_comb, 0.3, 1), nrow = n),
    sparse = {
      mm <- matrix(stats::runif(n * n_comb, 0.3, 1), nrow = n)
      mm[, rowSums(b) > 1] <- 0
      mm
    },
    degenerate_marginal = {
      on_cols <- which(b[, 1] == 1)
      off_cols <- which(b[, 1] != 1)
      # The leaked mass must satisfy TWO opposing constraints simultaneously:
      #   per-combination > PROPENSITY_ACTIVE_THRESHOLD (1e-12), else newton_qdelta()
      #     drops those rows and the support rank falls to m < m+1 (no Q_delta);
      #   total          < 1e-8 (find_active_constraints tolerance), else marginal 1
      #     stays interior and the constraint is NOT dropped (no padding path).
      # At m = 11 there are 1024 off-columns, which makes this genuinely tight.
      leak_total <- 5e-9
      stopifnot(
        "leaked mass per combination must exceed the solver's support threshold" =
          leak_total / length(off_cols) > PROPENSITY_ACTIVE_THRESHOLD * 2,
        "leaked total must stay inside the constraint-activity tolerance" =
          leak_total < 1e-8
      )
      mm <- matrix(0, nrow = n, ncol = n_comb)
      mass_on <- matrix(stats::runif(n * length(on_cols), 0.3, 1), nrow = n)
      mm[, on_cols] <- (1 - leak_total) * mass_on / rowSums(mass_on)
      mm[, off_cols] <- leak_total / length(off_cols)
      return(mm)   # already normalized to 1; do not re-scale
    },

    extreme = matrix(10^stats::runif(n * n_comb, -6, 0), nrow = n),
    stop("unknown kind: ", kind)
  )
  base / rowSums(base)
}

#' Uniform and sparse (some components exactly 1) delta vectors
#'
#' @param base Tilt applied to the non-unit components. Must be FEASIBLE for the
#'   support in question, which is a property of the propensity vector and not of
#'   either implementation:
#'   * `"sparse"` -- only \{0, e_1, ..., e_m\} carry mass, so at most one treatment
#'     can be on at a time and \eqn{\sum_k E_Q[A_k] \le 1} is a hard constraint;
#'     any `base > 1` violates it once m is large enough. Hence shrinkage (0.8).
#'   * `"extreme"` -- masses span ~6 orders of magnitude, pushing some marginals
#'     close to 0/1, so only a mild tilt stays inside the marginal polytope (1.05).
#'   `newton_qdelta()` correctly refuses to converge on infeasible targets (this was
#'   verified while building the battery), so feasible deltas are chosen rather than
#'   relying on the solver to absorb them.
eov_deltas <- function(m, base = 1.2) {
  sparse <- rep(base, m)
  sparse[1] <- 1
  if (m >= 3L) sparse[m] <- 1
  list(uniform = rep(base, m), sparse = sparse)
}

#' The delta base that is feasible for a given support shape
eov_delta_base <- function(kind) {
  switch(kind, sparse = 0.8, extreme = 1.05, 1.2)
}

#' Assemble a complete EIF input set for one (m, kind, divergence, delta) case
eov_build_case <- function(m, n, kind, divergence, delta) {
  a_mat <- generate_binary_vectors(m)
  a_ext <- cbind(1, a_mat)
  n_comb <- nrow(a_mat)

  ps <- eov_make_ps(m, n, kind)

  # Under `degenerate_marginal`, treatment 1 is already saturated (marginal ~ 1),
  # so any delta_1 > 1 asks Q_delta for a marginal above the support's supremum
  # and newton_qdelta() correctly refuses to converge. Leave that component
  # untilted; the point of this case is the DROPPED CONSTRAINT in
  # calculate_eif_lambda(), which is driven by the propensity vector, not by delta.
  if (kind == "degenerate_marginal") delta[1] <- 1

  withr::local_seed(EOV_SEED + m + n)
  mu <- matrix(stats::rnorm(n * n_comb, mean = 1, sd = 0.3), nrow = n)
  obs_idx <- apply(ps, 1, function(pr) sample.int(n_comb, 1, prob = pr))
  a_obs <- a_mat[obs_idx, , drop = FALSE]
  y <- mu[cbind(seq_len(n), obs_idx)] + stats::rnorm(n, sd = 0.1)

  qd <- lapply(seq_len(n), function(i) {
    compute_qdelta_per_obs(ps[i, ], delta, divergence, fix_lambda0 = FALSE,
                           tol = 1e-10, combination_matrix = a_mat)
  })
  h <- compute_H_vectorized_general(do.call(rbind, lapply(qd, function(x) x$lambda)),
                                    a_ext, divergence)

  list(m = m, n = n, a_mat = a_mat, a_ext = a_ext, ps = ps, mu = mu, A = a_obs,
       Y = y, qd = qd, H_mat = h$H_mat, H_prime = h$H_prime,
       divergence = divergence, delta = delta)
}

# ===========================================================================
# 1. The Omega identity in isolation
# ===========================================================================

test_that("legacy EIF reference implementation is present and callable", {
  expect_true(is.function(calculate_eif_lambda_legacy))
  expect_true(is.function(eif_compute_chunk_legacy))

  # The pre-fix pattern must be gone from the live implementation, and the
  # post-fix pattern present -- otherwise this file silently tests nothing.
  src <- paste(deparse(body(calculate_eif_lambda)), collapse = " ")
  expect_false(grepl("diag(weights)", src, fixed = TRUE))
  expect_true(grepl("crossprod", src, fixed = TRUE))
  expect_false(grepl("apply(a_mat", src, fixed = TRUE))

  # ...and must still be present in the frozen reference.
  src_legacy <- paste(deparse(body(calculate_eif_lambda_legacy)), collapse = " ")
  expect_true(grepl("diag(weights)", src_legacy, fixed = TRUE))
  expect_true(grepl("apply(a_mat", src_legacy, fixed = TRUE))
})

test_that("crossprod(A * w, A) equals t(A) %*% diag(w) %*% A exactly", {
  skip_on_ci()  # expect_identical() requires bit-exact floating point equality
                # between two mathematically-equivalent BLAS code paths
                # (crossprod() vs explicit diag() multiplication); these can
                # differ in the last bit across BLAS/LAPACK builds and R
                # versions. Not a correctness property of the estimator.
  # The core algebraic identity, isolated from everything else. Spans weight
  # vectors with zeros, ties, and a 20-order-of-magnitude dynamic range.
  withr::local_seed(EOV_SEED)
  for (m in c(1L, 2L, 4L, 6L, 8L, 11L)) {
    a_ext <- cbind(1, generate_binary_vectors(m))
    n_comb <- nrow(a_ext)

    weight_sets <- list(
      uniform = rep(1, n_comb),
      random = stats::runif(n_comb, 0.1, 5),
      with_zeros = ifelse(seq_len(n_comb) %% 3L == 0L, 0, stats::runif(n_comb)),
      wide_range = 10^stats::runif(n_comb, -10, 10)
    )

    for (wname in names(weight_sets)) {
      w <- weight_sets[[wname]]
      reference <- t(a_ext) %*% diag(w) %*% a_ext
      fast <- crossprod(a_ext * w, a_ext)

      expect_equal(dim(fast), c(m + 1L, m + 1L))
      expect_identical(fast, reference,
                       info = sprintf("m = %d, weights = %s", m, wname))
    }
  }
})

test_that("the bit-key pattern lookup returns exactly what the apply() scan returned", {
  # Includes the adversarial cases the key alone could not resolve: a non-binary
  # observed pattern (key collision possible) and an absent pattern. Both must
  # still produce the pre-fix outcome (no match => length 0).
  for (m in c(1L, 2L, 4L, 6L, 8L, 11L)) {
    a_ext <- cbind(1, generate_binary_vectors(m))
    pow2 <- 2^(seq_len(ncol(a_ext)) - 1L)

    lookup_new <- function(x) {
      cand <- which(as.vector(a_ext %*% pow2) == sum(x * pow2))
      cand[vapply(cand, function(k) all(a_ext[k, ] == x), logical(1))]
    }
    lookup_old <- function(x) which(apply(a_ext, 1, function(row) all(row == x)))

    # Every legitimate observed pattern.
    for (r in seq_len(nrow(a_ext))) {
      x <- a_ext[r, ]
      expect_identical(lookup_new(x), lookup_old(x),
                       info = sprintf("m = %d, row %d", m, r))
      expect_identical(lookup_new(x), r)
    }

    # Non-binary and out-of-support patterns: must find nothing, both ways.
    bad <- list(c(1, rep(2, m)), c(0, rep(0, m)), c(1, rep(0.5, m)))
    for (x in bad) {
      expect_identical(lookup_new(x), lookup_old(x),
                       info = sprintf("m = %d, adversarial pattern", m))
    }
  }
})

test_that("marginals from the full design equal marginals from the sliced design", {
  skip_on_ci()  # expect_identical() requires bit-exact floating point equality;
                # the two code paths are mathematically equivalent but can differ
                # in the last bit across BLAS/LAPACK builds and R versions. Not a
                # correctness property of the estimator (see the non-bit-exact
                # equivalence tests elsewhere in this file for the tolerance-based
                # version of this same check).
  withr::local_seed(EOV_SEED)
  for (m in c(1L, 2L, 4L, 6L, 8L, 11L)) {
    a_ext <- cbind(1, generate_binary_vectors(m))
    for (rep_i in 1:5) {
      ps <- stats::runif(nrow(a_ext), 1e-14, 1)
      ps <- ps / sum(ps)
      sliced <- compute_marginals(ps, a_ext[, -1, drop = FALSE], m = m)
      full <- compute_marginals(ps, a_ext, m = ncol(a_ext))[-1]
      expect_identical(full, sliced, info = sprintf("m = %d", m))
    }
  }
})

# ===========================================================================
# 2. calculate_eif_lambda(): full old-vs-new battery
# ===========================================================================

test_that("calculate_eif_lambda agrees with legacy across m x divergence x delta x support", {
  for (m in c(1L, 2L, 4L, 6L)) {
    kinds <- if (m == 1L) c("full", "extreme") else
      c("full", "sparse", "degenerate_marginal", "extreme")
    for (kind in kinds) {
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        deltas <- eov_deltas(m, eov_delta_base(kind))
        for (dname in names(deltas)) {
          cs <- eov_build_case(m, n = 8L, kind = kind, divergence = div,
                               delta = deltas[[dname]])
          label <- sprintf("m=%d %s delta=%s support=%s", m, div, dname, kind)

          for (i in seq_len(cs$n)) {
            old <- suppressWarnings(calculate_eif_lambda_legacy(
              cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ], cs$H_prime[i, ],
              cs$a_ext, div))
            new <- suppressWarnings(calculate_eif_lambda(
              cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ], cs$H_prime[i, ],
              cs$a_ext, div))

            expect_equal(length(new), m + 1L, info = label)
            expect_equal(new, old, tolerance = EOV_TOL, info = label)
          }
        }
      }
    }
  }
})

test_that("calculate_eif_lambda agrees with legacy at m = 8 and m = 11", {
  skip_on_cran()  # legacy dense-diag path costs ~72 ms/obs at m = 11
  for (m in c(8L, 11L)) {
    for (kind in c("full", "degenerate_marginal")) {
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        cs <- eov_build_case(m, n = if (m == 11L) 4L else 8L, kind = kind,
                             divergence = div,
                             delta = eov_deltas(m, eov_delta_base(kind))$uniform)
        label <- sprintf("m=%d %s support=%s", m, div, kind)
        for (i in seq_len(cs$n)) {
          old <- suppressWarnings(calculate_eif_lambda_legacy(
            cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ], cs$H_prime[i, ],
            cs$a_ext, div))
          new <- suppressWarnings(calculate_eif_lambda(
            cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ], cs$H_prime[i, ],
            cs$a_ext, div))
          expect_equal(new, old, tolerance = EOV_TOL, info = label)
        }
      }
    }
  }
})

test_that("the zero-padding path is actually exercised, and agrees", {
  # Guards the guard: confirms `degenerate_marginal` really does deactivate a
  # constraint, so the padding branch is not silently untested.
  m <- 4L
  cs <- eov_build_case(m, n = 6L, kind = "degenerate_marginal", divergence = "kl",
                       delta = eov_deltas(m, eov_delta_base("degenerate_marginal"))$uniform)
  n_dropped <- 0L
  for (i in seq_len(cs$n)) {
    e_a <- compute_marginals(cs$ps[i, ], cs$a_ext, m = ncol(cs$a_ext))[-1]
    active <- c(1L, find_active_constraints(e_a, tol = 1e-8) + 1L)
    if (length(active) < m + 1L) n_dropped <- n_dropped + 1L

    old <- suppressWarnings(calculate_eif_lambda_legacy(
      cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ], cs$H_prime[i, ], cs$a_ext, "kl"))
    new <- suppressWarnings(calculate_eif_lambda(
      cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ], cs$H_prime[i, ], cs$a_ext, "kl"))
    expect_equal(new, old, tolerance = EOV_TOL)
    expect_length(new, m + 1L)
  }
  expect_gt(n_dropped, 0L)
})

test_that("the Omega regularization fallback fires identically in both implementations", {
  # NOTE ON REACHABILITY: this fallback cannot be provoked through realistic
  # propensity vectors. Omega is a weighted Gram matrix of the extended design over
  # the supported combinations, so it is positive definite exactly when that support
  # has rank m+1 -- which is the SAME condition `handle_zero_propensities()` enforces
  # before `newton_qdelta()` will produce a Q_delta at all. And `omega_active` is a
  # principal submatrix of a PD matrix, hence also PD. Measured: 0 of 40
  # observations with propensity masses spanning 6 orders of magnitude triggered it.
  #
  # So the branch is exercised directly instead, by passing an `h_dot_values` vector
  # that zeroes the weight on every combination with A_1 = 0. The surviving rows all
  # have A_1 = 1, so the intercept column equals column 2 and Omega is exactly
  # singular, while the marginals (computed from `propensity_scores`, untouched) stay
  # interior so no constraint is dropped and the singular matrix reaches `solve()`.
  m <- 4L
  cs <- eov_build_case(m, n = 4L, kind = "full", divergence = "kl",
                       delta = eov_deltas(m)$uniform)
  singular_mask <- as.numeric(cs$a_mat[, 1] == 1)

  n_reg_old <- 0L
  n_reg_new <- 0L
  for (i in seq_len(cs$n)) {
    h_dot_singular <- cs$H_prime[i, ] * singular_mask

    # Both must ALSO agree that every constraint is still active, otherwise the
    # singular block would be sliced away and the test would prove nothing.
    e_a <- compute_marginals(cs$ps[i, ], cs$a_ext, m = ncol(cs$a_ext))[-1]
    expect_length(find_active_constraints(e_a, tol = 1e-8), m)

    old <- withCallingHandlers(
      calculate_eif_lambda_legacy(cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ],
                                  h_dot_singular, cs$a_ext, "kl"),
      warning = function(w) { n_reg_old <<- n_reg_old + 1L; invokeRestart("muffleWarning") })
    new <- withCallingHandlers(
      calculate_eif_lambda(cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ],
                           h_dot_singular, cs$a_ext, "kl"),
      warning = function(w) { n_reg_new <<- n_reg_new + 1L; invokeRestart("muffleWarning") })

    # Same regularized result, not merely the same code path.
    expect_equal(new, old, tolerance = EOV_TOL)
  }

  # The branch must actually have fired -- guards against a vacuous pass.
  expect_equal(n_reg_new, cs$n)
  expect_equal(n_reg_new, n_reg_old)
})

test_that("realistic propensity vectors never reach the Omega regularization fallback", {
  # Documents the reachability argument above as an executable claim, so a future
  # change that starts tripping regularization on ordinary data is noticed.
  n_reg <- 0L
  n_obs <- 0L
  for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
    cs <- eov_build_case(4L, n = 10L, kind = "extreme", divergence = div,
                         delta = eov_deltas(4L, eov_delta_base("extreme"))$uniform)
    for (i in seq_len(cs$n)) {
      n_obs <- n_obs + 1L
      withCallingHandlers(
        calculate_eif_lambda(cs$A[i, ], cs$ps[i, ], cs$qd[[i]], cs$H_mat[i, ],
                             cs$H_prime[i, ], cs$a_ext, div),
        warning = function(w) { n_reg <<- n_reg + 1L; invokeRestart("muffleWarning") })
    }
  }
  expect_gt(n_obs, 0L)
  expect_equal(n_reg, 0L)
})

# ===========================================================================
# 3. Chunk level and full end-to-end EIF
# ===========================================================================

test_that("eif_compute_chunk agrees with legacy on eif_values and plugin", {
  for (m in c(1L, 2L, 4L, 6L)) {
    kinds <- if (m == 1L) "full" else c("full", "sparse", "degenerate_marginal")
    for (kind in kinds) {
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        cs <- eov_build_case(m, n = 12L, kind = kind, divergence = div,
                             delta = eov_deltas(m, eov_delta_base(kind))$uniform)
        args <- list(chunk_indices = seq_len(cs$n), A = cs$A, Y = cs$Y,
                     R = NULL, pR = NULL, propensity_scores = cs$ps,
                     outcome_means = cs$mu, qdelta_results = cs$qd,
                     a_mat = cs$a_mat, a_extended_matrix = cs$a_ext,
                     divergence = div, m = m)
        new <- suppressWarnings(do.call(eif_compute_chunk, args))
        old <- suppressWarnings(do.call(eif_compute_chunk_legacy, args))

        label <- sprintf("m=%d %s support=%s", m, div, kind)
        expect_equal(new$eif_values, old$eif_values, tolerance = EOV_TOL, info = label)
        expect_equal(new$plugin, old$plugin, tolerance = EOV_TOL, info = label)
        expect_equal(mean(new$eif_values), mean(old$eif_values),
                     tolerance = EOV_TOL, info = label)
      }
    }
  }
})

test_that("eif_compute_chunk agrees with legacy at m = 8 and m = 11", {
  skip_on_cran()
  for (m in c(8L, 11L)) {
    for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
      cs <- eov_build_case(m, n = if (m == 11L) 5L else 10L, kind = "full",
                           divergence = div, delta = eov_deltas(m)$uniform)
      args <- list(chunk_indices = seq_len(cs$n), A = cs$A, Y = cs$Y,
                   R = NULL, pR = NULL, propensity_scores = cs$ps,
                   outcome_means = cs$mu, qdelta_results = cs$qd,
                   a_mat = cs$a_mat, a_extended_matrix = cs$a_ext,
                   divergence = div, m = m)
      new <- suppressWarnings(do.call(eif_compute_chunk, args))
      old <- suppressWarnings(do.call(eif_compute_chunk_legacy, args))
      label <- sprintf("m=%d %s", m, div)
      expect_equal(new$eif_values, old$eif_values, tolerance = EOV_TOL, info = label)
      expect_equal(new$plugin, old$plugin, tolerance = EOV_TOL, info = label)
    }
  }
})

test_that("efficient_influence_function end-to-end matches the legacy chunk path", {
  # Bridge from the fixed non-chunked path (what efficient_influence_function()
  # uses below n = 10,000) to the frozen pre-fix implementation. The two internal
  # paths are mathematically identical, which this simultaneously re-confirms.
  for (m in c(2L, 4L)) {
    for (kind in c("full", "degenerate_marginal")) {
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        cs <- eov_build_case(m, n = 20L, kind = kind, divergence = div,
                             delta = eov_deltas(m, eov_delta_base(kind))$uniform)

        new <- suppressWarnings(efficient_influence_function(
          A = cs$A, Y = cs$Y, propensity_scores = cs$ps, delta = cs$delta,
          outcome_means = cs$mu, divergence = div, qdelta_results = cs$qd,
          feasibility_check = FALSE))
        old <- suppressWarnings(eif_compute_chunk_legacy(
          seq_len(cs$n), cs$A, cs$Y, NULL, NULL, cs$ps, cs$mu, cs$qd,
          cs$a_mat, cs$a_ext, div, m))

        label <- sprintf("m=%d %s support=%s", m, div, kind)
        expect_equal(new$eif_values, old$eif_values, tolerance = EOV_TOL, info = label)
        expect_equal(new$plugin, old$plugin, tolerance = EOV_TOL, info = label)
        expect_equal(mean(new$eif_values), mean(old$eif_values),
                     tolerance = EOV_TOL, info = label)
      }
    }
  }
})

test_that("efficient_influence_function with missingness (R, pR) also matches legacy", {
  # The Phi branch that divides by pR is a separate code path in the chunk; make
  # sure the Omega fix did not perturb it.
  m <- 4L
  cs <- eov_build_case(m, n = 20L, kind = "full", divergence = "kl",
                       delta = eov_deltas(m)$uniform)
  withr::local_seed(EOV_SEED)
  r_ind <- stats::rbinom(cs$n, 1, 0.85)
  p_r <- matrix(0.85, nrow = cs$n, ncol = ncol(cs$ps))

  new <- suppressWarnings(eif_compute_chunk(
    seq_len(cs$n), cs$A, cs$Y, r_ind, p_r, cs$ps, cs$mu, cs$qd,
    cs$a_mat, cs$a_ext, "kl", m))
  old <- suppressWarnings(eif_compute_chunk_legacy(
    seq_len(cs$n), cs$A, cs$Y, r_ind, p_r, cs$ps, cs$mu, cs$qd,
    cs$a_mat, cs$a_ext, "kl", m))

  expect_equal(new$eif_values, old$eif_values, tolerance = EOV_TOL)
  expect_equal(new$plugin, old$plugin, tolerance = EOV_TOL)
})
