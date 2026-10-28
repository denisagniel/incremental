#' @file test-qdelta-vectorized-equivalence.R
#' @title Old-vs-new numerical equivalence for the vectorized Q_delta solver
#'
#' Spec: `quality_reports/specs/2026-08-10_qdelta-solver-performance-fix.md`,
#' MUST item "Numerical equivalence, old vs. new" (stage 1: vectorization).
#'
#' "Old" = `R/fast_qdelta_solver_legacy.R`, a verbatim in-tree snapshot of the
#' post-bugfix scalar-loop solver (both prerequisite fixes of 2026-08-10 already
#' applied). "New" = the vectorized `R/fast_qdelta_solver.R`. Nothing here reads
#' git history.
#'
#' Coverage required by the spec and implemented below:
#'   * m in {2, 4, 6, 8, 11}
#'   * all four divergences: kl, reverse_kl, chi2, hellinger
#'   * delta: (a) uniform, (b) sparse (some components exactly 1)
#'   * propensity vectors with sub-PROPENSITY_ACTIVE_THRESHOLD combinations
#'     (exercising `active_idx` construction) and rank-deficient support
#'     (exercising `handle_zero_propensities()`'s qr rank check). Per the spec's
#'     explicit correction of an earlier draft: this is driven by the PROPENSITY
#'     VECTOR, never by delta -- delta only moves `target_marginals`.
#'   * compared quantities: converged `lambda`, `Q_delta` weights, the
#'     EIF-relevant intermediates h/h' (`compute_H_vectorized_general()`, which is
#'     what `eif.R` actually consumes downstream of `lambda`), and
#'     `psi_delta`-style estimates (plugin `sum(q_delta * mu)` plus a full
#'     `efficient_influence_function()` pass driven by each implementation's own
#'     `qdelta_results`).
#'   * two documented exceptions to the 1e-6 bar (chi2 boundary relaxation;
#'     discrete branch flips) -- see the two dedicated blocks at the end.

VIE_TOL <- 1e-6          # repo's point-estimate tolerance (.claude/rules/quality-gates.md)
VIE_SEED <- 20260811L

# ---------------------------------------------------------------------------
# Case construction
# ---------------------------------------------------------------------------

#' Build a deterministic propensity vector of a requested support shape
#'
#' @param m Integer, number of treatments
#' @param kind One of:
#'   * `"full"` -- every one of the 2^m combinations is above
#'     PROPENSITY_ACTIVE_THRESHOLD (the production regime).
#'   * `"sparse_support"` -- only the m+1 combinations {0, e_1, ..., e_m} carry
#'     mass; the other 2^m - m - 1 are exactly 0, i.e. far below the active
#'     threshold. Rank is exactly m+1, the minimum `handle_zero_propensities()`
#'     accepts, so this is the tightest non-degenerate `active_idx` case.
#'   * `"half_support"` -- a middling fraction active (full rank, many inactive).
#'   * `"rank_deficient"` -- mass only on combinations with A_1 = 1, so column 1
#'     of the design is constant and equals the intercept: rank m < m+1. Both
#'     implementations must ERROR identically here.
#' @return Numeric vector of length 2^m summing to 1
vie_make_ps <- function(m, kind) {
  withr::local_seed(VIE_SEED + m)
  b <- generate_binary_vectors(m)
  n_comb <- nrow(b)
  p <- stats::runif(n_comb, 0.3, 1)

  if (kind == "full") {
    # keep all
  } else if (kind == "sparse_support") {
    p[rowSums(b) > 1] <- 0
  } else if (kind == "half_support") {
    keep <- rowSums(b) <= 1 | (seq_len(n_comb) %% 3L == 0L)
    p[!keep] <- 0
  } else if (kind == "rank_deficient") {
    p[b[, 1] != 1] <- 0
  } else {
    stop("unknown kind: ", kind)
  }

  p / sum(p)
}

#' Uniform and sparse delta vectors for a given m
#'
#' Sparse delta sets some components to exactly 1 (no tilt on those treatments),
#' which the spec calls out as a distinct case from uniform delta.
#'
#' @param m Integer, number of treatments
#' @param base Numeric, the tilt applied to the non-unit components. Must be
#'   chosen so the induced target marginals are FEASIBLE for the support in
#'   question: under `"sparse_support"` only \{0, e_1, ..., e_m\} carry mass, so at
#'   most one treatment can be on at a time and \eqn{\sum_k E_Q[A_k] \le 1} is a
#'   hard constraint -- any `base > 1` violates it once `m` is large enough.
#'   Callers therefore pass `base < 1` (odds shrinkage) for that support. This is
#'   a property of the sparse support, not of either implementation, and both
#'   implementations fail identically on infeasible deltas (see
#'   `vie_expect_equivalent()`, which asserts paired errors).
vie_deltas <- function(m, base = 1.2) {
  sparse <- rep(base, m)
  sparse[1] <- 1
  if (m >= 3L) sparse[m] <- 1
  list(uniform = rep(base, m), sparse = sparse)
}

#' The delta base that is feasible for a given support shape
vie_delta_base <- function(kind) {
  if (kind == "sparse_support") 0.8 else 1.2
}

# ---------------------------------------------------------------------------
# Comparison machinery
# ---------------------------------------------------------------------------

# Accumulator for cases that flipped a DISCRETE solver branch between the two
# implementations (spec exception (ii)). These must be logged and inspected, not
# silently passed and not silently failed, so they are collected here and
# asserted on in a dedicated block at the end of the file.
vie_flips <- new.env(parent = emptyenv())
vie_flips$rows <- list()

vie_log_flip <- function(label, what, old_val, new_val) {
  vie_flips$rows[[length(vie_flips$rows) + 1L]] <- sprintf(
    "%s: %s differs (legacy = %s, vectorized = %s)",
    label, what, format(old_val), format(new_val)
  )
}

#' Relative difference used for the 1e-6 bar
#'
#' Relative to the legacy magnitude, floored at 1 so that quantities legitimately
#' near zero (e.g. a lambda component of a no-tilt treatment) are compared on an
#' absolute scale rather than blowing up the ratio.
vie_rel_diff <- function(new_val, old_val) {
  max(abs(new_val - old_val) / pmax(abs(old_val), 1))
}

#' Compare one (m, divergence, delta, support) case, old vs. new
#'
#' Asserts equality of every spec-listed quantity within `VIE_TOL`, EXCEPT for
#' the two documented exceptions, which are routed to their own handling:
#'   (i)  chi2 cases the legacy run itself accepted via boundary relaxation
#'        (`constraint_norm < 0.1`): the bar becomes "the vectorized run does at
#'        least as well as the legacy `constraint_norm` it achieved", not 1e-6 on
#'        lambda, because 0.1 is the legacy contract for those cases.
#'   (ii) discrete branch flips: logged for inspection via `vie_log_flip()`.
#'
#' @param m,divergence,delta,ps The case
#' @param label Character label used in failure messages
vie_expect_equivalent <- function(m, divergence, delta, ps, label) {
  old <- try(newton_qdelta_legacy(ps, delta, divergence, fix_lambda0 = FALSE, tol = 1e-10),
             silent = TRUE)
  new <- try(newton_qdelta(ps, delta, divergence, fix_lambda0 = FALSE, tol = 1e-10),
             silent = TRUE)

  # --- Paired failure is itself an equivalence result ---------------------
  # Some (support, delta) combinations are genuinely infeasible (e.g. a minimal
  # rank-(m+1) support cannot realise target marginals summing above 1). Both
  # implementations must then refuse in the same way; a rewrite that "fixes" an
  # infeasible case, or breaks a feasible one, is a behaviour change.
  expect_identical(inherits(old, "try-error"), inherits(new, "try-error"),
                   info = paste(label, "- one implementation errored and the other did not"))
  if (inherits(old, "try-error")) {
    expect_identical(conditionMessage(attr(new, "condition")),
                     conditionMessage(attr(old, "condition")),
                     info = paste(label, "- error messages differ"))
    return(invisible(NULL))
  }

  # --- Structural contract ------------------------------------------------
  expect_equal(length(new$lambda), length(old$lambda), info = label)
  expect_equal(length(new$q_delta_weights), length(old$q_delta_weights), info = label)

  # --- Exception (ii): discrete branch flips ------------------------------
  # Iteration COUNTS need not match (floating-point reordering can shift the
  # path), but a differing `converged` flag means a different discrete
  # convergence branch fired. Log it.
  if (!identical(isTRUE(old$converged), isTRUE(new$converged))) {
    vie_log_flip(label, "converged flag", old$converged, new$converged)
  }

  # --- Exception (i): chi2 boundary relaxation ---------------------------
  # Legacy accepted convergence despite constraint_norm >= tol because the
  # solution sat on the chi2 domain boundary and the norm was < 0.1.
  legacy_used_chi2_relaxation <-
    divergence == "chi2" && isTRUE(old$converged) &&
    !is.na(old$constraint_norm) && old$constraint_norm >= 1e-10

  # ...or legacy did not converge at all: then the contract is likewise
  # "at least as good", since there is no converged value to match.
  legacy_unconverged <- !isTRUE(old$converged)

  if (legacy_used_chi2_relaxation || legacy_unconverged) {
    expect_true(
      is.na(old$constraint_norm) ||
        (!is.na(new$constraint_norm) &&
           new$constraint_norm <= old$constraint_norm * (1 + 1e-6)),
      info = paste(label, "- vectorized constraint_norm must be <= legacy's",
                   "achieved norm (legacy =", format(old$constraint_norm),
                   ", vectorized =", format(new$constraint_norm), ")")
    )
    return(invisible(NULL))
  }

  # --- The 1e-6 bar: converged lambda ------------------------------------
  expect_lt(vie_rel_diff(new$lambda, old$lambda), VIE_TOL)

  # --- The 1e-6 bar: Q_delta weights ------------------------------------
  expect_lt(vie_rel_diff(new$q_delta_weights, old$q_delta_weights), VIE_TOL)

  # --- The 1e-6 bar: EIF-relevant intermediates -------------------------
  # eif.R consumes lambda only through compute_H_vectorized_general(), so this
  # is the quantity that actually propagates into the influence function.
  a_ext <- cbind(1, generate_binary_vectors(m))
  h_old <- compute_H_vectorized_general(matrix(old$lambda, nrow = 1), a_ext, divergence)
  h_new <- compute_H_vectorized_general(matrix(new$lambda, nrow = 1), a_ext, divergence)
  expect_lt(vie_rel_diff(as.numeric(h_new$H_mat), as.numeric(h_old$H_mat)), VIE_TOL)
  expect_lt(vie_rel_diff(as.numeric(h_new$H_prime), as.numeric(h_old$H_prime)), VIE_TOL)

  # --- The 1e-6 bar: psi_delta-style plugin estimate ---------------------
  # Deterministic outcome-mean surface so this is a pure function of q_delta.
  mu <- 0.5 + 0.3 * rowSums(a_ext[, -1, drop = FALSE]) / m
  psi_old <- sum(old$q_delta_weights * mu)
  psi_new <- sum(new$q_delta_weights * mu)
  expect_lt(abs(psi_new - psi_old) / max(abs(psi_old), 1), VIE_TOL)

  # Constraint norm must not degrade.
  expect_true(
    !is.na(new$constraint_norm) &&
      new$constraint_norm <= max(old$constraint_norm * (1 + 1e-6), 1e-9),
    info = paste(label, "- constraint_norm degraded:",
                 format(old$constraint_norm), "->", format(new$constraint_norm))
  )

  invisible(NULL)
}

# ===========================================================================
# 1. Core rewrite: the Jacobian assembly itself
# ===========================================================================

test_that("legacy scalar-loop reference implementation is present and callable", {
  # Guards against the legacy snapshot being deleted or renamed, which would
  # silently turn this whole file into a no-op.
  expect_true(is.function(newton_qdelta_legacy))
  expect_true(is.function(compute_jacobian_analytical_internal_legacy))
  expect_true(is.function(compute_marginal_jacobian_internal_legacy))
  expect_true(is.function(handle_zero_propensities_legacy))

  # And that the vectorized version really did gain the reuse hooks.
  expect_true(all(c("b_ext_active", "moment_scatter", "row_perm") %in%
                    names(formals(compute_jacobian_analytical_internal))))
  expect_true("combination_matrix" %in% names(formals(newton_qdelta)))
  expect_true("combination_matrix" %in% names(formals(handle_zero_propensities)))
  expect_true("combination_matrix" %in% names(formals(compute_qdelta_per_obs)))
})

test_that("vectorized (m+1)x(m+1) Jacobian equals the legacy (i,j)-loop Jacobian", {
  # Tightest test of the rewrite: identical inputs, no Newton iteration in the
  # way, so any algebra error in the crossprod/row-permutation form shows up
  # immediately rather than being masked by re-convergence.
  for (m in c(2L, 4L, 6L)) {
    vectors <- generate_binary_vectors(m)
    p0 <- vie_make_ps(m, "full")

    lambda_grid <- list(
      kl         = list(c(0, rep(0.1, m)), c(0.3, seq(-0.2, 0.2, length.out = m))),
      reverse_kl = list(c(-1, rep(-0.2, m)), c(-2.5, rep(-0.05, m))),
      chi2       = list(c(-0.5, rep(0.1, m)), c(-0.2, rep(-0.03, m))),
      hellinger  = list(c(-0.1, rep(0.05, m)), c(-0.6, rep(-0.02, m)))
    )

    for (div in names(lambda_grid)) {
      for (lambda in lambda_grid[[div]]) {
        legacy <- compute_jacobian_analytical_internal_legacy(lambda, p0, vectors, div)
        vect <- compute_jacobian_analytical_internal(lambda, p0, vectors, div)

        expect_equal(dim(vect), c(m + 1L, m + 1L))
        expect_equal(vect, legacy, tolerance = 1e-10,
                     info = sprintf("m = %d, divergence = %s", m, div))
      }
    }
  }
})

test_that("single-cell Jacobian accessors agree with the assembled matrix and with legacy", {
  # The cell-level helpers are retained for diagnostics; they must still return
  # the same numbers as before and as the assembled matrix.
  m <- 4L
  vectors <- generate_binary_vectors(m)
  p0 <- vie_make_ps(m, "full")
  lambda <- c(-0.3, seq(-0.15, 0.15, length.out = m))

  for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
    assembled <- compute_jacobian_analytical_internal(lambda, p0, vectors, div)
    for (i in seq_len(m)) {
      for (j in seq_len(m + 1L)) {
        expect_equal(
          compute_marginal_jacobian_internal(i, j, lambda, p0, vectors, div),
          compute_marginal_jacobian_internal_legacy(i, j, lambda, p0, vectors, div),
          tolerance = 1e-10
        )
        expect_equal(
          compute_marginal_jacobian_internal(i, j, lambda, p0, vectors, div),
          assembled[i, j], tolerance = 1e-10
        )
      }
    }
    for (j in seq_len(m + 1L)) {
      expect_equal(
        compute_normalization_jacobian_internal(j, lambda, p0, vectors, div),
        compute_normalization_jacobian_internal_legacy(j, lambda, p0, vectors, div),
        tolerance = 1e-10
      )
      expect_equal(
        compute_normalization_jacobian_internal(j, lambda, p0, vectors, div),
        assembled[m + 1L, j], tolerance = 1e-10
      )
    }
  }
})

test_that("chi2's moment scatter is lambda-invariant, justifying the Newton-loop hoist", {
  # The hoist is only valid because h'(x) == 0.5 for every x under chi2. If a
  # future change makes chi2's h' non-constant, this test fails loudly instead
  # of the solver silently using a stale Jacobian.
  m <- 5L
  vectors <- generate_binary_vectors(m)
  b_ext <- cbind(1, vectors)
  p0 <- vie_make_ps(m, "full")

  reference <- crossprod(b_ext, as.vector(p0 * 0.5) * b_ext)

  for (lambda in list(c(0, rep(0, m)), c(-0.5, rep(0.1, m)), c(-1.9, rep(-0.01, m)))) {
    expect_equal(compute_moment_scatter_internal(lambda, p0, b_ext, "chi2"),
                 reference, tolerance = 1e-12)
  }

  # ...and the non-constant-h' divergences must NOT be lambda-invariant, i.e.
  # the hoist must not be extended to them by accident.
  for (div in c("kl", "reverse_kl", "hellinger")) {
    s_a <- compute_moment_scatter_internal(c(-0.5, rep(0.1, m)), p0, b_ext, div)
    s_b <- compute_moment_scatter_internal(c(-1.0, rep(-0.1, m)), p0, b_ext, div)
    expect_false(isTRUE(all.equal(s_a, s_b, tolerance = 1e-8)))
  }
})

test_that("vectorized external-normalization (fix_lambda0) Jacobian equals legacy", {
  # fix_lambda0 = FALSE is the only path production uses, but fix_lambda0 is a
  # documented argument, so its m x m Jacobian was vectorized too and must match.
  for (m in c(2L, 4L)) {
    vectors <- generate_binary_vectors(m)
    p0 <- vie_make_ps(m, "full")
    for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
      lambda_active <- switch(div,
        kl = rep(0.12, m), reverse_kl = rep(-0.2, m),
        chi2 = rep(0.05, m), hellinger = rep(0.04, m))
      lambda0 <- if (div == "reverse_kl") -1.0 else if (div == "hellinger") -0.3 else 0.0

      expect_equal(
        compute_jacobian_analytical_m(lambda_active, p0, vectors, div, lambda0),
        compute_jacobian_analytical_m_legacy(lambda_active, p0, vectors, div, lambda0),
        tolerance = 1e-10, info = sprintf("m = %d, divergence = %s", m, div)
      )
    }
  }
})

test_that("vectorized domain checks agree with the legacy scalar early-return loops", {
  for (m in c(2L, 4L, 6L)) {
    vectors <- generate_binary_vectors(m)
    p0 <- vie_make_ps(m, "full")
    # Deliberately spans feasible, marginal and infeasible lambdas per divergence.
    lambdas <- list(
      c(0, rep(0, m)), c(-1, rep(-0.2, m)), c(-0.5, rep(0.1, m)),
      c(0.9, rep(0.05, m)), c(-2, rep(-0.5, m)), c(1.5, rep(-0.1, m)),
      c(-3, rep(0.4, m))
    )
    for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
      for (lambda in lambdas) {
        expect_identical(
          check_constraint_satisfaction_internal(lambda, p0, vectors, div),
          check_constraint_satisfaction_internal_legacy(lambda, p0, vectors, div),
          info = sprintf("m = %d, %s, lambda0 = %g", m, div, lambda[1])
        )
        expect_identical(
          check_constraint_satisfaction_m(lambda[-1], p0, vectors, div, lambda[1]),
          check_constraint_satisfaction_m_legacy(lambda[-1], p0, vectors, div, lambda[1]),
          info = sprintf("m = %d, %s (external), lambda0 = %g", m, div, lambda[1])
        )
      }
    }
  }
})

test_that("vectorized initialization agrees with the legacy sapply initialization", {
  for (m in c(2L, 4L, 6L)) {
    for (kind in c("full", "sparse_support")) {
      ps <- vie_make_ps(m, kind)
      active <- which(ps > PROPENSITY_ACTIVE_THRESHOLD)
      vectors <- generate_binary_vectors(m)[active, , drop = FALSE]
      p0 <- ps[active]
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        for (dl in vie_deltas(m)) {
          expect_equal(
            get_initial_lambda_internal(dl, p0, vectors, div),
            get_initial_lambda_internal_legacy(dl, p0, vectors, div),
            tolerance = 1e-12
          )
          expect_equal(
            get_initial_lambda_active(dl, p0, vectors, div),
            get_initial_lambda_active_legacy(dl, p0, vectors, div),
            tolerance = 1e-12
          )
        }
      }
    }
  }
})

# ===========================================================================
# 2. Full-solve equivalence battery
# ===========================================================================

test_that("full solves agree within 1e-6: m in {2, 4, 6} x 4 divergences x delta x support", {
  for (m in c(2L, 4L, 6L)) {
    supports <- if (m == 2L) "full" else c("full", "sparse_support", "half_support")
    for (kind in supports) {
      ps <- vie_make_ps(m, kind)
      deltas <- vie_deltas(m, vie_delta_base(kind))
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        for (dname in names(deltas)) {
          label <- sprintf("m=%d %s delta=%s support=%s", m, div, dname, kind)
          vie_expect_equivalent(m, div, deltas[[dname]], ps, label)
        }
      }
    }
  }
})

test_that("full solves agree within 1e-6 at m = 8 (all four divergences)", {
  skip_on_cran()  # legacy scalar path costs ~0.5-1 s per solve at 2^8 combinations
  m <- 8L
  for (kind in c("full", "sparse_support")) {
    ps <- vie_make_ps(m, kind)
    deltas <- vie_deltas(m, vie_delta_base(kind))
    for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
      for (dname in names(deltas)) {
        label <- sprintf("m=8 %s delta=%s support=%s", div, dname, kind)
        vie_expect_equivalent(m, div, deltas[[dname]], ps, label)
      }
    }
  }
})

test_that("full solves agree within 1e-6 at m = 11 (all four divergences)", {
  skip_on_cran()  # legacy scalar path costs ~7-14 s per full-support solve at 2^11
  m <- 11L

  # Sparse support (n_active = m + 1 = 12) is cheap in the legacy path, so the
  # full divergence x delta grid runs there.
  ps_sparse <- vie_make_ps(m, "sparse_support")
  deltas_sparse <- vie_deltas(m, vie_delta_base("sparse_support"))
  for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
    for (dname in names(deltas_sparse)) {
      vie_expect_equivalent(m, div, deltas_sparse[[dname]], ps_sparse,
                            sprintf("m=11 %s delta=%s support=sparse_support", div, dname))
    }
  }

  # Full support (all 2048 combinations active) is the production-shaped case and
  # the one the speedup claim is about; one delta per divergence keeps the legacy
  # reference affordable.
  ps_full <- vie_make_ps(m, "full")
  delta_full <- vie_deltas(m, vie_delta_base("full"))$uniform
  for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
    vie_expect_equivalent(m, div, delta_full, ps_full,
                          sprintf("m=11 %s delta=uniform support=full", div))
  }
})

test_that("rank-deficient support errors identically in both implementations", {
  # Spec: this path is reached via the PROPENSITY VECTOR (constant treatment
  # column => design rank m < m+1), never via delta.
  for (m in c(2L, 4L, 6L, 8L)) {
    ps <- vie_make_ps(m, "rank_deficient")
    for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
      for (delta in vie_deltas(m)) {
        expect_error(newton_qdelta_legacy(ps, delta, div, fix_lambda0 = FALSE),
                     "Under-identification")
        expect_error(newton_qdelta(ps, delta, div, fix_lambda0 = FALSE),
                     "Under-identification")
      }
    }
  }
})

test_that("all-zero propensity vector errors identically in both implementations", {
  m <- 3L
  ps <- rep(0, 2^m)
  ps[1] <- 1  # sums to 1 but only one active combination => rank 1 < m + 1
  expect_error(newton_qdelta_legacy(ps, rep(1.2, m), "kl"), "Under-identification")
  expect_error(newton_qdelta(ps, rep(1.2, m), "kl"), "Under-identification")
})

# ===========================================================================
# 3. Combination-matrix threading (reuse, not caching)
# ===========================================================================

test_that("threading combination_matrix changes nothing numerically", {
  for (m in c(2L, 4L, 6L)) {
    b <- generate_binary_vectors(m)
    for (kind in c("full", "sparse_support")) {
      ps <- vie_make_ps(m, kind)
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        delta <- vie_deltas(m, vie_delta_base(kind))$uniform
        internal <- newton_qdelta(ps, delta, div, fix_lambda0 = FALSE, tol = 1e-10)
        threaded <- newton_qdelta(ps, delta, div, fix_lambda0 = FALSE, tol = 1e-10,
                                  combination_matrix = b)
        expect_identical(threaded, internal,
                         info = sprintf("m = %d, %s, %s", m, div, kind))

        # ...and through the public per-observation wrapper.
        w_internal <- compute_qdelta_per_obs(ps, delta, div, tol = 1e-10)
        w_threaded <- compute_qdelta_per_obs(ps, delta, div, tol = 1e-10,
                                             combination_matrix = b)
        expect_identical(w_threaded, w_internal)
      }
    }
  }
})

test_that("a malformed combination_matrix is rejected rather than silently used", {
  m <- 3L
  ps <- vie_make_ps(m, "full")
  expect_error(newton_qdelta(ps, rep(1.2, m), "kl",
                             combination_matrix = generate_binary_vectors(m + 1L)),
               "combination_matrix")
  expect_error(newton_qdelta(ps, rep(1.2, m), "kl",
                             combination_matrix = generate_binary_vectors(m)[, 1:2]),
               "combination_matrix")
  expect_error(newton_qdelta(ps, rep(1.2, m), "kl",
                             combination_matrix = as.numeric(1:8)),
               "combination_matrix")
})

test_that("handle_zero_propensities agrees with legacy with and without threading", {
  for (m in c(2L, 4L, 6L)) {
    b <- generate_binary_vectors(m)
    for (kind in c("full", "sparse_support", "half_support")) {
      ps <- vie_make_ps(m, kind)
      for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
        legacy <- handle_zero_propensities_legacy(ps, div, m)
        expect_equal(handle_zero_propensities(ps, div, m), legacy, tolerance = 0)
        expect_equal(handle_zero_propensities(ps, div, m, combination_matrix = b),
                     legacy, tolerance = 0)
      }
    }
  }
})

# ===========================================================================
# 4. End-to-end: EIF values and psi_delta
# ===========================================================================

test_that("efficient_influence_function gives the same psi_delta from legacy and vectorized Q_delta", {
  # The strongest behaviour-preservation check available without touching any
  # caller: run the real EIF pass twice, once fed by each implementation's own
  # qdelta_results, and compare the influence function and the point estimate.
  withr::local_seed(VIE_SEED)
  n <- 40L

  for (m in c(2L, 4L)) {
    b <- generate_binary_vectors(m)
    n_comb <- nrow(b)

    ps <- matrix(stats::runif(n * n_comb, 0.3, 1), nrow = n)
    ps <- ps / rowSums(ps)
    mu <- matrix(stats::rnorm(n * n_comb, mean = 1, sd = 0.3), nrow = n)
    obs_idx <- apply(ps, 1, function(pr) sample.int(n_comb, 1, prob = pr))
    a_obs <- b[obs_idx, , drop = FALSE]
    y <- mu[cbind(seq_len(n), obs_idx)] + stats::rnorm(n, sd = 0.1)

    for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
      delta <- vie_deltas(m)$uniform

      qd_old <- lapply(seq_len(n), function(i) {
        r <- newton_qdelta_legacy(ps[i, ], delta, div, fix_lambda0 = FALSE, tol = 1e-10)
        list(q_delta = r$q_delta_weights, lambda = r$lambda, converged = r$converged,
             constraint_norm = r$constraint_norm, iterations = r$iterations, delta = delta)
      })
      qd_new <- lapply(seq_len(n), function(i) {
        compute_qdelta_per_obs(ps[i, ], delta, div, fix_lambda0 = FALSE, tol = 1e-10,
                               combination_matrix = b)
      })

      eif_args <- list(A = a_obs, Y = y, propensity_scores = ps, delta = delta,
                       outcome_means = mu, divergence = div, feasibility_check = FALSE)
      eif_old <- do.call(efficient_influence_function, c(eif_args, list(qdelta_results = qd_old)))
      eif_new <- do.call(efficient_influence_function, c(eif_args, list(qdelta_results = qd_new)))

      label <- sprintf("m = %d, divergence = %s", m, div)
      expect_lt(vie_rel_diff(eif_new$eif_values, eif_old$eif_values), VIE_TOL)
      expect_lt(vie_rel_diff(eif_new$plugin, eif_old$plugin), VIE_TOL)
      expect_lt(abs(mean(eif_new$eif_values) - mean(eif_old$eif_values)) /
                  max(abs(mean(eif_old$eif_values)), 1), VIE_TOL)
      expect_equal(eif_new$convergence_summary$n_converged,
                   eif_old$convergence_summary$n_converged, info = label)
    }
  }
})

# ===========================================================================
# 5. Documented exception (i): chi2 boundary relaxation
# ===========================================================================

test_that("chi2 boundary-relaxation cases are compared against legacy's own achieved norm", {
  # These deltas push chi2's linear terms onto the domain boundary
  # (linear_term = -2), where h(x) = x/2 + 1 = 0 and Newton cannot improve
  # further. The legacy contract for such cases is `constraint_norm < 0.1`
  # (fast_qdelta_solver.R's chi2-boundary-convergence block), NOT 1e-10, so the
  # equivalence bar here is "the vectorized run does at least as well", exactly
  # as the spec prescribes.
  cases <- list(
    list(m = 2L, delta = 5),    # legacy: converged = TRUE via relaxation
    list(m = 3L, delta = 2),    # legacy: converged = TRUE via relaxation
    list(m = 2L, delta = 8),    # legacy: converged = FALSE, norm ~0.12
    list(m = 4L, delta = 2)     # legacy: converged = FALSE, norm ~0.16
  )

  n_relaxation_cases <- 0L
  for (case in cases) {
    m <- case$m
    ps <- vie_make_ps(m, "full")
    delta <- rep(case$delta, m)

    old <- newton_qdelta_legacy(ps, delta, "chi2", fix_lambda0 = FALSE, tol = 1e-10)
    new <- newton_qdelta(ps, delta, "chi2", fix_lambda0 = FALSE, tol = 1e-10)

    label <- sprintf("chi2 boundary case m = %d, delta = %g", m, case$delta)

    # Confirm the case really is a boundary case (guards the test itself against
    # drifting into the ordinary-convergence regime and quietly proving nothing).
    lt_old <- as.vector(cbind(1, generate_binary_vectors(m)) %*% old$lambda)
    expect_true(any(lt_old <= -2 + 1e-6), info = paste(label, "- expected boundary contact"))
    expect_gt(old$constraint_norm, 1e-10)
    n_relaxation_cases <- n_relaxation_cases + 1L

    # The spec's exception: compare against legacy's achieved norm, not 1e-6.
    expect_lte(new$constraint_norm, old$constraint_norm * (1 + 1e-6))

    # The relaxation VERDICT must still match: if legacy accepted at < 0.1, so
    # must the vectorized run (a flip here is a discrete branch flip, logged).
    if (!identical(isTRUE(old$converged), isTRUE(new$converged))) {
      vie_log_flip(label, "chi2 relaxation verdict", old$converged, new$converged)
    }
  }

  expect_equal(n_relaxation_cases, length(cases))
})

# ===========================================================================
# 6. Documented exception (ii): discrete branch flips must be logged
# ===========================================================================

test_that("no discrete solver branch flipped between implementations (and any flip is logged)", {
  # Spec exception (ii): BLAS reordering could in principle flip one of the
  # solver's discrete branches (the Jacobian-mode auto-switch at
  # constraint_norm > prev_constraint_norm * 1.5, the periodic perturbation, or
  # the chi2 boundary-relaxation verdict) for a case sitting within
  # floating-point distance of the threshold. That is a measure-zero coincidence
  # rather than a rewrite bug, but it must not silently pass either. Every case
  # above records such flips; this block surfaces them.
  flips <- vie_flips$rows
  if (length(flips) > 0L) {
    message("Discrete branch flips requiring manual inspection:\n  ",
            paste(unlist(flips), collapse = "\n  "))
  }
  expect_equal(length(flips), 0L,
               info = paste("branch flips detected:", paste(unlist(flips), collapse = "; ")))
})
