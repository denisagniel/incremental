#' @file test-omega-reduced-reverse-kl-bugfix.R
#' @title Regression test: compute_omega_reduced()'s reverse_kl h' was wrong
#'
#' [BUGFIX 2026-08-11] `compute_omega_reduced()` (`R/treatment_subset_utils.R`)
#' computed `h'(z) = 1/(1-z)^2` for `reverse_kl`, claiming `h(z) = 1/(1-z)`. The
#' canonical, paper-derived conjugate for reverse KL is `h(x) = -1/x`,
#' `h'(x) = 1/x^2` on the domain `x < 0` (`divergence_h_pair()`,
#' `R/divergence_h_functions.R`, which every other h/h' consumer in the package
#' uses since the 2026-08-10 unification -- this was the one site that pass
#' missed).
#'
#' `linear_terms` in `compute_omega_reduced()` is `as.vector(a_bar_matrix %*% lambda)`,
#' the identical quantity `x` the canonical function takes, so the two formulas are
#' directly comparable. They disagreed by 2.25x (at x = -2) to 121x (at x = -0.1)
#' across the *entire* region the solver reaches -- this was not a boundary-only
#' defect.
#'
#' This function feeds the active-set EIF path (`R/eif_active.R`) only, which is
#' separate from the main `efficient_influence_function()` pipeline. Blast radius
#' was investigated: no stored application result was ever computed with
#' `divergence = "reverse_kl"` through this path (all three reachable call sites ran
#' with `kl`, whose branch was already correct). See
#' `quality_reports/2026-08-11_omega-reduced-reverse-kl-bugfix-report.md`.
#'
#' Scope guard: the fix touched ONLY the `reverse_kl` branch. `kl`, `chi2` and
#' `hellinger` are pinned below against explicitly-written pre-fix reference
#' formulas so an accidental edit to them would fail here.

# ---------------------------------------------------------------------------
# Reference implementations
# ---------------------------------------------------------------------------

#' Omega assembled from a caller-supplied h' vector
#'
#' Mirrors `compute_omega_reduced()`'s assembly step exactly, so any difference
#' between it and the function under test is attributable to h' alone.
orkl_omega_from_hprime <- function(propensity_scores, a_bar_matrix, h_prime) {
  weights <- propensity_scores * h_prime
  omega <- t(a_bar_matrix * weights) %*% a_bar_matrix
  # compute_omega_reduced() adds ridge regularization when near-singular; mirror it
  # so the comparison is apples-to-apples rather than tripping on that branch.
  if (min(eigen(omega, only.values = TRUE)$values) < 1e-8) {
    omega <- omega + diag(1e-6, nrow(omega))
  }
  omega
}

#' The h' formulas EXACTLY as they stood before the 2026-08-11 fix
#'
#' `reverse_kl` is the buggy one (kept here so the test can prove the function no
#' longer returns it); the other three are the pre-fix references the fix must not
#' have disturbed.
orkl_prefix_h_prime <- function(linear_terms, divergence) {
  n <- length(linear_terms)
  if (divergence == "kl") {
    exp(linear_terms)
  } else if (divergence == "chi2") {
    rep(0.5, n)
  } else if (divergence == "reverse_kl") {
    denominator <- (1 - linear_terms)^2          # <- THE BUG
    denominator <- pmax(denominator, 1e-10)
    1 / denominator
  } else if (divergence == "hellinger") {
    denominator <- (1 - linear_terms)^3
    if (any(abs(denominator) < 1e-10)) {
      denominator <- pmax(abs(denominator), 1e-10) * sign(denominator)
    }
    2 / denominator
  }
}

#' A reverse_kl-valid case: all linear terms strictly negative
#'
#' Mirrors the m = 3 case used in the stage-2 investigation that measured the 69%
#' Omega discrepancy. `lambda = c(-1, -0.2, -0.2, -0.2)` puts the linear terms in
#' [-1.6, -1], comfortably inside the domain `x < 0` and clear of the clipping
#' buffer, so nothing is saturated and the comparison isolates the formula.
orkl_case <- function(m = 3L, seed = 20260811L) {
  withr::local_seed(seed)
  a_bar <- cbind(1, generate_binary_vectors(m))
  ps <- stats::runif(nrow(a_bar), 0.3, 1)
  ps <- ps / sum(ps)
  list(a_bar = a_bar, ps = ps, lambda = c(-1, rep(-0.2, m)), m = m)
}

# ===========================================================================
# 1. The bug itself
# ===========================================================================

test_that("compute_omega_reduced's reverse_kl Omega now equals the canonical divergence_h_pair Omega", {
  for (m in c(1L, 2L, 3L, 4L, 6L)) {
    cs <- orkl_case(m)
    linear_terms <- as.vector(cs$a_bar %*% cs$lambda)

    # Sanity: the case really is inside the reverse_kl domain and unclipped, so a
    # disagreement can only come from the formula, not from saturation.
    expect_true(all(linear_terms < -REVERSE_KL_DOMAIN_BUFFER),
                info = sprintf("m = %d: case must be strictly inside x < 0", m))
    expect_equal(divergence_h_pair(linear_terms, "reverse_kl")$n_clipped, 0L)

    canonical <- orkl_omega_from_hprime(
      cs$ps, cs$a_bar, divergence_h_pair(linear_terms, "reverse_kl")$h_prime)
    actual <- compute_omega_reduced(cs$ps, cs$lambda, cs$a_bar, "reverse_kl")

    expect_equal(actual, canonical, tolerance = 1e-12,
                 info = sprintf("m = %d", m))
  }
})

test_that("the pre-fix reverse_kl formula is demonstrably NOT what is returned any more", {
  # Guards against a silent revert. Also records the magnitude of the old error so
  # the severity is not lost if this test is ever read in isolation.
  cs <- orkl_case(3L)
  linear_terms <- as.vector(cs$a_bar %*% cs$lambda)

  buggy <- orkl_omega_from_hprime(cs$ps, cs$a_bar,
                                  orkl_prefix_h_prime(linear_terms, "reverse_kl"))
  actual <- compute_omega_reduced(cs$ps, cs$lambda, cs$a_bar, "reverse_kl")

  rel_diff <- max(abs(actual - buggy) / pmax(abs(buggy), 1e-12))
  # The stage-2 investigation measured ~0.69 relative on this shape; require a
  # substantial, unmistakable difference rather than pinning an exact number.
  expect_gt(rel_diff, 0.5)
})

test_that("reverse_kl h' equals 1/x^2, not 1/(1-x)^2, across the solver-reachable domain", {
  # Formula-level check, independent of the Omega assembly: recovers h' from a
  # one-hot propensity vector, where Omega[1, 1] = pi_k * h'(x_k) = h'(x_k).
  m <- 2L
  a_bar <- cbind(1, generate_binary_vectors(m))
  for (target in c(-0.1, -0.25, -0.5, -1, -2, -5)) {
    # lambda = (target, 0, 0) => every linear term equals `target`
    lambda <- c(target, rep(0, m))
    ps <- rep(1 / nrow(a_bar), nrow(a_bar))
    omega <- compute_omega_reduced(ps, lambda, a_bar, "reverse_kl")

    # Omega[1,1] = sum_k pi_k * h'(target) * 1 * 1 = h'(target)
    expect_equal(omega[1, 1], 1 / target^2, tolerance = 1e-10,
                 info = sprintf("x = %g: expected canonical 1/x^2", target))
    expect_false(isTRUE(all.equal(omega[1, 1], 1 / (1 - target)^2, tolerance = 1e-6)))
  }
})

test_that("reverse_kl clipping warns and saturates per the canonical convention", {
  # Domain violation (x >= 0) must be pinned to the interior by divergence_h_pair()
  # and surfaced, rather than silently producing a negative or nonsense h'.
  m <- 2L
  a_bar <- cbind(1, generate_binary_vectors(m))
  ps <- rep(1 / nrow(a_bar), nrow(a_bar))
  lambda <- c(0.5, 0, 0)   # all linear terms = +0.5, outside x < 0

  expect_warning(
    omega <- compute_omega_reduced(ps, lambda, a_bar, "reverse_kl"),
    "clipped to the reverse_kl domain"
  )
  # h' saturates at 1 / buffer^2, matching divergence_h_pair() exactly.
  expect_equal(omega[1, 1], 1 / REVERSE_KL_DOMAIN_BUFFER^2, tolerance = 1e-6)
  expect_gt(omega[1, 1], 0)
})

# ===========================================================================
# 2. Scope guard: kl / chi2 / hellinger must be untouched
# ===========================================================================

test_that("kl, chi2 and hellinger Omega are UNCHANGED by the reverse_kl fix", {
  # Pinned against the pre-fix formulas written out explicitly above, so this
  # fails if a future edit routes them through divergence_h_pair() (which would
  # change them in their clipping regimes -- see the note in the test below).
  #
  # Interior lambdas only: `kl` and `hellinger` agree with canonical there but
  # legitimately differ from it once their guards engage, so the interior is where
  # "unchanged" and "canonical" coincide and the assertion is unambiguous.
  for (m in c(1L, 2L, 3L, 4L, 6L)) {
    withr::local_seed(20260811L + m)
    a_bar <- cbind(1, generate_binary_vectors(m))
    ps <- stats::runif(nrow(a_bar), 0.3, 1)
    ps <- ps / sum(ps)

    lambdas <- list(
      kl        = c(0.1, rep(0.05, m)),    # exp() well away from overflow
      chi2      = c(-0.5, rep(0.1, m)),
      hellinger = c(-0.3, rep(0.05, m))    # linear terms < 1, guard idle
    )

    for (div in names(lambdas)) {
      lambda <- lambdas[[div]]
      linear_terms <- as.vector(a_bar %*% lambda)

      prefix <- orkl_omega_from_hprime(ps, a_bar,
                                       orkl_prefix_h_prime(linear_terms, div))
      actual <- compute_omega_reduced(ps, lambda, a_bar, div)

      expect_equal(actual, prefix, tolerance = 0,
                   info = sprintf("m = %d, divergence = %s must be bit-identical to pre-fix",
                                  m, div))
    }
  }
})

test_that("kl, chi2 and hellinger also still agree with canonical in the interior", {
  # The other half of the scope claim: these three were already correct, and the
  # fix did not make them wrong. (They diverge from canonical only in their
  # clipping regimes -- kl above x = 700, hellinger as x -> 1 and for x > 1 -- a
  # separate, documented, deliberately-unfixed inconsistency.)
  m <- 3L
  withr::local_seed(20260811L)
  a_bar <- cbind(1, generate_binary_vectors(m))
  ps <- stats::runif(nrow(a_bar), 0.3, 1)
  ps <- ps / sum(ps)

  lambdas <- list(kl = c(0.1, rep(0.05, m)),
                  chi2 = c(-0.5, rep(0.1, m)),
                  hellinger = c(-0.3, rep(0.05, m)))

  for (div in names(lambdas)) {
    lambda <- lambdas[[div]]
    linear_terms <- as.vector(a_bar %*% lambda)
    canonical <- orkl_omega_from_hprime(
      ps, a_bar, divergence_h_pair(linear_terms, div)$h_prime)
    expect_equal(compute_omega_reduced(ps, lambda, a_bar, div), canonical,
                 tolerance = 1e-12, info = div)
  }
})

test_that("all four divergences still return a positive definite Omega", {
  # The function's documented contract; must survive the fix.
  for (m in c(2L, 3L, 4L)) {
    cs <- orkl_case(m)
    for (div in c("kl", "chi2", "hellinger", "reverse_kl")) {
      omega <- compute_omega_reduced(cs$ps, cs$lambda, cs$a_bar, div)
      expect_equal(dim(omega), c(m + 1L, m + 1L))
      expect_true(all(eigen(omega, only.values = TRUE)$values > 0),
                  info = sprintf("m = %d, %s", m, div))
      expect_equal(omega, t(omega), tolerance = 1e-12)
    }
  }
})

# ===========================================================================
# 3. The fix at the level that actually matters: the active-set EIF path
# ===========================================================================

test_that("the active-set EIF path now uses canonical reverse_kl h' end to end", {
  # compute_omega_reduced() exists to serve calculate_eif_lambda_active(). Verify
  # the corrected h' propagates all the way to EIF values, and that a reverse_kl
  # active-set run differs from what the buggy Omega would have produced.
  withr::local_seed(20260811L)
  m <- 3L
  n <- 25L
  a_mat <- generate_binary_vectors(m)
  a_ext <- cbind(1, a_mat)
  n_comb <- nrow(a_mat)

  ps <- matrix(stats::runif(n * n_comb, 0.3, 1), nrow = n)
  ps <- ps / rowSums(ps)
  mu <- matrix(stats::rnorm(n * n_comb, mean = 1, sd = 0.3), nrow = n)
  obs_idx <- apply(ps, 1, function(pr) sample.int(n_comb, 1, prob = pr))
  A <- a_mat[obs_idx, , drop = FALSE]
  Y <- mu[cbind(seq_len(n), obs_idx)] + stats::rnorm(n, sd = 0.1)
  delta <- rep(1.2, m)

  # `active_sets` is a required argument; build it the way the package's own
  # active-set entry points do. A small epsilon keeps every treatment active, so
  # the reduced system is the full (m+1) one and Omega is directly comparable to
  # the canonical construction below.
  active_sets <- compute_active_sets(extract_marginal_probabilities(ps, m = m),
                                     epsilon = 0.05)
  expect_true(all(active_sets$sufficient_active))

  res <- suppressWarnings(efficient_influence_function_active(
    A = A, Y = Y, propensity_scores = ps, delta = delta, outcome_means = mu,
    divergence = "reverse_kl", active_sets = active_sets,
    feasibility_check = FALSE))

  expect_length(res$eif_values, n)
  expect_true(all(is.finite(res$eif_values)))

  # Per-observation cross-check: the Omega the active path builds must match the
  # canonical one for every observation.
  qd <- lapply(seq_len(n), function(i) {
    compute_qdelta_per_obs(ps[i, ], delta, "reverse_kl", tol = 1e-10,
                           combination_matrix = a_mat)
  })
  for (i in seq_len(n)) {
    lambda <- qd[[i]]$lambda
    linear_terms <- as.vector(a_ext %*% lambda)
    canonical <- orkl_omega_from_hprime(
      ps[i, ], a_ext, divergence_h_pair(linear_terms, "reverse_kl")$h_prime)
    actual <- suppressWarnings(
      compute_omega_reduced(ps[i, ], lambda, a_ext, "reverse_kl"))
    expect_equal(actual, canonical, tolerance = 1e-10,
                 info = sprintf("observation %d", i))
  }
})

test_that("existing omega diagnostic equivalence (KL) is preserved", {
  # test-omega-diagnostic.R asserts compute_omega_reduced() reproduces
  # calculate_eif_lambda()'s Omega for KL. Re-assert it here against the
  # post-stage-2 (crossprod) form of calculate_eif_lambda()'s assembly, so this
  # file also guards the cross-path agreement the original bug hunt started from.
  withr::local_seed(20260811L)
  m <- 3L
  a_mat <- generate_binary_vectors(m)
  a_ext <- cbind(1, a_mat)
  ps <- stats::runif(nrow(a_mat), 0.3, 1)
  ps <- ps / sum(ps)
  qd <- compute_qdelta_per_obs(ps, rep(1.5, m), "kl", tol = 1e-8,
                               combination_matrix = a_mat)

  linear_terms <- as.vector(a_ext %*% qd$lambda)
  h_prime <- divergence_h_pair(linear_terms, "kl")$h_prime
  eif_form <- crossprod(a_ext * as.vector(ps * h_prime), a_ext)   # calculate_eif_lambda()'s form

  expect_equal(compute_omega_reduced(ps, qd$lambda, a_ext, "kl"), eif_form,
               tolerance = 1e-10)
})
