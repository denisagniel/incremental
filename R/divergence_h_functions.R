#' @file divergence_h_functions.R
#' @title Canonical Divergence h-Function and Derivative
#' @description
#' Single source of truth for the divergence conjugate function h(x) and its
#' derivative h'(x), shared by the Newton solver (`fast_qdelta_solver.R`), the
#' EIF computation (`eif.R`, `eif_active.R`), and the estimation helper
#' (`estimation.R`).
#'
#' `[LEARN 2026-08-10]` Before this file existed, h(x)/h'(x) were reimplemented
#' independently at ~6 call sites across 4 files, with 4-5 different
#' boundary-clipping conventions for reverse_kl and hellinger (chi2 and kl
#' were unaffected -- chi2 has no singularity, kl's only guard is overflow
#' protection). Two conventions disagreed by 18 orders of magnitude on the
#' override value at the exact same domain-violation condition (1e-6 vs.
#' 1e12), and within a single EIF computation, eif1 used one convention
#' while eif2/eif3 used another. The solver's analytical Jacobian
#' additionally evaluated h' at the UNCLIPPED linear term while the moment
#' sum evaluated h at the CLIPPED linear term -- so J != dF/dlambda whenever
#' clipping was active, a genuine (if rare -- clipping only engages outside
#' the feasible region enforced by check_constraint_satisfaction_*) Newton
#' defect. See session_notes/2026-08-10.md and MEMORY.md.
#'
#' Fix: one function, called everywhere. Clip the linear term ONCE, then
#' evaluate both h and h' at the SAME clipped value. The chosen constants are
#' exactly the pre-existing solver-weights convention
#' (`compute_q_delta_weights_internal()`), re-expressed as a symmetric clip so
#' the old two-line "clip + override near boundary" pattern collapses to one
#' clip + one formula, with no numerical change in any regime the solver's
#' weights function already covered (verified: for reverse_kl, the old
#' override value 1e12 is exactly `-1 / (-REVERSE_KL_DOMAIN_BUFFER)`, i.e.
#' what the clipped formula evaluates to at the boundary, by construction of
#' the constants in constants.R -- the now-deleted H_VALUE_MAX/H_VALUE_SMALL
#' constants were the old, inconsistent override values).
#' @author Denis Agniel

#' Divergence h(x) and h'(x), evaluated at a single clipped linear term
#'
#' @description
#' Computes h(x) and h'(x) for the specified divergence's constrained-tilt
#' problem, where x is the linear term `lambda[1] + a' lambda[2:(m+1)]`. x is
#' clipped to the interior of h's domain before EITHER h or h' is evaluated,
#' so both are guaranteed to be functions of the identical (possibly clipped)
#' argument -- required for `h'` to be the correct Jacobian entry for the
#' residual actually computed from `h`.
#'
#' @param x Numeric vector or matrix of linear terms. Any shape; computation
#'   is elementwise, so this works for a single scalar (per-observation
#'   solver loops), a vector (one observation, all `2^m` combinations), or a
#'   matrix (many observations x combinations, e.g. `eif.R`'s `LAMBDA %*%
#'   t(a_extended_matrix)`).
#' @param divergence Character, one of "kl", "reverse_kl", "chi2", "hellinger".
#' @return List with:
#'   \item{h}{h(x) at the clipped x, same shape as `x`}
#'   \item{h_prime}{h'(x) at the SAME clipped x, same shape as `x`}
#'   \item{n_clipped}{integer count of elements where clipping was active
#'     (diagnostic; nonzero means an iterate came within the numerical
#'     safety buffer of the domain boundary -- for a converged solution this
#'     should essentially never happen, since
#'     `check_constraint_satisfaction_*()` rejects iterates well before this
#'     buffer is reached)}
#' @keywords internal
divergence_h_pair <- function(x, divergence) {
  if (divergence == "kl") {
    # h(x) = exp(x), h'(x) = exp(x). Clip only guards exp() overflow.
    xs <- pmin(x, KL_LINEAR_TERM_MAX)
    h <- exp(xs)
    h_prime <- h
    n_clipped <- sum(xs != x)
  } else if (divergence == "reverse_kl") {
    # h(x) = -1/x, h'(x) = 1/x^2. Domain: x < 0. Symmetric clip: floor
    # guards the (benign, h -> 0) end; ceiling guards the (singular, h ->
    # +Inf) end and reproduces the old override exactly (-1/(-buf) = 1/buf).
    xs <- pmin(pmax(x, REVERSE_KL_LINEAR_TERM_MIN), -REVERSE_KL_DOMAIN_BUFFER)
    h <- -1 / xs
    h_prime <- 1 / xs^2
    n_clipped <- sum(xs != x)
  } else if (divergence == "chi2") {
    # h(x) = x/2 + 1, h'(x) = 1/2 (constant). No singularity: h -> 0 at the
    # domain boundary x = -2 is the correct, finite limit; no clip needed.
    h <- x / 2 + 1
    h_prime <- if (is.matrix(x)) matrix(0.5, nrow(x), ncol(x)) else rep(0.5, length(x))
    n_clipped <- 0L
  } else if (divergence == "hellinger") {
    # h(x) = 1/(1-x)^2, h'(x) = 2/(1-x)^3. Domain: x < 1. Clip the DISTANCE
    # to the boundary (d = 1-x), not x itself: clipping x near 1 and then
    # computing 1-x loses precision to cancellation, whereas clipping d
    # directly is exact.
    d <- pmax(1 - x, HELLINGER_DOMAIN_BUFFER)
    h <- 1 / d^2
    h_prime <- 2 / d^3
    n_clipped <- sum(d != (1 - x))
  } else {
    stop("divergence must be one of 'kl', 'reverse_kl', 'chi2', 'hellinger', not '",
         divergence, "'", call. = FALSE)
  }

  list(h = h, h_prime = h_prime, n_clipped = n_clipped)
}
