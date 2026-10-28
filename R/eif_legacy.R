#' @file eif_legacy.R
#' @title Frozen pre-optimization EIF pass (equivalence reference ONLY)
#' @description
#' Verbatim snapshot of `calculate_eif_lambda()` and `eif_compute_chunk()` as they
#' stood immediately before the 2026-08-11 stage-2 Omega-construction fix, kept in
#' the tree so old-vs-new comparison never depends on a git checkout. This mirrors
#' the discipline established for the solver in `R/fast_qdelta_solver_legacy.R`
#' (spec MUST item "Numerical equivalence, old vs. new").
#'
#' What stage 2 changed in `R/eif.R`, and therefore what this file preserves:
#'   1. **Omega construction.** Was
#'      `weight_mat <- diag(weights); omega_matrix <- t(a_mat) %*% weight_mat %*% a_mat`
#'      -- materializing a dense `2^m x 2^m` diagonal matrix (33.6 MB at m = 11)
#'      per observation and multiplying through it (~50M flops instead of ~0.3M).
#'      Now `crossprod(a_mat * weights, a_mat)`, the exact algebraic identity.
#'      Profiled at 95.8% of the whole EIF pass at m = 11 before the fix.
#'   2. **Observed-pattern lookup.** Was
#'      `apply(a_mat, 1, function(row) all(row == a_prime_extended))` -- `2^m`
#'      interpreted closure calls per observation. Now an exact integer bit-key
#'      prefilter (one BLAS dgemv) plus elementwise verification of the surviving
#'      candidates, so no-match and duplicate-match still reach the same `stop()`.
#'   3. **Marginals.** Was `compute_marginals(propensity_scores, a_mat[, -1, drop = FALSE], ...)`
#'      -- copying a `2^m x m` matrix per observation (184 KB at m = 11). Now
#'      `compute_marginals(propensity_scores, a_mat, ...)[-1]`, dropping the
#'      intercept from the result instead of slicing the input (verified bitwise
#'      identical for m in 1..11).
#' Each of (2) and (3) only became visible as a dominant cost after the preceding
#' fix removed the term above it; all three are exact, not approximations.
#'
#' Nothing else differs: the `tryCatch`/`solve()` regularization fallback, the
#' active-constraint selection, the `g` vector, the zero-padding and every returned
#' value are untouched.

#'
#' Both functions here are suffixed `_legacy`, are called by nothing in the package
#' except `tests/testthat/test-eif-omega-vectorized-equivalence.R`, and are NOT
#' exported, NOT a fallback and NOT maintained. Do not "fix" anything in this file.
#' If it needs to change, the optimization is wrong.
#' @author Denis Agniel
#' @keywords internal
#' @name eif_legacy
NULL

#' Legacy (dense-diag) per-observation EIF-lambda calculation
#'
#' @inheritParams calculate_eif_lambda
#' @return Numeric vector of length m+1
#' @keywords internal
calculate_eif_lambda_legacy <- function(A, propensity_scores, qdelta_result, h_values, h_dot_values, a_mat, divergence) {
  # Get treatment assignment for this observation
  #treatment_combinations <- generate_treatment_combinations(m = length(A))

  # Get Q_delta result for this observation
  lambda <- qdelta_result$lambda
  q_delta <- qdelta_result$q_delta
  delta <- qdelta_result$delta

  # Determine active constraints for this observation.
  # drop = FALSE keeps the treatment design a matrix when m = 1 (single eligible
  # treatment), where a_mat[, -1] would otherwise collapse to a vector and break
  # compute_marginals (nrow(vectors) becomes NULL).
  E_A <- compute_marginals(propensity_scores, a_mat[, -1, drop = FALSE], m = length(A))
  active_constraints <- find_active_constraints(E_A, tol = 1e-8)
  active_constraints <- c(1, (active_constraints + 1))

  # Calculate the Omega matrix (e.g. E[h(A)h(A)'] for KL)
  #a_extended_matrix <- cbind(1, treatment_combinations)
  n_combinations <- length(propensity_scores)
  weights <- as.numeric(propensity_scores * h_dot_values)
  # Optimized: use diag() directly instead of creating zero matrix
  weight_mat <- diag(weights)
  omega_matrix <- (t(a_mat) %*% weight_mat %*% a_mat)
  omega_active <- omega_matrix[active_constraints, active_constraints]

  # Number of active constraints (including intercept)
  n_active <- length(active_constraints)

  # Compute inverse with regularization if needed; is there a faster function than solve?
  regularization_applied <- FALSE
  tryCatch({
    omega_inv <- solve(omega_active)
    # Check for Inf/NaN in result (can happen with near-singular matrices)
    if (any(is.infinite(omega_inv)) || any(is.na(omega_inv))) {
      stop("Matrix inverse contains Inf/NaN")
    }
  }, error = function(e) {
    # Add stronger regularization for near-singular matrices
    regularization_applied <<- TRUE
    omega_inv <<- solve(omega_active + OMEGA_REGULARIZATION_EPS * diag(n_active))
  })

  if (regularization_applied) {
    warning(
      "Omega matrix regularization applied in EIF lambda calculation. ",
      "This indicates the Omega matrix is near-singular, which may occur with: ",
      "(1) extreme propensity scores, (2) redundant constraints, (3) numerical instability. ",
      "Results may be less reliable.",
      call. = FALSE
    )
  }

  # Compute g(bar_A; lambda) from paper equation (199)
  # g = h(bar_A) * bar_A - E_Q(bar_A|X) - Delta^2 * (bar_A - E(bar_A|X))

  # bar_A = (1, A) - observed treatment extended
  a_prime_extended <- c(1, A)

  # Find h value for observed pattern bar_A
  pattern_match <- apply(a_mat, 1, function(row) all(row == a_prime_extended))
  pattern_idx <- which(pattern_match)
  if (length(pattern_idx) != 1) {
    stop("Could not uniquely identify observed treatment pattern in calculate_eif_lambda")
  }
  h_A <- h_values[pattern_idx]  # h(bar_A; lambda) - scalar

  # E_Q(bar_A | X) - extended marginals under Q_delta
  E_QA <- as.numeric(q_delta %*% a_mat)

  # E(bar_A | X) - extended marginals under P_0
  E_A <- c(1, E_A)

  # Diagonal of \bar{Delta} (definition: Delta_{kk} = delta_k / (1 + P_0(A_k=1|X)(delta_k - 1)))
  delta_extended <- c(1, delta)  # (1, δ_1, ..., δ_m)
  Delta_diag <- delta_extended / (1 + E_A * (delta_extended - 1))
  Delta2_diag <- Delta_diag^2

  # From implicit function theorem derivation:
  # ∂ζ_k/∂π_a = (Δ_{kk}^2 / δ_k) · a_k
  # So g should have term: Δ^2/δ * (bar_A - E(bar_A))
  Delta2_over_delta <- Delta2_diag / delta_extended  # Δ^2 / δ

  # g(bar_A; lambda) from my derivation
  g <- h_A * a_prime_extended - E_QA - Delta2_over_delta * (a_prime_extended - E_A)

  # EIF_lambda = Omega^{-1} * g (for active constraints)
  # Note: minus sign applied in EIF computation (line 338: eif_values = eif1 + eif2 - eif3)
  eif_active <- as.numeric(omega_inv %*% g[active_constraints])

  # Pad to full dimension (m+1) with zeros for inactive constraints
  # This ensures consistent dimensions when rbind-ing across observations
  m_plus_1 <- ncol(a_mat)
  if (length(active_constraints) < m_plus_1) {
    eif_full <- rep(0, m_plus_1)
    eif_full[active_constraints] <- eif_active
    return(eif_full)
  }
  eif_active
}

#' Legacy EIF chunk driver (calls the legacy per-observation EIF-lambda)
#'
#' Identical to `eif_compute_chunk()` except that it calls
#' `calculate_eif_lambda_legacy()`. Exists so the chunk-level A/B (profiling and
#' benchmarking) runs both implementations in one process with no session-level
#' confounds.
#'
#' @inheritParams eif_compute_chunk
#' @return List with eif_values and plugin for the chunk
#' @keywords internal
eif_compute_chunk_legacy <- function(chunk_indices, A, Y, R, pR, propensity_scores, outcome_means,
                                     qdelta_results, a_mat, a_extended_matrix, divergence, m) {

  n_chunk <- length(chunk_indices)
  n_combinations <- 2^m

  # Extract chunk data
  A_chunk <- A[chunk_indices, , drop = FALSE]
  Y_chunk <- Y[chunk_indices]
  propensity_scores_chunk <- propensity_scores[chunk_indices, , drop = FALSE]
  outcome_means_chunk <- outcome_means[chunk_indices, , drop = FALSE]
  qdelta_results_chunk <- qdelta_results[chunk_indices]

  if (!is.null(R)) {
    R_chunk <- R[chunk_indices]
    pR_chunk <- pR[chunk_indices, , drop = FALSE]
  } else {
    R_chunk <- NULL
    pR_chunk <- NULL
  }

  # Compute indicator matrix for chunk using optimized bit-pattern matching
  A_patterns <- as.vector(A_chunk %*% (2^(0:(m-1))))
  a_patterns <- as.vector(a_mat %*% (2^(0:(m-1))))
  I_chunk <- outer(A_patterns, a_patterns, "==") * 1.0

  # Extract lambda and compute H matrices for chunk
  LAMBDA_chunk <- as.matrix(do.call(rbind, lapply(qdelta_results_chunk, function(x) x$lambda)))
  H_result <- compute_H_vectorized_general(LAMBDA_chunk, a_extended_matrix, divergence)
  H_mat_chunk <- H_result$H_mat
  H_prime_chunk <- H_result$H_prime

  # Compute EIF_lambda for each observation in chunk
  EIF_LAMBDA_chunk <- list()
  for (i in 1:n_chunk) {
    EIF_LAMBDA_chunk[[i]] <- calculate_eif_lambda_legacy(
      A_chunk[i, ], propensity_scores_chunk[i, ],
      qdelta_results_chunk[[i]],
      H_mat_chunk[i, ], H_prime_chunk[i, ],
      a_extended_matrix, divergence
    )
  }
  EIF_LAMBDA_chunk <- as.matrix(do.call(rbind, EIF_LAMBDA_chunk))

  # Extract Q_delta values for chunk
  Q_chunk <- as.matrix(do.call(rbind, lapply(qdelta_results_chunk, function(x) x$q_delta)))

  h_chunk <- H_mat_chunk

  # Expand Y_chunk to matrix for broadcasting (n_chunk x n_combinations)
  Y_mat <- matrix(Y_chunk, nrow = n_chunk, ncol = n_combinations)

  # Compute Phi component using h instead of dividing by pi
  if (is.null(R_chunk) && is.null(pR_chunk)) {
    Phi_chunk <- h_chunk * (I_chunk * (Y_mat - outcome_means_chunk) +
                            propensity_scores_chunk * outcome_means_chunk)
  } else {
    safe_pR <- pmax(pR_chunk, 1e-100)
    Phi_chunk <- h_chunk * (R_chunk * I_chunk * (Y_mat - outcome_means_chunk) / safe_pR +
                            propensity_scores_chunk * outcome_means_chunk)
  }

  # Compute residual component
  HResid_chunk <- (I_chunk - propensity_scores_chunk) * H_mat_chunk

  # Compute EIF components
  eif1_chunk <- rowSums(Phi_chunk)
  eif2_chunk <- rowSums(outcome_means_chunk * HResid_chunk)
  eif3_chunk <- rowSums(outcome_means_chunk * propensity_scores_chunk * H_prime_chunk *
                        (EIF_LAMBDA_chunk %*% t(a_extended_matrix)))

  plugin_chunk <- rowSums(outcome_means_chunk * Q_chunk)
  eif_values_chunk <- eif1_chunk + eif2_chunk - eif3_chunk

  return(list(
    eif_values = eif_values_chunk,
    plugin = plugin_chunk
  ))
}
