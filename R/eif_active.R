#' Compute EIF Lambda with Pre-Specified Active Set
#'
#' Modified version of calculate_eif_lambda that uses a pre-specified active set
#' instead of computing it on the fly. Used for threshold-based active set selection.
#'
#' @param A Numeric vector (length m) of observed treatment for this observation.
#' @param propensity_scores Numeric vector (length 2^m) of propensity scores for this observation.
#' @param qdelta_result List with lambda, q_delta, delta from Q_delta solver.
#' @param h_values Numeric vector (length 2^m) of h(a; lambda) values.
#' @param h_prime Numeric vector (length 2^m) of h'(a; lambda) values.
#' @param a_mat Matrix (2^m x (m+1)) of treatment combinations with intercept.
#' @param divergence Character, divergence type.
#' @param active_indices Integer vector of active treatment indices (1 to m).
#'   Empty vector means no active treatments (use fallback).
#'
#' @return Numeric vector of EIF_lambda values (length m+1 for active constraints).
#'
#' @details
#' This function differs from the original calculate_eif_lambda in that it:
#' - Takes active set as input rather than computing it
#' - Uses compute_omega_reduced() from treatment_subset_utils.R
#' - Handles empty active sets (returns zero EIF_lambda)
#'
#' If active_indices is empty, returns zero vector (no active constraints,
#' so lambda derivative is zero).
#'
#' @keywords internal
calculate_eif_lambda_active <- function(A, propensity_scores, qdelta_result,
                                        h_values, h_prime, a_mat, divergence,
                                        active_indices) {
  m <- length(A)

  # Extract Q_delta results
  lambda <- qdelta_result$lambda
  q_delta <- qdelta_result$q_delta
  delta <- qdelta_result$delta

  # Handle empty active set: no constraints, so EIF_lambda = 0
  if (length(active_indices) == 0) {
    return(numeric(m + 1))
  }

  # Active constraints include intercept + active treatments
  active_constraints <- c(1, active_indices + 1)
  n_active <- length(active_constraints)

  # Extract reduced lambda for active constraints only
  lambda_active <- lambda[active_constraints]

  # Extract active treatment combinations from full a_mat
  a_mat_active <- a_mat[, active_constraints, drop = FALSE]

  # Compute Ω for active set using shared utility
  # Note: compute_omega_reduced expects REDUCED lambda (length = |active|+1)
  omega_active <- compute_omega_reduced(
    propensity_scores = propensity_scores,
    lambda = lambda_active,
    a_bar_matrix = a_mat_active,
    divergence = divergence
  )

  # Invert Ω (with regularization if needed)
  regularization_applied <- FALSE
  tryCatch({
    omega_inv <- solve(omega_active)
    # Check for Inf/NaN in result
    if (any(is.infinite(omega_inv)) || any(is.na(omega_inv))) {
      stop("Matrix inverse contains Inf/NaN")
    }
  }, error = function(e) {
    regularization_applied <<- TRUE
    omega_inv <<- solve(omega_active + OMEGA_REGULARIZATION_EPS * diag(n_active))
  })

  if (regularization_applied) {
    warning(
      "Omega matrix regularization applied in EIF lambda calculation with active sets. ",
      "This may indicate extreme propensity scores or numerical instability.",
      call. = FALSE
    )
  }

  # Compute g(bar_A; lambda) - same as original
  a_prime_extended <- c(1, A)

  # Find h value for observed pattern
  pattern_match <- apply(a_mat, 1, function(row) all(row == a_prime_extended))
  pattern_idx <- which(pattern_match)
  if (length(pattern_idx) != 1) {
    stop("Could not uniquely identify observed treatment pattern in calculate_eif_lambda_active")
  }
  h_A <- h_values[pattern_idx]

  # E_Q(bar_A | X) - extended marginals under Q_delta
  E_QA <- as.numeric(q_delta %*% a_mat)

  # E(bar_A | X) - extended marginals under P_0
  # Extract marginals for all treatments (not just active)
  marginal_probs <- extract_marginal_probabilities(
    matrix(propensity_scores, nrow = 1),
    m = m
  )[1, ]
  E_A <- c(1, marginal_probs)

  # Diagonal of \bar{Delta}
  delta_extended <- c(1, delta)
  Delta_diag <- delta_extended / (1 + E_A * (delta_extended - 1))
  Delta2_diag <- Delta_diag^2
  Delta2_over_delta <- Delta2_diag / delta_extended

  # g(bar_A; lambda)
  g <- h_A * a_prime_extended - E_QA - Delta2_over_delta * (a_prime_extended - E_A)

  # EIF_lambda = Omega^{-1} * g (for active constraints only)
  eif_lambda_active <- as.numeric(omega_inv %*% g[active_constraints])

  # Expand to full dimension (m+1) with zeros for inactive constraints
  eif_lambda_full <- numeric(m + 1)
  eif_lambda_full[active_constraints] <- eif_lambda_active

  eif_lambda_full
}


#' Efficient Influence Function with Active Set Selection
#'
#' Computes EIF for vector incremental effects with threshold-based active set
#' selection. Uses reduced-dimension Ω based on active constraints per observation.
#'
#' @inheritParams efficient_influence_function
#' @param active_sets List returned by \code{\link{compute_active_sets}} containing:
#'   \itemize{
#'     \item active_matrix: Logical matrix (n x m) indicating active treatments
#'     \item active_indices: List of length n with active treatment indices per obs
#'     \item n_active: Integer vector with number of active treatments per obs
#'     \item sufficient_active: Logical vector indicating sufficient active treatments
#'   }
#'
#' @return Same structure as \code{\link{efficient_influence_function}}, with
#'   additional components:
#'   \item{active_sets}{The input active_sets list}
#'   \item{n_insufficient_active}{Number of observations with insufficient active treatments}
#'
#' @details
#' **Difference from standard EIF:**
#' - Uses pre-computed active sets `K_ε(X) = {k: m_k(X) ∈ [ε, 1-ε]}`
#' - Computes reduced-dimension Ω based on active constraints only
#' - Observations with insufficient active treatments handled specially
#'
#' **Active set motivation:**
#' Treatments with marginal propensities near 0 or 1 provide little information
#' and can make Ω singular. Restricting to active set:
#' - Improves numerical stability
#' - Reduces effective dimension
#' - Satisfies margin conditions for asymptotic theory
#'
#' **Insufficient active treatments:**
#' If `active_sets$sufficient_active[i]` is `FALSE`, the EIF for that observation
#' is set to zero (no contribution). The user should verify this is acceptable
#' or adjust the threshold ε.
#'
#' @seealso
#' \code{\link{estimate_qdelta_mean_stratified}} for eligibility-based stratification
#' (complementary approach using eligibility indicators rather than propensity thresholds)
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 3, p = 2, seed = 123)
#'
#' # Estimate nuisances
#' ps <- estimate_propensity_scores(sim_data$data, method = "glm")
#' outcomes <- estimate_outcome_models(sim_data$data, method = "glm")
#'
#' # Extract marginal propensities
#' A_mat <- as.matrix(sim_data$data[, c("A1", "A2", "A3")])
#' marginal_probs <- extract_marginal_probabilities(ps, m = 3)
#'
#' # Compute vanishing threshold
#' eps_n <- compute_vanishing_threshold(n = 500, gamma = 0.3)
#'
#' # Identify active sets
#' active <- compute_active_sets(marginal_probs, epsilon = eps_n)
#'
#' # Compute EIF with active sets
#' eif_result <- efficient_influence_function_active(
#'   A = A_mat,
#'   Y = sim_data$data$Y,
#'   propensity_scores = ps,
#'   outcome_means = outcomes,
#'   delta = c(0.1, 0.1, 0.1),
#'   divergence = "kl",
#'   active_sets = active
#' )
#'
#' # Check how many observations had insufficient active treatments
#' print(eif_result$n_insufficient_active)
#' }
efficient_influence_function_active <- function(X = NULL, A, Y = NULL, R = NULL, pR = NULL,
                                                propensity_scores,
                                                delta,
                                                outcome_means = NULL,
                                                divergence = "kl",
                                                qdelta_results = NULL,
                                                active_sets,
                                                feasibility_check = FALSE) {

  # Determine dimensions
  if (!is.null(X)) {
    n <- nrow(X)
  } else if (!is.null(A)) {
    n <- nrow(A)
  } else {
    n <- nrow(propensity_scores)
  }

  m <- ncol(A)
  n_combinations <- 2^m

  # Validate inputs (same as standard EIF)
  if (sum(c(is.null(R), is.null(pR))) == 1) {
    stop("Either both R and pR must be null or non-null")
  }
  if (is.null(outcome_means)) stop("Outcome means must be provided")
  if (is.null(Y)) stop("Y must be provided")

  if (any(!is.finite(propensity_scores))) {
    stop("Propensity scores contain non-finite values")
  }

  ps_row_sums <- rowSums(propensity_scores)
  if (!all(abs(ps_row_sums - 1) < PROPENSITY_SUM_TOL)) {
    max_deviation <- max(abs(ps_row_sums - 1))
    stop("Propensity scores do not sum to 1 (max deviation: ", sprintf("%.2e", max_deviation), ")")
  }

  if (any(!is.finite(outcome_means))) {
    stop("Outcome means contain non-finite values")
  }

  if (any(!is.finite(Y))) {
    stop("Outcome Y contains non-finite values")
  }

  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a positive numeric vector")
  }

  # Validate active_sets structure
  if (!is.list(active_sets) ||
      !all(c("active_matrix", "active_indices", "n_active", "sufficient_active") %in% names(active_sets))) {
    stop("active_sets must be a list from compute_active_sets() with required components")
  }

  if (nrow(active_sets$active_matrix) != n || ncol(active_sets$active_matrix) != m) {
    stop("active_sets dimensions do not match data (expected ", n, " x ", m, ")")
  }

  # Count observations with insufficient active treatments
  n_insufficient <- sum(!active_sets$sufficient_active)
  if (n_insufficient > 0) {
    warning(sprintf(
      "%d out of %d observations (%.1f%%) have insufficient active treatments.\nTheir EIF contributions will be zero. Consider adjusting threshold epsilon.",
      n_insufficient, n, 100 * n_insufficient / n
    ), call. = FALSE)
  }

  # Feasibility check if requested (with early stopping for fast failure)
  if (feasibility_check) {
    feasibility_analysis <- analyze_delta_feasibility_all_obs(propensity_scores, delta,
                                                              early_stop = TRUE)
    if (feasibility_analysis$n_feasible < n) {
      # Check if early stopped
      if (!is.null(feasibility_analysis$early_stopped) && feasibility_analysis$early_stopped) {
        stop(sprintf(
          "Delta infeasible for observation %d (first infeasible found).\nOnly checked %d / %d observations before stopping early.\nConsider using a smaller delta (all components < %.2f).",
          feasibility_analysis$first_infeasible,
          feasibility_analysis$first_infeasible,
          n,
          min(delta) * 0.9
        ), call. = FALSE)
      }

      # Full check: try to suggest alternative
      suggested_result <- suggest_feasible_delta(propensity_scores, delta)
      if (suggested_result$n_obs_feasible == n) {
        stop(paste("Delta infeasible for", n - feasibility_analysis$n_feasible, "observations.",
                   "\nSuggested delta:", paste(round(suggested_result$suggested_delta, 4), collapse = ", ")))
      } else {
        stop(paste("Delta infeasible for", n - feasibility_analysis$n_feasible, "observations.",
                   "\nNo fully feasible alternative found. Reduce delta or check propensity scores."))
      }
    }
  }

  # Prepare treatment combinations
  a_mat <- generate_binary_vectors(m)
  a_extended_matrix <- cbind(1, a_mat)

  # Convert to matrices
  propensity_scores <- as.matrix(propensity_scores)
  outcome_means <- as.matrix(outcome_means)
  A <- as.matrix(A)
  Y <- as.numeric(Y)

  # Compute or use provided Q_delta results
  if (is.null(qdelta_results)) {
    qdelta_results <- vector("list", n)
    for (i in 1:n) {
      # a_mat depends only on m: pass the matrix built above rather than letting
      # each per-observation solve regenerate it.
      qdelta_results[[i]] <- compute_qdelta_per_obs(
        propensity_scores[i, ],
        delta,
        divergence = divergence,
        combination_matrix = a_mat
      )
    }
  }

  # Extract lambda and compute H matrices
  LAMBDA <- as.matrix(do.call(rbind, lapply(qdelta_results, function(x) x$lambda)))
  H_result <- compute_H_vectorized_general(LAMBDA, a_extended_matrix, divergence)
  H_mat <- H_result$H_mat
  H_prime <- H_result$H_prime

  # Compute indicator matrix (which treatment combination was observed)
  A_patterns <- as.vector(A %*% (2^(0:(m-1))))
  a_patterns <- as.vector(a_mat %*% (2^(0:(m-1))))
  I_mat <- outer(A_patterns, a_patterns, "==") * 1.0

  # Compute EIF_lambda for each observation (using active sets)
  EIF_LAMBDA <- list()
  for (i in 1:n) {
    if (active_sets$sufficient_active[i]) {
      # Use active set version
      EIF_LAMBDA[[i]] <- calculate_eif_lambda_active(
        A[i, ],
        propensity_scores[i, ],
        qdelta_results[[i]],
        H_mat[i, ],
        t(I_mat[i, ]) %*% H_prime[i, ],
        a_extended_matrix,
        divergence,
        active_indices = active_sets$active_indices[[i]]
      )
    } else {
      # Insufficient active treatments: EIF_lambda = 0
      EIF_LAMBDA[[i]] <- numeric(m + 1)
    }
  }
  EIF_LAMBDA <- as.matrix(do.call(rbind, EIF_LAMBDA))

  # Extract Q_delta values
  Q_delta <- as.matrix(do.call(rbind, lapply(qdelta_results, function(x) x$q_delta)))

  # h_mat is the same h(u) already computed above as H_mat via
  # compute_H_vectorized_general() (single canonical source of h/h', see
  # divergence_h_functions.R). Previously this was recomputed a second time
  # here with its own (differently-clipped) formula -- redundant work, and
  # a source of eif1-vs-eif2/eif3 inconsistency near domain boundaries for
  # reverse_kl/hellinger. See [LEARN 2026-08-10] in divergence_h_functions.R.
  h_mat <- H_mat

  # Expand Y to matrix for broadcasting
  Y_mat <- matrix(Y, nrow = n, ncol = n_combinations)

  # Compute Phi component
  if (is.null(R) && is.null(pR)) {
    # No missingness
    Phi <- h_mat * (I_mat * (Y_mat - outcome_means) +
                    propensity_scores * outcome_means)
  } else {
    # With missingness
    safe_pR <- pmax(pR, 1e-100)
    Phi <- h_mat * (R * I_mat * (Y_mat - outcome_means) / safe_pR +
                    propensity_scores * outcome_means)
  }

  # Compute residual component
  HResid <- (I_mat - propensity_scores) * H_mat

  # Compute EIF components
  eif1 <- rowSums(Phi)
  eif2 <- rowSums(outcome_means * HResid)
  eif3 <- rowSums(outcome_means * propensity_scores * H_prime *
                  (EIF_LAMBDA %*% t(a_extended_matrix)))

  plugin <- rowSums(outcome_means * Q_delta)
  eif_values <- eif1 + eif2 - eif3

  # Observations with insufficient active treatments (empty active set): ALL
  # coordinates are near-boundary, so the intervention is the no-op (delta -> 1),
  # Q_delta = P0, and the target contributes E[Y^{P0}] = E[Y]. The efficient
  # influence function for E[Y] is simply Y_i (the AIPW/one-step contribution: in
  # the no-tilt limit h=1, H=1, EIF_lambda=0, so eif1+eif2-eif3 collapses to Y_i).
  # PREVIOUS BEHAVIOR set these to 0, which biases the estimate of E[Y^Q] toward
  # 0 in proportion to the empty-active fraction (the near-boundary regime). Using
  # Y_i is the correct P0 plug-in / one-step contribution and keeps the aggregate
  # a valid estimator of the mixed estimand (active obs: psi(delta_K); all-inactive
  # obs: E[Y|X] integrated under P0). See plan 2026-07-10_adaptive-empty-active-set-fix.
  insufficient <- !active_sets$sufficient_active
  if (any(insufficient)) {
    eif_values[insufficient] <- Y[insufficient]
    plugin[insufficient]     <- rowSums(outcome_means[insufficient, , drop = FALSE] *
                                        propensity_scores[insufficient, , drop = FALSE])
  }

  # Convergence diagnostics
  converged <- sapply(qdelta_results, function(x) x$converged)
  iterations <- sapply(qdelta_results, function(x) x$iterations)
  constraint_norms <- sapply(qdelta_results, function(x) x$constraint_norm)

  convergence_info <- list(
    converged = converged,
    iterations = iterations,
    constraint_norms = constraint_norms
  )

  convergence_summary <- list(
    n_converged = sum(converged),
    convergence_rate = mean(converged),
    mean_iterations = mean(iterations),
    max_iterations = max(iterations),
    mean_constraint_norm = mean(constraint_norms),
    max_constraint_norm = max(constraint_norms)
  )

  # Warning if convergence rate low
  if (convergence_summary$convergence_rate < MIN_CONVERGENCE_RATE) {
    warning(sprintf(
      "Low Q_delta convergence rate: %.1f%% (%d/%d converged).\n",
      100 * convergence_summary$convergence_rate, convergence_summary$n_converged, n,
      "Check propensity scores and delta feasibility."
    ))
  }

  return(list(
    eif_values = eif_values,
    centered_eif = eif_values - mean(eif_values),
    plugin = plugin,
    qdelta_results = qdelta_results,
    convergence_info = convergence_info,
    convergence_summary = convergence_summary,
    active_sets = active_sets,
    n_insufficient_active = n_insufficient
  ))
}
