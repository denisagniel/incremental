#' @file adaptive_delta.R
#' @title Adaptive Delta Computation for Fixed-Threshold Estimator
#' @description
#' Functions for computing observation-specific delta_k(X) values for the
#' adaptive estimator (fixed epsilon, adaptive delta). Uses localized margin
#' conditions to adapt delta based on propensity patterns.

#' Compute Adaptive Delta for Single Observation
#'
#' Computes observation-specific delta_k(X) based on marginal propensities
#' and active set membership. For the adaptive estimator with fixed epsilon.
#'
#' @param marginal_probs Numeric vector (length m) of marginal propensities for this observation
#' @param delta_target Numeric vector (length m) of target deltas (user-requested)
#' @param active_indices Integer vector of active treatment indices (from active set)
#' @param epsilon Numeric, the fixed threshold for active set
#' @param scaling_method Character, how to scale delta for inactive treatments:
#'   - "zero": Set to 1.0 (no shift, **paper's canonical definition**)
#'   - "proportional": Scale proportionally to distance from boundary (smooth extension)
#'   - "fixed": Use fixed fraction of target delta (extension)
#' @param inactive_fraction Numeric, fraction of target delta for inactive treatments
#'   (used with scaling_method = "fixed")
#'
#' @return Numeric vector (length m) of adaptive deltas for this observation
#'
#' @details
#' **Adaptive delta rationale:**
#'
#' The adaptive estimator uses fixed epsilon but adapts delta_k(X) based on
#' propensity patterns. Intuition:
#' - Active treatments (m_k(X) ∈ [ε, 1-ε]): Use target delta
#' - Inactive treatments (m_k(X) < ε or > 1-ε): Set to 1.0 (paper) or scale (extensions)
#'
#' This satisfies the localized margin condition while allowing flexible
#' constraints per observation.
#'
#' **Scaling methods:**
#' - **zero**: δ_k(X) = 1.0 for inactive treatments (**paper's definition**)
#'   Binary active/inactive: δ_k for active, 1.0 (no shift) for inactive
#' - **proportional**: δ_k(X) = δ_target * min(m_k / ε, (1-m_k) / ε, 1)
#'   Smoothly scales delta based on distance from boundary (extension)
#' - **fixed**: δ_k(X) = δ_target * inactive_fraction for inactive treatments
#'   Simple scaling by fixed fraction (extension)
#'
#' @keywords internal
compute_adaptive_delta_single <- function(marginal_probs, delta_target, active_indices,
                                          epsilon, scaling_method = "proportional",
                                          inactive_fraction = 0.5) {

  m <- length(marginal_probs)

  # Validate inputs
  if (length(delta_target) != m) {
    stop("delta_target must have length m")
  }
  if (any(delta_target <= 0)) {
    stop("delta_target must be positive")
  }
  if (epsilon <= 0 || epsilon >= 0.5) {
    stop("epsilon must be in (0, 0.5)")
  }

  # Initialize adaptive delta with target values
  delta_adaptive <- delta_target

  # Identify inactive treatments
  inactive_idx <- setdiff(1:m, active_indices)

  if (length(inactive_idx) == 0) {
    # All treatments active: use target delta
    return(delta_adaptive)
  }

  # Scale delta for inactive treatments based on method
  if (scaling_method == "proportional") {
    # Proportional scaling based on distance from boundary. CRITICAL: delta = 1 is
    # the no-op (multiply odds by 1). Turning the intervention OFF as m_k -> boundary
    # means interpolating delta toward 1, NOT toward 0. Scaling the delta VALUE toward
    # 0 (the previous behavior) produces an EXTREME odds-shrinking shift on the most
    # degenerate coordinate -- an ill-conditioned KL solve, the opposite of the intent.
    for (k in inactive_idx) {
      m_k <- marginal_probs[k]
      # Distance scaling factor in [0,1]: 1 at the active edge, 0 at the boundary.
      # - If m_k < ε: scale by m_k/ε (closer to 0 → closer to no-op)
      # - If m_k > 1-ε: scale by (1-m_k)/ε (closer to 1 → closer to no-op)
      if (m_k < epsilon) {
        scale_factor <- m_k / epsilon
      } else if (m_k > 1 - epsilon) {
        scale_factor <- (1 - m_k) / epsilon
      } else {
        # Should not happen (inactive but in [ε, 1-ε])
        scale_factor <- 1
      }
      # Interpolate toward 1 (no shift) as scale_factor -> 0.
      delta_adaptive[k] <- 1 + (delta_target[k] - 1) * max(scale_factor, 0)
    }
  } else if (scaling_method == "fixed") {
    # Fixed fraction scaling
    delta_adaptive[inactive_idx] <- delta_target[inactive_idx] * inactive_fraction
  } else if (scaling_method == "zero") {
    # Set inactive deltas to 1.0 (no shift, paper's canonical definition)
    # δ_k = 1.0 means multiply odds by 1 (leave unchanged)
    delta_adaptive[inactive_idx] <- 1.0
  } else {
    stop("Unknown scaling_method: ", scaling_method)
  }

  return(delta_adaptive)
}


#' Compute Adaptive Deltas for All Observations
#'
#' Computes observation-specific delta_k(X) for all observations in dataset.
#'
#' @inheritParams compute_adaptive_delta_single
#' @param marginal_probs Matrix (n x m) of marginal propensities
#' @param active_sets List from \code{\link{compute_active_sets}}
#'
#' @return Matrix (n x m) of adaptive deltas
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data
#' n <- 100
#' m <- 3
#' ps <- matrix(runif(n * 2^m), nrow = n)
#' ps <- ps / rowSums(ps)
#'
#' # Extract marginal propensities
#' marginal_probs <- extract_marginal_probabilities(ps, m = m)
#'
#' # Compute active sets with fixed epsilon
#' active <- compute_active_sets(marginal_probs, epsilon = 0.1)
#'
#' # Compute adaptive deltas
#' delta_target <- c(0.1, 0.1, 0.1)
#' delta_adaptive <- compute_adaptive_deltas(
#'   marginal_probs = marginal_probs,
#'   delta_target = delta_target,
#'   active_sets = active,
#'   epsilon = 0.1,
#'   scaling_method = "proportional"
#' )
#'
#' # Check: active treatments get target delta, inactive get scaled
#' head(delta_adaptive)
#' }
compute_adaptive_deltas <- function(marginal_probs, delta_target, active_sets,
                                    epsilon, scaling_method = "proportional",
                                    inactive_fraction = 0.5) {

  n <- nrow(marginal_probs)
  m <- ncol(marginal_probs)

  if (length(delta_target) != m) {
    stop("delta_target must have length m")
  }

  # Initialize result matrix
  delta_adaptive_matrix <- matrix(0, nrow = n, ncol = m)

  # Compute adaptive delta for each observation
  for (i in 1:n) {
    delta_adaptive_matrix[i, ] <- compute_adaptive_delta_single(
      marginal_probs = marginal_probs[i, ],
      delta_target = delta_target,
      active_indices = active_sets$active_indices[[i]],
      epsilon = epsilon,
      scaling_method = scaling_method,
      inactive_fraction = inactive_fraction
    )
  }

  return(delta_adaptive_matrix)
}
