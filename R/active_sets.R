#' Compute Vanishing Threshold for Active Set Selection
#'
#' Computes the vanishing threshold `ε_n = C * n^{-γ}` for active set selection.
#' The threshold decreases with sample size at rate controlled by γ.
#'
#' @param n Integer. Sample size.
#' @param gamma Numeric. Rate parameter controlling threshold decay (`ε_n = O(n^{-γ})`).
#' @param constant Numeric. Multiplicative constant (default = 1). Controls base level.
#'
#' @return Numeric. The computed threshold ε_n.
#'
#' @details
#' The vanishing threshold is used to define the active set:
#'   `K_ε(X) = {k : m_k(X) ∈ [ε_n, 1-ε_n]}`
#'
#' where m_k(X) = P(A_k = 1 | X) is the marginal propensity for treatment k.
#'
#' The threshold depends only on the sample size `n` and the user-chosen rate
#' `gamma`. The asymptotic-normality theory (Corollary in the paper) requires
#' `gamma > 1/(2(1+β))`, where β is the margin exponent of the data-generating
#' law — an object that is not estimable from data in practice. That condition
#' governs the *asymptotic* guarantee, not whether the estimate is computable, so
#' it is deliberately NOT enforced here: the estimator is well-defined for any
#' valid `gamma`, and larger `gamma` (faster-vanishing threshold) is the practical
#' knob a user controls.
#'
#' @examples
#' # With n=1000, γ=0.3: ε_n ≈ 0.050
#' compute_vanishing_threshold(1000, gamma = 0.3)
#'
#' @export
compute_vanishing_threshold <- function(n, gamma, constant = 1) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("n must be a positive integer")
  }
  if (!is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("gamma must be a positive number")
  }
  if (!is.numeric(constant) || length(constant) != 1 || constant <= 0) {
    stop("constant must be a positive number")
  }

  # Compute threshold
  epsilon_n <- constant * n^(-gamma)

  # Ensure threshold is in (0, 0.5)
  if (epsilon_n >= 0.5) {
    warning(sprintf(
      "Threshold ε_n = %.4f >= 0.5 (n=%d, γ=%.3f). Consider larger γ or smaller constant.",
      epsilon_n, n, gamma
    ))
    epsilon_n <- 0.49  # Cap at just below 0.5
  }

  return(epsilon_n)
}


#' Compute Active Sets Based on Marginal Propensities
#'
#' Identifies which treatments are "active" (have propensities away from boundaries)
#' for each observation. A treatment k is active if `m_k(X) ∈ [ε, 1-ε]`.
#'
#' @param marginal_probs Numeric matrix of size n x m. Marginal propensity scores.
#'   Element `[i,k]` is P(A_k = 1 | X_i).
#' @param epsilon Numeric. Threshold defining active set (scalar or vector length n).
#'   Treatments with marginal probabilities in `[ε, 1-ε]` are considered active.
#' @param min_active Integer. Minimum number of active treatments required per observation
#'   (default = 1). Observations with fewer active treatments will have all treatments
#'   marked as inactive (active_matrix will be all FALSE for that row).
#'
#' @return List with components:
#'   \item{active_matrix}{Logical matrix (n x m). TRUE if treatment k is active for obs i.}
#'   \item{n_active}{Integer vector (length n). Number of active treatments per observation.}
#'   \item{epsilon}{Numeric. The threshold used (scalar or vector).}
#'   \item{sufficient_active}{Logical vector (length n). TRUE if observation has ≥ min_active
#'     active treatments.}
#'   \item{active_indices}{List of length n. Element i contains integer vector of active
#'     treatment indices for observation i.}
#'
#' @details
#' **Active set definition:**
#'   `K_ε(X_i) = {k : m_k(X_i) ∈ [ε, 1-ε]}`
#'
#' Treatments near the boundary (m_k close to 0 or 1) are excluded because:
#' - Nearly deterministic treatments provide little information
#' - Constraint Jacobian Ω may become singular
#' - Requires stronger overlap assumptions
#'
#' **Handling insufficient active treatments:**
#' If an observation has fewer than `min_active` active treatments, we mark
#' all treatments as inactive for that observation. This allows calling code
#' to handle such observations specially (e.g., use P_0 distribution, warn user).
#'
#' @examples
#' # Simulate marginal propensities for 100 obs, 4 treatments
#' set.seed(123)
#' m_probs <- matrix(runif(400, 0.1, 0.9), nrow = 100, ncol = 4)
#'
#' # Compute active sets with ε = 0.05
#' active <- compute_active_sets(m_probs, epsilon = 0.05)
#'
#' # Check how many observations have all 4 treatments active
#' table(active$n_active)
#'
#' # Check which observations have insufficient active treatments
#' sum(!active$sufficient_active)
#'
#' @export
compute_active_sets <- function(marginal_probs, epsilon, min_active = 1) {
  # Input validation
  if (!is.matrix(marginal_probs) || !is.numeric(marginal_probs)) {
    stop("marginal_probs must be a numeric matrix")
  }
  n <- nrow(marginal_probs)
  m <- ncol(marginal_probs)

  if (any(marginal_probs < 0 | marginal_probs > 1)) {
    stop("marginal_probs must contain probabilities in [0,1]")
  }

  if (!is.numeric(epsilon) || any(epsilon <= 0) || any(epsilon >= 0.5)) {
    stop("epsilon must be in (0, 0.5)")
  }

  if (length(epsilon) != 1 && length(epsilon) != n) {
    stop("epsilon must be scalar or vector of length n")
  }

  if (!is.numeric(min_active) || length(min_active) != 1 ||
      min_active < 1 || min_active != floor(min_active)) {
    stop("min_active must be a positive integer")
  }

  if (min_active > m) {
    stop(sprintf("min_active (%d) cannot exceed number of treatments (%d)", min_active, m))
  }

  # Broadcast epsilon if scalar
  if (length(epsilon) == 1) {
    epsilon_vec <- rep(epsilon, n)
  } else {
    epsilon_vec <- epsilon
  }

  # Compute active set indicators
  # Treatment k is active for obs i if m_k(X_i) ∈ [ε_i, 1-ε_i]
  active_matrix <- matrix(FALSE, nrow = n, ncol = m)
  for (i in 1:n) {
    eps_i <- epsilon_vec[i]
    active_matrix[i, ] <- (marginal_probs[i, ] >= eps_i) &
                          (marginal_probs[i, ] <= (1 - eps_i))
  }

  # Count active treatments per observation
  n_active <- rowSums(active_matrix)

  # Check sufficiency
  sufficient_active <- n_active >= min_active

  # For observations with insufficient active treatments, mark all as inactive
  if (any(!sufficient_active)) {
    active_matrix[!sufficient_active, ] <- FALSE
    n_active[!sufficient_active] <- 0
  }

  # Extract active indices for each observation
  active_indices <- vector("list", n)
  for (i in 1:n) {
    if (sufficient_active[i]) {
      active_indices[[i]] <- which(active_matrix[i, ])
    } else {
      active_indices[[i]] <- integer(0)
    }
  }

  return(list(
    active_matrix = active_matrix,
    n_active = n_active,
    epsilon = epsilon,
    sufficient_active = sufficient_active,
    active_indices = active_indices
  ))
}


#' Summarize Active Set Statistics
#'
#' Provides descriptive statistics about active sets across observations.
#' Useful for diagnosing whether threshold is appropriate and checking
#' whether sufficient treatments remain active.
#'
#' @param active_sets List returned by \code{\link{compute_active_sets}}.
#' @param marginal_probs Numeric matrix (n x m) of marginal propensities (optional).
#'   If provided, computes additional diagnostics about propensities.
#'
#' @return List with components:
#'   \item{n_obs}{Integer. Number of observations.}
#'   \item{m_treatments}{Integer. Number of treatments.}
#'   \item{epsilon}{Numeric. Threshold used.}
#'   \item{n_active_summary}{Numeric vector. Summary statistics (min, Q1, median, Q3, max)
#'     of number of active treatments per observation.}
#'   \item{n_insufficient}{Integer. Number of observations with insufficient active treatments.}
#'   \item{prop_insufficient}{Numeric. Proportion of observations with insufficient active.}
#'   \item{active_by_treatment}{Integer vector (length m). Number of observations where
#'     each treatment is active.}
#'   \item{prop_active_by_treatment}{Numeric vector (length m). Proportion of observations
#'     where each treatment is active.}
#'   \item{min_propensity}{Numeric (length m, optional). Minimum marginal propensity per treatment.}
#'   \item{max_propensity}{Numeric (length m, optional). Maximum marginal propensity per treatment.}
#'
#' @examples
#' # Generate data
#' set.seed(123)
#' m_probs <- matrix(runif(400, 0.1, 0.9), nrow = 100, ncol = 4)
#' active <- compute_active_sets(m_probs, epsilon = 0.05)
#'
#' # Summarize
#' summarize_active_sets(active, marginal_probs = m_probs)
#'
#' @export
summarize_active_sets <- function(active_sets, marginal_probs = NULL) {
  n <- length(active_sets$n_active)
  m <- ncol(active_sets$active_matrix)

  # Basic counts
  n_active_summary <- summary(active_sets$n_active)
  n_insufficient <- sum(!active_sets$sufficient_active)
  prop_insufficient <- n_insufficient / n

  # Extract individual summary statistics
  min_n_active <- as.numeric(n_active_summary["Min."])
  max_n_active <- as.numeric(n_active_summary["Max."])
  mean_n_active <- as.numeric(n_active_summary["Mean"])
  median_n_active <- as.numeric(n_active_summary["Median"])

  # Proportion with all treatments active
  n_all_active <- sum(active_sets$n_active == m)
  pct_all_active <- n_all_active / n

  # By treatment
  active_by_treatment <- colSums(active_sets$active_matrix)
  prop_active_by_treatment <- active_by_treatment / n

  result <- list(
    n_obs = n,
    m = m,  # Short name for user convenience
    epsilon = active_sets$epsilon,
    n_active_summary = n_active_summary,  # Keep full summary for backwards compat
    min_n_active = min_n_active,
    max_n_active = max_n_active,
    mean_n_active = mean_n_active,
    median_n_active = median_n_active,
    pct_all_active = pct_all_active,
    n_insufficient = n_insufficient,
    prop_insufficient = prop_insufficient,
    active_by_treatment = active_by_treatment,
    prop_active_by_treatment = prop_active_by_treatment
  )

  # Add propensity diagnostics if provided
  if (!is.null(marginal_probs)) {
    result$min_propensity <- apply(marginal_probs, 2, min)
    result$max_propensity <- apply(marginal_probs, 2, max)
  }

  return(result)
}
