#' @file treatment_subset_utils.R
#' @title Utility Functions for Working with Treatment Subsets
#' @description
#' Shared utilities for working with reduced treatment sets. Used by both
#' eligibility stratification (where treatments are subset by eligibility pattern)
#' and active set thresholding (where treatments are subset by marginal propensity
#' boundaries). Provides consistent infrastructure for treatment vector reduction,
#' constraint Jacobian computation, and EIF extension.
#'
#' @section Architectural Separation:
#' These utilities operate at the **matrix level** (treatment combinations, propensities):
#' \itemize{
#'   \item Binary treatment matrices: 2^m rows × m columns
#'   \item Propensity vectors/matrices: n observations × 2^m combinations
#'   \item Constraint Jacobians: (|active|+1) × (|active|+1)
#' }
#'
#' For **data.frame level** operations (observation subsetting, column renaming),
#' see \code{\link{estimate_pattern_effect}} in estimation.R, which handles
#' eligibility-based data preparation before calling estimation functions.
#'
#' @author Denis Agniel

#' Create reduced treatment combination matrix
#'
#' @description
#' Extracts columns from full treatment combination matrix corresponding to
#' active/eligible treatments. Used when estimating on a subset of treatments.
#'
#' @param vectors Matrix (2^m x m), full treatment combinations in binary order
#' @param active_indices Integer vector, indices of active treatments (e.g., c(1,3,5))
#'
#' @return Matrix (2^m x |active|) containing only columns for active treatments
#'
#' @details
#' When working with a subset of treatments (either due to eligibility restrictions
#' or active set thresholding), we need to extract the relevant columns from the
#' full treatment combination matrix. This function does not change the number of
#' rows (still all 2^m combinations), but reduces columns to only those needed.
#'
#' @examples
#' \dontrun{
#' # Full treatment matrix for m=4
#' vectors <- generate_binary_vectors(4)  # 16 x 4
#'
#' # Only interested in treatments 2 and 4
#' active_indices <- c(2, 4)
#' vectors_reduced <- create_reduced_treatment_matrix(vectors, active_indices)
#' # Result: 16 x 2 matrix with columns for A2 and A4 only
#' }
#'
#' @export
create_reduced_treatment_matrix <- function(vectors, active_indices) {
  # Input validation
  if (!is.matrix(vectors)) {
    stop("vectors must be a matrix")
  }

  if (!is.numeric(active_indices) || any(active_indices < 1) ||
      any(active_indices > ncol(vectors))) {
    stop("active_indices must be integers between 1 and ncol(vectors)")
  }

  if (length(active_indices) == 0) {
    stop("active_indices cannot be empty")
  }

  # Extract columns
  vectors_reduced <- vectors[, active_indices, drop = FALSE]

  # Preserve column names if they exist
  if (!is.null(colnames(vectors))) {
    colnames(vectors_reduced) <- colnames(vectors)[active_indices]
  }

  return(vectors_reduced)
}


#' Create intercept plus reduced treatment matrix (a-bar)
#'
#' @description
#' Creates the "a-bar" matrix used in constraint Jacobian computation: an intercept
#' column (all 1s) plus columns for active treatments only. This is the matrix that
#' appears in Omega = sum_a pi_a * h'(a;lambda) * a_bar * a_bar^T.
#'
#' @param vectors Matrix (2^m x m), full treatment combinations
#' @param active_indices Integer vector, indices of active treatments
#'
#' @return Matrix (2^m x (|active|+1)) with intercept as first column, then active treatments
#'
#' @details
#' The constraint Jacobian Omega is computed using "a-bar" vectors that include
#' an intercept term plus the active treatment indicators. This function constructs
#' that matrix for all treatment combinations. The first column is always 1
#' (intercept), followed by columns for each active treatment.
#'
#' Dimension: If |active| = k, returns (2^m x (k+1)) matrix
#'
#' @examples
#' \dontrun{
#' vectors <- generate_binary_vectors(3)  # 8 x 3
#' active_indices <- c(1, 3)  # Only A1 and A3 active
#'
#' a_bar <- create_a_bar_matrix(vectors, active_indices)
#' # Result: 8 x 3 matrix
#' # Column 1: all 1s (intercept)
#' # Column 2: A1 values
#' # Column 3: A3 values
#' }
#'
#' @export
create_a_bar_matrix <- function(vectors, active_indices) {
  # Get reduced treatment matrix
  vectors_reduced <- create_reduced_treatment_matrix(vectors, active_indices)

  # Add intercept as first column
  n_combinations <- nrow(vectors_reduced)
  intercept <- rep(1, n_combinations)

  a_bar <- cbind(intercept, vectors_reduced)

  # Set column names
  colnames(a_bar) <- c("intercept",
                       paste0("A", active_indices))

  return(a_bar)
}


#' Extract marginal treatment probabilities from joint propensities
#'
#' @description
#' Computes marginal probabilities m_k(X) = P(A_k=1|X) from joint propensity
#' scores pi_a(X) = P(A=a|X). Used to determine which treatments are "active"
#' based on whether m_k(X) is near 0 or 1.
#'
#' @param propensity_scores Matrix (n x 2^m), joint propensity scores for all
#'   treatment combinations. Rows are observations, columns are combinations.
#'   Must sum to 1 across each row.
#' @param m Integer, number of treatments
#'
#' @return Matrix (n x m), marginal probabilities where [i,k] = P(A_k=1|X_i)
#'
#' @details
#' For each treatment k, the marginal probability is computed by summing the
#' joint probabilities over all treatment combinations where A_k = 1:
#'
#' m_k(X_i) = sum_{a: a_k=1} pi_a(X_i)
#'
#' This is used in active set thresholding to identify which treatments have
#' marginal propensities in (epsilon, 1-epsilon) and are therefore "active" for
#' that observation.
#'
#' @examples
#' \dontrun{
#' # For m=2, n=100
#' ps <- matrix(c(0.4, 0.3, 0.2, 0.1), nrow=100, ncol=4, byrow=TRUE)
#' marginals <- extract_marginal_probabilities(ps, m=2)
#'
#' # marginals[,1] = P(A1=1) = ps[,2] + ps[,4] (combos (0,1) and (1,1))
#' # marginals[,2] = P(A2=1) = ps[,3] + ps[,4] (combos (1,0) and (1,1))
#' }
#'
#' @export
extract_marginal_probabilities <- function(propensity_scores, m) {
  # Input validation
  if (!is.matrix(propensity_scores)) {
    stop("propensity_scores must be a matrix")
  }

  n <- nrow(propensity_scores)
  expected_cols <- 2^m

  if (ncol(propensity_scores) != expected_cols) {
    stop(sprintf("propensity_scores should have 2^m = %d columns, has %d",
                 expected_cols, ncol(propensity_scores)))
  }

  # Check rows sum to 1 (within tolerance)
  row_sums <- rowSums(propensity_scores)
  if (!all(abs(row_sums - 1) < 1e-6)) {
    max_deviation <- max(abs(row_sums - 1))
    warning(sprintf("Propensity scores do not sum to 1 (max deviation: %.2e)",
                    max_deviation))
  }

  # Generate treatment combinations to map columns to marginals
  vectors <- generate_binary_vectors(m)

  # Compute marginals: m_k = sum over combinations where a_k = 1
  marginals <- matrix(0, nrow = n, ncol = m)

  for (k in 1:m) {
    # Find which combinations have A_k = 1
    combos_with_k <- which(vectors[, k] == 1)

    # Sum propensities for those combinations
    marginals[, k] <- rowSums(propensity_scores[, combos_with_k, drop = FALSE])
  }

  # Set column names
  colnames(marginals) <- paste0("m", 1:m)

  return(marginals)
}


#' Compute constraint Jacobian (Omega) for reduced treatment set
#'
#' @description
#' Computes the constraint Jacobian matrix Omega for a single observation when
#' working with a reduced set of treatments. Used in EIF computation when some
#' treatments are dropped from the active set.
#'
#' @param propensity_scores Numeric vector (length 2^m), propensity scores for
#'   all treatment combinations for ONE observation
#' @param lambda Numeric vector (length |active|+1), Lagrange multipliers for
#'   this observation (intercept + active treatments)
#' @param a_bar_matrix Matrix (2^m x (|active|+1)), from create_a_bar_matrix()
#' @param divergence Character, divergence type ("kl", "chi2", "hellinger", "reverse_kl")
#'
#' @return Matrix ((|active|+1) x (|active|+1)), constraint Jacobian
#'
#' @details
#' The constraint Jacobian is:
#'
#' Omega = sum_{a in A} pi_a * h'(a_bar; lambda) * a_bar * a_bar^T
#'
#' where h'() is the derivative of the h-function with respect to the linear term.
#' The dimension depends on the number of active constraints: if |active| = k,
#' then Omega is (k+1) x (k+1).
#'
#' This matrix must be inverted to compute the EIF. It is guaranteed to be
#' positive definite under the interiority and no-redundant-constraints assumptions.
#'
#' @examples
#' \dontrun{
#' # Single observation with m=3, active set = {1,3}
#' ps_obs <- c(0.1, 0.2, 0.15, 0.25, 0.05, 0.1, 0.1, 0.05)  # length 8
#' lambda_obs <- c(0.5, 0.3, -0.2)  # intercept + 2 active treatments
#'
#' vectors <- generate_binary_vectors(3)
#' a_bar <- create_a_bar_matrix(vectors, c(1, 3))  # 8 x 3
#'
#' omega <- compute_omega_reduced(ps_obs, lambda_obs, a_bar, "kl")
#' # Result: 3 x 3 positive definite matrix
#' }
#'
#' @export
compute_omega_reduced <- function(propensity_scores, lambda, a_bar_matrix, divergence) {
  # Input validation
  if (!is.numeric(propensity_scores) || !is.vector(propensity_scores)) {
    stop("propensity_scores must be a numeric vector")
  }

  if (!is.numeric(lambda) || !is.vector(lambda)) {
    stop("lambda must be a numeric vector")
  }

  if (!is.matrix(a_bar_matrix)) {
    stop("a_bar_matrix must be a matrix")
  }

  n_combinations <- length(propensity_scores)
  dim_lambda <- length(lambda)

  if (nrow(a_bar_matrix) != n_combinations) {
    stop("a_bar_matrix rows must equal length of propensity_scores")
  }

  if (ncol(a_bar_matrix) != dim_lambda) {
    stop("a_bar_matrix columns must equal length of lambda")
  }

  if (!divergence %in% c("kl", "chi2", "hellinger", "reverse_kl")) {
    stop("divergence must be one of: kl, chi2, hellinger, reverse_kl")
  }

  # Compute linear terms: lambda^T * a_bar for each combination
  linear_terms <- as.vector(a_bar_matrix %*% lambda)

  # Compute h'(linear_term) for each combination based on divergence
  if (divergence == "kl") {
    # h(z) = exp(z), so h'(z) = exp(z)
    h_prime <- exp(linear_terms)

  } else if (divergence == "chi2") {
    # h(z) = z/2 + 1, so h'(z) = 1/2
    h_prime <- rep(0.5, n_combinations)

  } else if (divergence == "reverse_kl") {
    # [BUGFIX 2026-08-11] This branch previously used h(z) = 1/(1-z),
    # h'(z) = 1/(1-z)^2 -- the WRONG conjugate for reverse KL. The canonical,
    # paper-derived form (see divergence_h_pair() in divergence_h_functions.R,
    # and compute_q_delta_weights_internal() in fast_qdelta_solver.R) is
    # h(x) = -1/x, h'(x) = 1/x^2 on the domain x < 0. `linear_terms` here is
    # exactly the same quantity x = lambda' a_bar that the canonical function
    # takes, so the two are directly comparable -- and they disagreed by factors
    # of 2.25x (x = -2) to 121x (x = -0.1) across the whole region the solver
    # actually reaches, i.e. this was wrong for every reverse_kl call, not just
    # near a boundary. Delegating to the canonical primitive rather than
    # rewriting the formula inline, so this site cannot drift again.
    # Blast radius was checked: no stored application result was computed with
    # divergence = "reverse_kl" through this function (see
    # quality_reports/2026-08-11_omega-reduced-reverse-kl-bugfix-report.md).
    hp <- divergence_h_pair(linear_terms, "reverse_kl")
    h_prime <- as.vector(hp$h_prime)

    # Preserve this function's diagnostic contract (it warned when its guard
    # engaged), but on the canonical clipping condition: x pinned to the
    # interior of the domain x < 0 because an iterate came within
    # REVERSE_KL_DOMAIN_BUFFER of the singularity at x = 0, or below
    # REVERSE_KL_LINEAR_TERM_MIN.
    if (hp$n_clipped > 0L) {
      warning("Linear terms clipped to the reverse_kl domain (x < 0) in ",
              hp$n_clipped, " of ", n_combinations, " combinations when ",
              "computing Omega; h' saturated there.", call. = FALSE)
    }

  } else if (divergence == "hellinger") {
    # h(z) = 1/(1-z)^2, so h'(z) = 2/(1-z)^3
    denominator <- (1 - linear_terms)^3

    # Guard against division issues
    if (any(abs(denominator) < 1e-10)) {
      warning("Near-singular h' values in hellinger; clamping denominators")
      denominator <- pmax(abs(denominator), 1e-10) * sign(denominator)
    }

    h_prime <- 2 / denominator
  }

  # Compute Omega = sum_a pi_a * h'_a * a_bar_a * a_bar_a^T
  # Use vectorized computation for efficiency
  weights <- propensity_scores * h_prime  # length n_combinations

  # Omega = t(a_bar_matrix) %*% diag(weights) %*% a_bar_matrix
  # More efficiently: t(a_bar_matrix * weights) %*% a_bar_matrix
  weighted_a_bar <- a_bar_matrix * weights  # broadcast weights across columns
  omega <- t(weighted_a_bar) %*% a_bar_matrix

  # Add small regularization for numerical stability (if needed)
  # Only if matrix appears near-singular
  min_eigenvalue <- min(eigen(omega, only.values = TRUE)$values)
  if (min_eigenvalue < 1e-8) {
    omega <- omega + diag(1e-6, nrow(omega))
  }

  return(omega)
}


#' Extend reduced-dimension EIF values to full dataset
#'
#' @description
#' When estimating within a subset of observations (e.g., one eligibility stratum),
#' the EIF is computed only for that subset. This function extends those values
#' to the full dataset dimension by filling in zeros for observations not in the subset.
#'
#' @param eif_reduced Numeric vector, EIF values for subset of observations
#' @param indices_in_full Integer vector, indices mapping eif_reduced to full dataset.
#'   Must have length equal to length(eif_reduced).
#' @param n_full Integer, size of full dataset
#'
#' @return Numeric vector (length n_full) with EIF values in correct positions and
#'   zeros elsewhere
#'
#' @details
#' Used in stratified estimation where we estimate within strata and need to
#' aggregate across strata. Each stratum produces EIF values for its observations;
#' this function places those values in the full dataset indexing.
#'
#' For observations not in the subset, EIF = 0, which is correct because those
#' observations don't contribute to the estimand for that stratum.
#'
#' @examples
#' \dontrun{
#' # Suppose we have n=1000, but only estimate on stratum of n=300
#' eif_stratum <- rnorm(300)  # EIF for 300 observations
#' stratum_indices <- c(50:349)  # These obs are in the stratum
#'
#' eif_full <- extend_eif_to_full(eif_stratum, stratum_indices, 1000)
#' # Result: length 1000 vector
#' #   eif_full[50:349] = eif_stratum
#' #   eif_full[elsewhere] = 0
#' }
#'
#' @export
extend_eif_to_full <- function(eif_reduced, indices_in_full, n_full) {
  # Input validation
  if (!is.numeric(eif_reduced)) {
    stop("eif_reduced must be numeric")
  }

  if (!is.numeric(indices_in_full) || any(indices_in_full < 1) ||
      any(indices_in_full > n_full)) {
    stop("indices_in_full must be integers between 1 and n_full")
  }

  if (length(eif_reduced) != length(indices_in_full)) {
    stop("eif_reduced and indices_in_full must have same length")
  }

  if (length(unique(indices_in_full)) != length(indices_in_full)) {
    stop("indices_in_full contains duplicates")
  }

  # Create full vector initialized to zero
  eif_full <- numeric(n_full)

  # Fill in values at specified indices
  eif_full[indices_in_full] <- eif_reduced

  return(eif_full)
}


#' Map treatment indices between full and reduced sets
#'
#' @description
#' Creates explicit mapping between original treatment indices (1 to m) and
#' reduced treatment indices (1 to |active|). Useful for clear documentation
#' and avoiding off-by-one errors when working with subsets.
#'
#' @param active_indices Integer vector, indices of active treatments in original indexing
#' @param m Integer, total number of treatments in full set (optional, for validation)
#'
#' @return Data frame with columns:
#'   - original_index: Original treatment number (1 to m)
#'   - reduced_index: Index in reduced set (1 to |active|)
#'   - treatment_name: Name (e.g., "A3")
#'
#' @details
#' When working with reduced treatment sets, it's helpful to have explicit mapping
#' between original and reduced indices. For example, if active treatments are {2,5,7},
#' the mapping is: 2→1, 5→2, 7→3.
#'
#' @examples
#' \dontrun{
#' active <- c(2, 5, 7)
#' mapping <- map_treatment_indices(active, m = 10)
#' # Result:
#' #   original_index reduced_index treatment_name
#' #              2             1            A2
#' #              5             2            A5
#' #              7             3            A7
#' }
#'
#' @export
map_treatment_indices <- function(active_indices, m = NULL) {
  # Input validation
  if (!is.numeric(active_indices) || any(active_indices < 1)) {
    stop("active_indices must be positive integers")
  }

  if (!is.null(m) && any(active_indices > m)) {
    stop("active_indices contains values > m")
  }

  # Create mapping
  mapping <- data.frame(
    original_index = active_indices,
    reduced_index = seq_along(active_indices),
    treatment_name = paste0("A", active_indices)
  )

  return(mapping)
}


#' Validate active set configuration
#'
#' @description
#' Checks that an active set configuration is valid and provides diagnostics.
#' Used to catch common errors when working with treatment subsets.
#'
#' @param active_indices Integer vector, proposed active indices
#' @param m Integer, total number of treatments
#' @param min_active Integer, minimum number of active treatments required (default: 1)
#'
#' @return List with validation results:
#'   - valid: Logical, whether configuration is valid
#'   - n_active: Number of active treatments
#'   - issues: Character vector of any problems found
#'   - warnings: Character vector of potential issues
#'
#' @details
#' Checks for:
#' - Empty active set (usually an error)
#' - Out-of-range indices
#' - Duplicates
#' - Very few active treatments (may be valid but worth flagging)
#'
#' @examples
#' \dontrun{
#' validation <- validate_active_set(c(2, 5), m = 10)
#' if (!validation$valid) {
#'   stop("Invalid active set: ", paste(validation$issues, collapse = "; "))
#' }
#' }
#'
#' @export
validate_active_set <- function(active_indices, m, min_active = 1) {
  issues <- character(0)
  warnings <- character(0)

  # Check empty
  if (length(active_indices) == 0) {
    issues <- c(issues, "Active set is empty")
  }

  # Check range
  if (any(active_indices < 1) || any(active_indices > m)) {
    issues <- c(issues, sprintf("Indices out of range [1, %d]", m))
  }

  # Check duplicates
  if (length(unique(active_indices)) != length(active_indices)) {
    issues <- c(issues, "Active set contains duplicates")
  }

  # Check minimum size
  if (length(active_indices) < min_active) {
    issues <- c(issues, sprintf("Fewer than %d active treatments", min_active))
  }

  # Warnings for edge cases
  if (length(active_indices) == m) {
    warnings <- c(warnings, "All treatments active (no reduction)")
  }

  if (length(active_indices) == 1) {
    warnings <- c(warnings, "Only 1 active treatment (minimal model)")
  }

  valid <- length(issues) == 0

  return(list(
    valid = valid,
    n_active = length(active_indices),
    issues = issues,
    warnings = warnings
  ))
}
