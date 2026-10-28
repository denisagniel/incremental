#' @file adaptive_constraints.R
#' @title Adaptive Constraint Management for Singular Systems
#' @description
#' Utilities for detecting and handling singular Jacobian matrices in Q_delta
#' optimization. When the constraint Jacobian Ω becomes singular, this module:
#' 1. Drops truly redundant constraints (`P(A_k=1) ∈ {0,1}`)
#' 2. Uses regularization for active constraints causing near-boundary singularity
#' This preserves the estimand while allowing estimation at/near the boundary.
#' @author Denis Agniel

#' Check if a marginal constraint is truly redundant
#'
#' A constraint k is redundant if P(A_k = 1) is deterministic (0 or 1).
#' In this case, λ_k is not identified and can be dropped without changing Q_δ.
#'
#' @param constraint_idx Index of constraint to check (1 = normalization, 2:(m+1) = marginals)
#' @param treatment_vectors Matrix of treatment combinations (n_active × m)
#' @param propensity_scores Vector of propensity scores for active combinations
#' @param epsilon Threshold for considering marginal as deterministic (default 1e-6)
#' @return TRUE if constraint is redundant, FALSE if active
is_constraint_redundant <- function(constraint_idx, treatment_vectors,
                                   propensity_scores, epsilon = 1e-6) {

  # Normalization constraint (idx=1) is never redundant
  if (constraint_idx == 1) {
    return(FALSE)
  }

  # For marginal constraints (idx = k+1 for treatment k)
  k <- constraint_idx - 1  # Convert to treatment index (1-based)

  if (k > ncol(treatment_vectors)) {
    warning("Constraint index exceeds number of treatments")
    return(FALSE)
  }

  # Compute marginal: E_P0[A_k] = Σ_{a: a_k=1} π_a
  marginal_k <- sum(propensity_scores[treatment_vectors[, k] == 1])

  # Check if marginal is deterministic (0 or 1)
  is_redundant <- (marginal_k < epsilon) || (marginal_k > 1 - epsilon)

  return(is_redundant)
}


#' Check if a matrix is numerically singular
#'
#' Uses condition number and determinant to detect singularity
#'
#' @param mat Square matrix to check
#' @param cond_threshold Condition number threshold (default 1e10)
#' @param det_threshold Determinant threshold (default 1e-12)
#' @return List with is_singular (logical), condition_number, determinant, rank
is_matrix_singular <- function(mat, cond_threshold = 1e10, det_threshold = 1e-12) {

  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    stop("Input must be a square matrix")
  }

  n <- nrow(mat)

  # Compute condition number using singular value decomposition
  svd_result <- svd(mat)
  singular_values <- svd_result$d

  # Condition number = σ_max / σ_min
  cond_num <- if (min(singular_values) > 1e-15) {
    max(singular_values) / min(singular_values)
  } else {
    Inf
  }

  # Determinant
  det_val <- det(mat)

  # Numerical rank (count singular values > threshold)
  rank_val <- sum(singular_values > 1e-10 * max(singular_values))

  # Is singular?
  is_singular <- (cond_num > cond_threshold) ||
                 (abs(det_val) < det_threshold) ||
                 (rank_val < n)

  return(list(
    is_singular = is_singular,
    condition_number = cond_num,
    determinant = det_val,
    rank = rank_val,
    full_rank = n,
    singular_values = singular_values
  ))
}


#' Identify which constraint to drop from a singular system
#'
#' Uses several heuristics to identify the most problematic constraint:
#' 1. Smallest singular value contribution
#' 2. Row norm (constraints with very small norms are likely redundant)
#' 3. Linear dependence check
#'
#' @param jacobian Constraint Jacobian matrix (m+1 × m+1 or m × m)
#' @param constraint_vector Constraint vector F(λ)
#' @param method Method to identify redundant constraint ("svd", "norm", "rank")
#' @return Index of constraint to drop (1-indexed, excluding intercept)
identify_redundant_constraint <- function(jacobian, constraint_vector,
                                         method = c("svd", "norm", "rank")) {

  method <- match.arg(method)
  n <- nrow(jacobian)

  if (method == "svd") {
    # Use SVD to find which constraint contributes least to smallest singular value
    # The right singular vector corresponding to smallest singular value
    # shows which constraint is most redundant
    svd_result <- svd(jacobian)

    # Right singular vector for smallest singular value
    v_min <- svd_result$v[, n]  # Last column = smallest singular value

    # Don't drop the normalization constraint (first row in internal normalization)
    # For internal: constraint 1 is normalization, 2:(m+1) are marginals
    v_min_safe <- abs(v_min)
    v_min_safe[1] <- -Inf  # Protect normalization constraint

    # Constraint with largest absolute component in v_min is most redundant
    # (contributes most to the null space)
    redundant_idx <- which.max(v_min_safe)

  } else if (method == "norm") {
    # Find constraint (row) with smallest norm
    # Small norm suggests the constraint is nearly satisfied automatically
    row_norms <- apply(jacobian, 1, function(row) sqrt(sum(row^2)))

    # Don't drop the normalization constraint (first row in internal normalization)
    # For internal: constraint 1 is normalization, 2:(m+1) are marginals
    # For simplicity, consider dropping from constraints 2:n if n > 1
    if (n > 1) {
      redundant_idx <- which.min(row_norms[2:n]) + 1  # +1 to account for starting at 2
    } else {
      redundant_idx <- 1
    }

  } else if (method == "rank") {
    # Find which row, when removed, maximizes the rank of remaining matrix
    # This is expensive but most principled
    best_idx <- 1
    best_rank <- 0

    for (i in 2:n) {  # Don't try removing first constraint (normalization)
      reduced_jacobian <- jacobian[-i, -i, drop = FALSE]
      if (nrow(reduced_jacobian) > 0) {
        rank_i <- qr(reduced_jacobian)$rank
        if (rank_i > best_rank) {
          best_rank <- rank_i
          best_idx <- i
        }
      }
    }

    redundant_idx <- best_idx
  }

  return(redundant_idx)
}


#' Solve linear system with hybrid adaptive approach
#'
#' If the Jacobian is singular:
#' 1. First, drop truly redundant constraints (`P(A_k=1) ∈ {0,1}`)
#' 2. If still singular, use regularization (Ω + ρI) instead of dropping active constraints
#'
#' This preserves the estimand while handling both redundant constraints and
#' near-boundary singularity.
#'
#' @param jacobian Constraint Jacobian matrix
#' @param constraint_vector Constraint vector
#' @param lambda Current lambda values
#' @param treatment_vectors Matrix of treatment combinations (for redundancy check)
#' @param propensity_scores Vector of propensity scores (for redundancy check)
#' @param max_drops Maximum number of constraints to drop (default 2)
#' @param cond_threshold Condition number threshold for singularity (default 1e10)
#' @param regularization_rho Regularization parameter for active constraints (default 1e-6)
#' @param verbose Print diagnostic messages (default FALSE)
#' @return List with:
#'   - lambda_update: Newton step `Δλ = -J^{-1} F(λ)`
#'   - constraints_dropped: Indices of dropped constraints
#'   - n_drops: Number of constraints dropped
#'   - used_regularization: Whether regularization was applied
#'   - condition_number: Final condition number
#'   - converged: Whether a solution was found
adaptive_solve <- function(jacobian, constraint_vector, lambda,
                          treatment_vectors = NULL, propensity_scores = NULL,
                          max_drops = 2, cond_threshold = 1e10,
                          regularization_rho = 1e-6,
                          verbose = FALSE) {

  n <- nrow(jacobian)
  dropped_indices <- integer(0)
  used_regularization <- FALSE

  # Check if system is already solvable
  singularity_check <- is_matrix_singular(jacobian, cond_threshold = cond_threshold)

  if (!singularity_check$is_singular) {
    # System is fine, solve normally
    lambda_update <- tryCatch({
      solve(jacobian, constraint_vector)
    }, error = function(e) {
      # Numerical issues even though condition number looks OK
      # Fall back to pseudo-inverse
      MASS::ginv(jacobian) %*% constraint_vector
    })

    return(list(
      lambda_update = lambda_update,
      constraints_dropped = integer(0),
      n_drops = 0,
      used_regularization = FALSE,
      condition_number = singularity_check$condition_number,
      converged = TRUE
    ))
  }

  # System is singular - begin hybrid adaptive approach
  if (verbose) {
    message(sprintf("Jacobian singular (κ = %.2e, det = %.2e, rank = %d/%d)",
                    singularity_check$condition_number,
                    singularity_check$determinant,
                    singularity_check$rank,
                    singularity_check$full_rank))
    message("Starting hybrid adaptive approach...")
  }

  # ============================================================================
  # STEP 1: Identify and drop ONLY truly redundant constraints
  # ============================================================================

  redundant_constraints <- integer(0)

  if (!is.null(treatment_vectors) && !is.null(propensity_scores)) {
    # Check each constraint for redundancy
    for (idx in 1:n) {
      if (is_constraint_redundant(idx, treatment_vectors, propensity_scores)) {
        redundant_constraints <- c(redundant_constraints, idx)
        if (verbose) {
          if (idx == 1) {
            message(sprintf("  Constraint %d (normalization) is redundant", idx))
          } else {
            message(sprintf("  Constraint %d (marginal A%d) is redundant: P(A%d=1) ∈ {0,1}",
                           idx, idx-1, idx-1))
          }
        }
      }
    }
  }

  # Drop redundant constraints if found
  if (length(redundant_constraints) > 0) {
    if (verbose) {
      message(sprintf("  Dropping %d redundant constraint(s): %s",
                     length(redundant_constraints),
                     paste(redundant_constraints, collapse = ", ")))
    }

    # Remove redundant constraints
    jacobian_reduced <- jacobian[-redundant_constraints, -redundant_constraints, drop = FALSE]
    constraint_reduced <- constraint_vector[-redundant_constraints]
    dropped_indices <- redundant_constraints

    # Check if system is now solvable
    if (nrow(jacobian_reduced) > 0) {
      singularity_check_reduced <- is_matrix_singular(jacobian_reduced,
                                                      cond_threshold = cond_threshold)

      if (!singularity_check_reduced$is_singular) {
        # System is now solvable after dropping redundant constraints
        lambda_update_reduced <- tryCatch({
          solve(jacobian_reduced, constraint_reduced)
        }, error = function(e) {
          MASS::ginv(jacobian_reduced) %*% constraint_reduced
        })

        # Map back to full lambda space (dropped constraints set to 0)
        lambda_update_full <- numeric(n)
        keep_indices <- setdiff(1:n, dropped_indices)
        lambda_update_full[keep_indices] <- lambda_update_reduced

        if (verbose) {
          message(sprintf("  Success after dropping redundant constraints! κ = %.2e",
                         singularity_check_reduced$condition_number))
        }

        return(list(
          lambda_update = lambda_update_full,
          constraints_dropped = dropped_indices,
          n_drops = length(dropped_indices),
          used_regularization = FALSE,
          condition_number = singularity_check_reduced$condition_number,
          converged = TRUE
        ))
      } else {
        # Still singular after dropping redundant constraints
        # Move to Step 2 (regularization)
        jacobian <- jacobian_reduced
        constraint_vector <- constraint_reduced
        if (verbose) {
          message(sprintf("  Still singular after dropping redundant constraints (κ = %.2e)",
                         singularity_check_reduced$condition_number))
          message("  Switching to regularization for active constraints...")
        }
      }
    }
  } else {
    if (verbose) {
      message("  No redundant constraints found")
      message("  All constraints are active - using regularization...")
    }
  }

  # ============================================================================
  # STEP 2: Use regularization for remaining active constraints
  # ============================================================================

  # Apply regularization: Ω + ρI
  n_current <- nrow(jacobian)
  jacobian_regularized <- jacobian + regularization_rho * diag(n_current)

  singularity_check_reg <- is_matrix_singular(jacobian_regularized,
                                               cond_threshold = cond_threshold)

  lambda_update_regularized <- tryCatch({
    solve(jacobian_regularized, constraint_vector)
  }, error = function(e) {
    # Even regularization failed - use pseudo-inverse
    MASS::ginv(jacobian_regularized) %*% constraint_vector
  })

  # Map back to full lambda space if we dropped redundant constraints
  if (length(dropped_indices) > 0) {
    lambda_update_full <- numeric(n)
    keep_indices <- setdiff(1:n, dropped_indices)
    lambda_update_full[keep_indices] <- lambda_update_regularized
    final_cond_num <- singularity_check_reg$condition_number
  } else {
    lambda_update_full <- lambda_update_regularized
    final_cond_num <- singularity_check_reg$condition_number
  }

  if (verbose) {
    message(sprintf("  Regularization applied (ρ = %.2e)", regularization_rho))
    message(sprintf("  Final condition number: κ = %.2e", final_cond_num))
  }

  return(list(
    lambda_update = lambda_update_full,
    constraints_dropped = dropped_indices,
    n_drops = length(dropped_indices),
    used_regularization = TRUE,
    regularization_rho = regularization_rho,
    condition_number = final_cond_num,
    converged = TRUE
  ))
}


#' Wrapper for newton_qdelta with adaptive constraints (STUB)
#'
#' Adds adaptive constraint dropping to the Newton-Raphson solver.
#' When enabled, will automatically drop redundant constraints when the
#' Jacobian becomes singular (e.g., at the boundary of the marginal polytope).
#'
#' NOTE: This is currently a stub function. When adaptive_constraints=FALSE,
#' it forwards to newton_qdelta. When adaptive_constraints=TRUE, it throws an
#' error. Future implementation will integrate adaptive constraint dropping.
#'
#' @param ... All arguments passed to newton_qdelta
#' @param adaptive_constraints Use adaptive constraint dropping (default FALSE)
#' @param max_constraint_drops Maximum constraints to drop per iteration (default 2)
#' @param condition_threshold Condition number threshold for singularity (default 1e10)
#' @param verbose Print diagnostic messages (default FALSE)
#' @return Same as newton_qdelta, plus:
#'   - constraints_dropped_history: List tracking which constraints were dropped each iteration
#'   - total_constraints_dropped: Total number of unique constraints dropped
#' @keywords internal
newton_qdelta_adaptive <- function(..., adaptive_constraints = FALSE,
                                   max_constraint_drops = 2,
                                   condition_threshold = 1e10,
                                   verbose = FALSE) {

  if (!adaptive_constraints) {
    # Standard Newton-Raphson without adaptive constraints
    # Forward to original function
    args <- list(...)
    return(do.call(newton_qdelta, args))
  }

  # TODO: Implement adaptive version that intercepts the Jacobian solve step
  # This requires modifying newton_qdelta or creating a parallel implementation
  # For now, return a message
  stop("Adaptive constraint dropping not yet fully implemented. Use adaptive_constraints=FALSE or modify newton_qdelta directly.")

  # Future implementation would:
  # 1. Call newton_qdelta with a modified solve step
  # 2. Replace line 147 solve(jacobian, constraint_vector) with adaptive_solve()
  # 3. Track which constraints were dropped across iterations
  # 4. Return extended result with dropping diagnostics
}
