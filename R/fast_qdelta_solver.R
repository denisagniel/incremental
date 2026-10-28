#' @file fast_qdelta_solver.R
#' @title Fast Per-Observation Q_delta Solver
#' @description
#' Newton-Raphson solver for computing Q_delta distributions via constrained
#' optimization. Solves the (m+1) x (m+1) system (internal normalization) or
#' m x m system (external normalization) with analytical Jacobians. Supports
#' adaptive constraint dropping for singular systems.
#'
#' @section Matrix form of the moment/Jacobian pass (read this before editing):
#'
#' The paper presents the constrained-tilt problem as a **sum over the \eqn{2^m}
#' treatment atoms**; this file implements the same quantities as **matrix
#' products**, and the two look nothing alike. The translation, once, so nobody
#' has to re-derive it:
#'
#' Let \eqn{B} be the active slice of the \eqn{2^m \times m} combination matrix
#' ([generate_binary_vectors()]) and \eqn{B_{ext} = [\mathbf{1}, B]} its extended
#' \eqn{n_a \times (m+1)} form. Under internal normalization the tilted weights
#' are \eqn{q_k = p_{0k}\, h(u_k)} with linear term
#' \eqn{u_k = \lambda_0 + \sum_j a_{kj}\lambda_j}, i.e.
#'
#' \deqn{u = B_{ext}\lambda,\qquad q = p_0 \odot h(u).}
#'
#' The \eqn{m+1} residuals are the \eqn{m} marginal constraints plus
#' normalization, i.e. \eqn{F = B_{ext}^\top q} reordered, minus the targets.
#' Differentiating, every Jacobian entry is a weighted second moment of
#' \eqn{B_{ext}}:
#'
#' \deqn{\frac{\partial F_i}{\partial \lambda_j}
#'       = \sum_k p_{0k} h'(u_k)\, r_i(k)\, B_{ext}[k,j],}
#'
#' where the row multiplier \eqn{r_i(k)} is \eqn{a_{ki} = B_{ext}[k, i+1]} for a
#' marginal row \eqn{i \le m} and \eqn{1 = B_{ext}[k, 1]} for the normalization
#' row \eqn{i = m+1}. So **the whole Jacobian is one symmetric scatter matrix
#' with its rows permuted**:
#'
#' \deqn{S_2 = B_{ext}^\top \,\mathrm{diag}(p_0 \odot h'(u))\, B_{ext},
#'       \qquad J = S_2[(2,\dots,m+1,1),\, \cdot\,].}
#'
#' That is `compute_moment_scatter_internal()` +
#' `compute_jacobian_analytical_internal()`. The row permutation
#' \eqn{(2,\dots,m+1,1)} is the *only* reason the code does not read as a plain
#' `crossprod`: it maps the scatter matrix's "intercept-first" column order onto
#' the residual vector's "marginals-first, normalization-last" row order used by
#' `compute_constraint_vector_internal()`.
#'
#' Two consequences exploited here:
#' \itemize{
#'   \item \eqn{S_2} is the same object for the marginal rows and the
#'     normalization row, so \eqn{h'} is evaluated **once per Newton iteration**
#'     over the whole active vector -- not once per \eqn{(i,j,k)} triple. Before
#'     2026-08-11 this file evaluated it \eqn{(m+1)^2} times per combination per
#'     iteration (1.77M scalar [divergence_h_pair()] calls per observation at
#'     m = 11), which profiling showed to be 99.8% of the solve. See
#'     `explorations/investigations/profile_qdelta_solver_pre_rewrite.md`.
#'   \item For `divergence = "chi2"`, \eqn{h'(x) \equiv 1/2} for all \eqn{x}, so
#'     \eqn{S_2} does not depend on \eqn{\lambda} at all and is hoisted out of
#'     the Newton loop entirely (computed once per solve).
#' }
#'
#' The pre-vectorization scalar-loop implementation is preserved verbatim in
#' `R/fast_qdelta_solver_legacy.R` as a numerical-equivalence reference; see
#' `tests/testthat/test-qdelta-vectorized-equivalence.R`.
#'
#' @author Denis Agniel
#
# Note: adaptive_constraints.R functions (adaptive_solve, is_matrix_singular, etc.)
# are available via package namespace. Load order managed by R's Collate field.

#' Newton-Raphson solver for per-observation Q_delta computation
#'
#' Uses the fixed internal normalization approach that matches the theoretical framework.
#' Solves the (m+1)×(m+1) system with proper constraint formulation.
#'
#' @param propensity_scores_i Vector of propensity scores p_0(a|X_i) for all 2^m combinations
#' @param delta Vector of target marginal changes (length m)
#' @param divergence Type of divergence ("kl", "reverse_kl", "chi2", "hellinger")
#' @param max_iter Maximum iterations (default 50)
#' @param tol Convergence tolerance (default 1e-8)
#' @param use_analytical_jacobian Whether to use analytical Jacobian (default TRUE)
#' @param fix_lambda0 Whether to use external normalization (default FALSE for internal)
#' @param adaptive_constraints Use adaptive constraint dropping for singular Jacobians (default FALSE)
#' @param max_constraint_drops Maximum constraints to drop per iteration (default 2)
#' @param condition_threshold Condition number threshold for singularity detection (default 1e10)
#' @param verbose_adaptive Print adaptive constraint diagnostics (default FALSE)
#' @param combination_matrix Optional `2^m x m` matrix of all binary treatment
#'   combinations, as returned by [generate_binary_vectors()]. It depends only on
#'   `m`, which is fixed for the whole solve, so callers that solve many
#'   observations at one `m` (e.g. `efficient_influence_function()`) should
#'   compute it once and pass it here. `NULL` (default) computes it internally,
#'   preserving backward compatibility for direct callers.
#' @return List with lambda solution, convergence info, and q_delta weights
newton_qdelta <- function(propensity_scores_i, delta, divergence,
                         max_iter = 200, tol = 1e-10, use_analytical_jacobian = TRUE,
                         fix_lambda0 = FALSE, adaptive_constraints = FALSE,
                         max_constraint_drops = 2, condition_threshold = 1e10,
                         verbose_adaptive = FALSE, combination_matrix = NULL) {

  # Input validation
  stopifnot(
    "propensity_scores_i must be numeric" = is.numeric(propensity_scores_i),
    "delta must be numeric" = is.numeric(delta),
    "propensity_scores_i must be a vector" = is.vector(propensity_scores_i),
    "delta must be a vector" = is.vector(delta),
    "propensity_scores_i length must equal 2^m" = length(propensity_scores_i) == 2^length(delta),
    "propensity scores must sum to approximately 1" = abs(sum(propensity_scores_i) - 1) < PROPENSITY_SUM_TOL,
    "divergence must be valid" = divergence %in% c("kl", "reverse_kl", "chi2", "hellinger"),
    "max_iter must be positive" = max_iter > 0,
    "tol must be positive" = tol > 0,
    "condition_threshold must be positive" = condition_threshold > 0,
    "max_constraint_drops must be non-negative" = max_constraint_drops >= 0
  )

  m <- length(delta)
  n_combinations <- length(propensity_scores_i)

  # Combination matrix reuse (NOT caching): it depends only on m, which is fixed
  # for the whole solve, so it is built at most once here and then threaded into
  # the rank check, the active slice and every Newton iteration. The
  # pre-2026-08-11 code regenerated it up to three times per observation (once
  # here, once inside handle_zero_propensities(), once inside
  # compute_target_marginals()).
  if (is.null(combination_matrix)) {
    combination_matrix <- generate_binary_vectors(m)
  } else {
    stopifnot(
      "combination_matrix must be a matrix" = is.matrix(combination_matrix),
      "combination_matrix must have 2^m rows" = nrow(combination_matrix) == n_combinations,
      "combination_matrix must have m columns" = ncol(combination_matrix) == m
    )
  }

  # Initialize adaptive constraint tracking
  constraints_dropped_history <- list()
  total_adaptive_drops <- 0
  adaptive_diagnostics <- list()

  # Handle zero/near-zero propensity scores
  # Now throws error if n_nonzeros < m + 1
  propensity_scores_i <- handle_zero_propensities(propensity_scores_i, divergence, m,
                                                 combination_matrix = combination_matrix)

  # Identify active subset (combinations with p0 > threshold)
  active_idx <- which(propensity_scores_i > PROPENSITY_ACTIVE_THRESHOLD)
  p0_active <- propensity_scores_i[active_idx]

  # Subset the (reused) combination matrix to the active support. b_ext_active is
  # the extended [1, B] design used by every matrix-form expression below; see
  # the "Matrix form" section in this file's header.
  treatment_combinations_all <- combination_matrix
  treatment_combinations <- treatment_combinations_all[active_idx, , drop = FALSE]
  b_ext_active <- cbind(1, treatment_combinations)

  # Target marginals (current + delta)
  target_marginals <- compute_target_marginals(propensity_scores_i, treatment_combinations_all, delta, m = m)

  # Closed-form fast path: a single eligible treatment (m = 1) fully determines
  # Q_delta from the marginal + normalization constraints, independent of the
  # divergence. No Newton iteration needed. Only the internal-normalization
  # convention is closed-form here (fix_lambda0 = FALSE, the only path in use).
  if (m == 1 && !fix_lambda0) {
    return(solve_qdelta_single_treatment(active_idx, p0_active, treatment_combinations,
                                         target_marginals, divergence, n_combinations))
  }

  if (fix_lambda0) {
    # External normalization approach (m×m system)
    if (divergence == "reverse_kl") {
      lambda0 <- REVERSE_KL_LAMBDA0_INIT
    } else {
      lambda0 <- KL_LAMBDA0_INIT
    }
    # Using active subset for initialization
    lambda_active <- get_initial_lambda_active(delta, p0_active, treatment_combinations, divergence)

    # Newton-Raphson iteration on m-dimensional system
    for (iter in 1:max_iter) {

      # Compute constraint vector F(λ) - m constraints only
      constraint_vector <- compute_constraint_vector_m(lambda_active, p0_active,
                                                    treatment_combinations, target_marginals, divergence, lambda0)

      # Check convergence
      if (sqrt(sum(constraint_vector^2)) < tol) {
        break
      }

      # Compute Jacobian (m×m)
      if (use_analytical_jacobian) {
        jacobian <- compute_jacobian_analytical_m(lambda_active, p0_active,
                                                treatment_combinations, divergence, lambda0)
      } else {
        jacobian <- compute_jacobian_numerical_m(lambda_active, p0_active,
                                               treatment_combinations, target_marginals, divergence, lambda0)
      }

      # Newton step: λ^(t+1) = λ^(t) - J^(-1) F(λ^(t))
      tryCatch({
        # Use adaptive constraint dropping if enabled
        if (adaptive_constraints) {
          adaptive_result <- adaptive_solve(
            jacobian, constraint_vector, lambda_active,
            treatment_vectors = treatment_combinations,
            propensity_scores = p0_active,
            max_drops = max_constraint_drops,
            cond_threshold = condition_threshold,
            verbose = verbose_adaptive
          )
          lambda_update <- adaptive_result$lambda_update

          # Track dropped constraints and regularization
          if (adaptive_result$n_drops > 0 || adaptive_result$used_regularization) {
            constraints_dropped_history[[length(constraints_dropped_history) + 1]] <- list(
              iteration = iter,
              constraints_dropped = adaptive_result$constraints_dropped,
              n_drops = adaptive_result$n_drops,
              used_regularization = adaptive_result$used_regularization,
              condition_number = adaptive_result$condition_number
            )
            total_adaptive_drops <- total_adaptive_drops + adaptive_result$n_drops
          }

          lambda_new <- lambda_active - lambda_update
        } else {
          # Standard solve
          lambda_new <- lambda_active - solve(jacobian, constraint_vector)
        }

        # Check for constraint violations (reverse_kl, chi2, hellinger)
        if (divergence %in% c("reverse_kl", "chi2", "hellinger")) {
          if (!check_constraint_satisfaction_m(lambda_new, p0_active,
                                             treatment_combinations, divergence, lambda0)) {
            # If constraints violated, use line search
            lambda_active <- line_search_newton_m(lambda_active, lambda_new, p0_active,
                                                treatment_combinations, target_marginals, divergence, lambda0)
          } else {
            lambda_active <- lambda_new
          }
        } else {
          lambda_active <- lambda_new
        }

      }, error = function(e) {
        # If Jacobian is singular and adaptive constraints not enabled, use gradient descent fallback
        if (verbose_adaptive) {
          message(sprintf("Jacobian solve failed at iteration %d: %s", iter, e$message))
          if (!adaptive_constraints) {
            message("Hint: Try adaptive_constraints=TRUE to handle singular Jacobians")
          }
        }
        lambda_active <<- lambda_active - 0.1 * constraint_vector
      })
    }

    # Final lambda with λ₀ = 0
    lambda <- c(lambda0, lambda_active)

  } else {
    # Internal normalization approach (m+1)×(m+1) system
    lambda <- get_initial_lambda_internal(delta, p0_active, treatment_combinations, divergence,
                                          b_ext_active = b_ext_active)

    # Row permutation mapping the symmetric moment-scatter matrix S2 (intercept
    # column first) onto the residual vector's row order (m marginals first,
    # normalization last). See the "Matrix form" section in this file's header.
    jacobian_row_perm <- c(seq_len(m) + 1L, 1L)

    # chi2 has h'(x) == 0.5 for EVERY x, so S2 does not depend on lambda: compute
    # it ONCE here rather than recomputing a constant matrix on every Newton
    # iteration.
    chi2_moment_scatter <- if (divergence == "chi2") {
      crossprod(b_ext_active, as.vector(p0_active * 0.5) * b_ext_active)
    } else {
      NULL
    }

    # Newton-Raphson iteration on (m+1)×(m+1) system
    prev_constraint_norm <- Inf
    use_levenberg_marquardt <- FALSE
    for (iter in 1:max_iter) {

      # Compute constraint vector F(λ) - (m+1) constraints
      constraint_vector <- compute_constraint_vector_internal(lambda, p0_active,
                                                           treatment_combinations, target_marginals, divergence)

      # Check convergence with adaptive tolerance
      constraint_norm <- sqrt(sum(constraint_vector^2))

      effective_tol <- tol

      # Handle NAs from numerical issues (can happen after regularization)
      if (is.na(constraint_norm)) {
        if (verbose_adaptive) {
          message(sprintf("Iteration %d: Numerical issues detected (constraint_norm is NA)", iter))
        }
        break  # Exit with current lambda
      }

      if (constraint_norm < effective_tol) {
        break
      }

      # Check for divergence
      if (iter > 5 && constraint_norm > prev_constraint_norm * 1.5) {
        # Switch to numerical Jacobian
        use_analytical_jacobian <- FALSE
      }


      # Compute Jacobian (m+1)×(m+1) with fallback
      tryCatch({
        if (use_analytical_jacobian) {
          jacobian <- compute_jacobian_analytical_internal(lambda, p0_active,
                                                         treatment_combinations, divergence,
                                                         b_ext_active = b_ext_active,
                                                         moment_scatter = chi2_moment_scatter,
                                                         row_perm = jacobian_row_perm)
        } else {
          jacobian <- compute_jacobian_numerical_internal(lambda, p0_active,
                                                         treatment_combinations, target_marginals, divergence)
        }
      }, error = function(e) {
        # Fallback to numerical Jacobian
        jacobian <<- compute_jacobian_numerical_internal(lambda, p0_active,
                                                       treatment_combinations, target_marginals, divergence)
      })

      # Newton step with damping
      tryCatch({
        # Use adaptive constraint dropping if enabled
        if (adaptive_constraints) {
          adaptive_result <- adaptive_solve(
            jacobian, constraint_vector, lambda,
            treatment_vectors = treatment_combinations,
            propensity_scores = p0_active,
            max_drops = max_constraint_drops,
            cond_threshold = condition_threshold,
            verbose = verbose_adaptive
          )
          lambda_update <- adaptive_result$lambda_update

          # Track dropped constraints and regularization
          if (adaptive_result$n_drops > 0 || adaptive_result$used_regularization) {
            constraints_dropped_history[[length(constraints_dropped_history) + 1]] <- list(
              iteration = iter,
              constraints_dropped = adaptive_result$constraints_dropped,
              n_drops = adaptive_result$n_drops,
              used_regularization = adaptive_result$used_regularization,
              condition_number = adaptive_result$condition_number,
              converged = adaptive_result$converged
            )
            total_adaptive_drops <- total_adaptive_drops + adaptive_result$n_drops

            if (verbose_adaptive) {
              reg_msg <- if (adaptive_result$used_regularization) " (regularized)" else ""
              message(sprintf("Iteration %d: Dropped %d constraint(s), κ = %.2e%s",
                              iter, adaptive_result$n_drops, adaptive_result$condition_number, reg_msg))
            }
          }

          lambda_new <- lambda - lambda_update
        } else if (divergence == "chi2" && use_levenberg_marquardt) {
          # Chi2-specific: Use gradient descent with backtracking line search
          # Since h''(x) = 0, the problem is essentially linear, making
          # Newton-Raphson ineffective. Gradient descent is more appropriate.

          # Gradient direction: -J^T F (steepest descent)
          gradient <- -t(jacobian) %*% constraint_vector

          # Normalize gradient
          grad_norm <- sqrt(sum(gradient^2))
          if (grad_norm > 0) {
            gradient <- gradient / grad_norm
          }

          # Backtracking line search with domain-aware step size
          alpha <- 1.0
          best_norm <- constraint_norm
          best_lambda <- lambda

          for (ls_iter in 1:10) {
            lambda_candidate <- lambda + alpha * gradient

            # Check domain constraints
            if (check_constraint_satisfaction_internal(lambda_candidate, p0_active,
                                                     treatment_combinations, divergence,
                                                     b_ext_active = b_ext_active)) {
              # Evaluate constraint at candidate
              F_candidate <- compute_constraint_vector_internal(lambda_candidate, p0_active,
                                                               treatment_combinations, target_marginals, divergence)
              norm_candidate <- sqrt(sum(F_candidate^2))

              # Accept if improvement
              if (norm_candidate < best_norm) {
                best_norm <- norm_candidate
                best_lambda <- lambda_candidate
              }
            }

            # Reduce step size
            alpha <- alpha * 0.5
          }

          lambda_new <- best_lambda

          # Diagnostic output
          if (verbose_adaptive && iter %% 20 == 0) {
            message(sprintf("  Iter %d: Gradient descent, constraint norm %.3e → %.3e",
                          iter, constraint_norm, best_norm))
          }
        } else {
          # Standard solve
          lambda_new <- lambda - solve(jacobian, constraint_vector)
        }

        # Check for constraint violations and apply line search if needed
        if (divergence %in% c("reverse_kl", "chi2", "hellinger")) {
          constraint_check <- tryCatch({
            check_constraint_satisfaction_internal(lambda_new, p0_active,
                                                  treatment_combinations, divergence,
                                                  b_ext_active = b_ext_active)
          }, error = function(e) {
            FALSE  # If check fails, assume constraints violated
          })

          if (isTRUE(constraint_check)) {
            lambda <- lambda_new
          } else {
            # Line search
            lambda <- line_search_newton_internal(lambda, lambda_new, p0_active,
                                                treatment_combinations, target_marginals, divergence,
                                                b_ext_active = b_ext_active)
          }
        } else {
          lambda <- lambda_new
        }

      }, error = function(e) {
        # If Jacobian is singular and adaptive constraints not enabled, use gradient descent fallback
        if (verbose_adaptive) {
          message(sprintf("Jacobian solve failed at iteration %d: %s", iter, e$message))
          if (!adaptive_constraints) {
            message("Hint: Try adaptive_constraints=TRUE to handle singular Jacobians")
          }
        }
        lambda <<- lambda - 0.1 * constraint_vector
      })

      # Early termination for Chi2 boundary cases: if stalled near boundary, stop
      # NOTE: compares against prev_constraint_norm from the END of the
      # PREVIOUS iteration (updated below, after both stall-checks). Bugfix
      # 2026-08-10: this update used to run immediately above, before this
      # check, which made `constraint_norm >= prev_constraint_norm * 0.999`
      # tautologically true (comparing constraint_norm to itself) for every
      # chi2 solve past iteration 5, regardless of whether progress had
      # actually stalled. See MEMORY.md [LEARN] entry.
      if (divergence == "chi2" && iter > 5 && !is.na(constraint_norm) &&
          !is.na(prev_constraint_norm) && constraint_norm >= prev_constraint_norm * 0.999) {
        # Check if at domain boundary (matrix form of the per-row linear term)
        boundary_linear <- as.vector(b_ext_active %*% lambda)
        if (any(boundary_linear <= -2 + 1e-4)) {
          break  # At boundary, can't improve further
        }
      }

      # Progress monitoring for long runs
      # NOTE: same prev_constraint_norm-timing dependency as above -- this
      # check is only meaningful because the update below now happens after
      # it, not before.
      if (iter %% 25 == 0 && iter < max_iter) {
        # Check if making progress (with NA protection)
        if (iter > 50 && !is.na(constraint_norm) && !is.na(prev_constraint_norm) &&
            constraint_norm > prev_constraint_norm * 0.95) {
          # Stalled: try to escape by perturbing lambda slightly (deterministic)
          # Skip for chi2 (boundary stalling is expected)
          if (divergence != "chi2") {
            lambda <- lambda + seq(-NEWTON_PERTURBATION_MAGNITUDE,
                                   NEWTON_PERTURBATION_MAGNITUDE,
                                   length.out = length(lambda))
          }
        }
      }

      # Update AFTER the stall checks above so they compare against the
      # previous iteration's norm, not this iteration's own (see bugfix note
      # above).
      prev_constraint_norm <- constraint_norm
    }
  }

  # Compute final q_delta weights (internal normalization)
  # q_delta_weights must have original length (2^m)
  q_delta_weights_full <- numeric(n_combinations)

  if (fix_lambda0) {
    # External normalization: normalize weights
    q_delta_weights_active <- compute_q_delta_weights(lambda, p0_active,
                                             treatment_combinations, divergence)
    q_delta_weights_active <- q_delta_weights_active / sum(q_delta_weights_active)

    # Final constraint check
    final_constraints <- compute_constraint_vector_m(lambda[2:(m+1)], p0_active,
                                                  treatment_combinations, target_marginals, divergence, lambda[1])
  } else {
    # Internal normalization: weights are already normalized
    q_delta_weights_active <- compute_q_delta_weights_internal(lambda, p0_active,
                                                     treatment_combinations, divergence)

    # Final constraint check
    final_constraints <- compute_constraint_vector_internal(lambda, p0_active,
                                                          treatment_combinations, target_marginals, divergence)
  }

  # Map active weights back to full vector
  q_delta_weights_full[active_idx] <- q_delta_weights_active

  # Check if convergence failed (with NA protection)
  final_constraint_norm <- sqrt(sum(final_constraints^2))
  final_tol <- tol
  final_converged <- !is.na(final_constraint_norm) && (final_constraint_norm < final_tol)

  # Chi2 boundary convergence: accept solution when at domain boundary

  # For Chi2, the solution may lie exactly on the boundary (linear_term = -2),
  # which is the best achievable. Accept if constraint norm is small.
  if (!final_converged && divergence == "chi2" && !fix_lambda0) {
    test_linear <- as.vector(b_ext_active %*% lambda)
    if (any(test_linear <= -2 + 1e-6) && !is.na(final_constraint_norm) &&
        final_constraint_norm < 0.1) {
      # At domain boundary — this is the best Chi2 can achieve
      final_converged <- TRUE
    }
  }

  if (!final_converged && iter >= max_iter) {
    norm_msg <- if (is.na(final_constraint_norm)) "NA (numerical issues)" else sprintf("%.2e", final_constraint_norm)

    # Chi2-specific diagnostics
    if (divergence == "chi2" && !fix_lambda0) {
      if (!exists("test_linear")) {
        test_linear <- as.vector(b_ext_active %*% lambda)
      }

      if (any(test_linear <= -1.9)) {
        stop(sprintf(
          paste0(
            "Chi2 convergence failure: Domain boundary issue.\n",
            "Linear terms approaching -2 boundary (min = %.3f).\n",
            "Constraint norm: %s after %d iterations.\n",
            "Suggestions:\n",
            "  1. Try smaller delta values (current: %s)\n",
            "  2. Check feasibility with analyze_delta_feasibility()\n",
            "  3. Consider using 'kl' divergence instead (more robust)"
          ),
          min(test_linear), norm_msg, iter, paste(delta, collapse=", ")),
          call. = FALSE
        )
      } else {
        stop(sprintf(
          paste0(
            "Chi2 convergence failure.\n",
            "Constraint norm: %s after %d iterations.\n",
            "Suggestions:\n",
            "  1. This may indicate infeasible delta or numerical instability\n",
            "  2. Try smaller delta values (current: %s)\n",
            "  3. Consider using 'kl' or 'reverse_kl' divergence (more robust)\n",
            "  4. Try adaptive_constraints=TRUE"
          ),
          norm_msg, iter, paste(delta, collapse=", ")),
          call. = FALSE
        )
      }
    }

    # Standard error message for other divergences
    stop(sprintf(
      "Newton-Raphson solver failed to converge for %s divergence after %d iterations (constraint norm = %s).\n",
      divergence, iter, norm_msg),
      "Possible solutions:\n",
      "  1. Check if delta is feasible for your propensity scores (use analyze_delta_feasibility)\n",
      "  2. Try a smaller delta value (closer to observed marginals)\n",
      "  3. Increase max_iter (currently ", max_iter, ")\n",
      "  4. Increase tolerance (currently ", tol, ")\n",
      "  5. Try a different divergence measure (e.g., 'kl' instead of '", divergence, "')\n",
      "  6. Check for severe positivity violations in your propensity scores"
    )
  }

  # Prepare adaptive constraint diagnostics
  if (adaptive_constraints && length(constraints_dropped_history) > 0) {
    adaptive_diagnostics <- list(
      constraints_dropped_history = constraints_dropped_history,
      total_constraints_dropped = total_adaptive_drops,
      n_iterations_with_drops = length(constraints_dropped_history)
    )
  } else {
    adaptive_diagnostics <- NULL
  }

  return(list(
    lambda = as.vector(lambda),
    q_delta_weights = q_delta_weights_full,
    iterations = iter,
    converged = final_converged,
    constraint_norm = final_constraint_norm,
    adaptive_diagnostics = adaptive_diagnostics
  ))
}

# =============================================================================
# INTERNAL NORMALIZATION FUNCTIONS (Support-Aware)
# =============================================================================

#' Better initialization for internal normalization
#'
#' @param delta Numeric vector, odds scaling factors (length m)
#' @param p0_active Numeric vector, active propensity scores
#' @param vectors_active Matrix, active treatment combinations (n_active x m)
#' @param divergence Character, divergence type
#' @param b_ext_active Optional pre-built extended design `cbind(1, vectors_active)`.
#'   Supplied by `newton_qdelta()`; `NULL` builds it locally.
#' @return Numeric vector of length m+1, initial lambda
#' @keywords internal
get_initial_lambda_internal <- function(delta, p0_active, vectors_active, divergence,
                                       b_ext_active = NULL) {

  m <- length(delta)

  # kl needs no domain feasibility check, so skip building the design entirely.
  if (divergence == "kl") {
    # For KL: start with small values, λ₀ can be 0
    return(c(0, delta * 0.1))
  }

  if (is.null(b_ext_active)) b_ext_active <- cbind(1, vectors_active)

  if (divergence == "reverse_kl") {
    # More conservative initialization for reverse_kl
    lambda <- c(-1.0, -abs(delta) * 0.2)

    # Verify all linear terms are safely negative (matrix form: u = B_ext lambda)
    test_linear_terms <- as.vector(b_ext_active %*% lambda)

    # Adjust to ensure max is at most -0.1
    if (max(test_linear_terms) >= -0.1) {
      lambda[1] <- lambda[1] + min(test_linear_terms) - 0.2
    }

  } else if (divergence == "chi2") {
    # For chi2: use conservative initialization
    # Since h(x) = x/2 + 1, the domain is linear_term > -2
    # Start with very small values to avoid boundary issues

    # Conservative approach: start near zero
    # Chi2's constant Jacobian means we don't benefit much from clever initialization
    lambda_k <- delta * 0.1

    # Set intercept to ensure all linear terms are safely above -2
    # With binary treatments, worst case is when all treatments are 1
    # So lambda[1] + sum(lambda[2:(m+1)]) > -2
    # Start conservatively at -0.5 to have plenty of room
    lambda_0 <- -0.5

    lambda <- c(lambda_0, lambda_k)

    # Verify domain safety: all linear terms must be > -2
    test_linear <- as.vector(b_ext_active %*% lambda)

    # Additional safety check: ensure the closest linear term stays clear of the
    # boundary so Newton steps have room
    if (min(test_linear) <= -1.5) {
      # Too close to boundary, scale everything down
      lambda <- lambda * 0.3
      lambda[1] <- -0.3
    }

  } else if (divergence == "hellinger") {
    # For hellinger: need all linear terms < 1
    lambda <- c(-0.1, delta * 0.05)

    # Verify constraints
    test_linear_terms <- as.vector(b_ext_active %*% lambda)

    if (any(test_linear_terms >= 1 - DOMAIN_VALIDATION_STRICT_TOL)) {
      lambda[1] <- -0.5  # More negative
    }
  }

  return(lambda)
}

#' Compute constraint vector for internal normalization (m+1 constraints)
compute_constraint_vector_internal <- function(lambda, p0_active, vectors_active,
                                              target_marginals, divergence) {

  m <- length(lambda) - 1

  # Compute q_delta weights (internal normalization - only for active set)
  q_delta_weights <- compute_q_delta_weights_internal(lambda, p0_active,
                                                    vectors_active, divergence)

  # Compute constraint vector (m+1 constraints) - vectorized
  # Marginal constraints (first m components)
  expected_A <- as.numeric(t(q_delta_weights) %*% vectors_active)
  constraints <- c(expected_A - target_marginals, sum(q_delta_weights) - 1)

  return(constraints)
}

#' Compute q_delta weights for internal normalization
compute_q_delta_weights_internal <- function(lambda, p0_active, vectors_active, divergence) {

  m <- length(lambda) - 1
  n_active <- length(p0_active)

  # Vectorized: compute linear terms for all combinations at once
  linear_terms <- lambda[1] + vectors_active %*% lambda[2:(m+1)]

  # Canonical h-function (see divergence_h_functions.R) -- clips once, shared
  # with the Jacobian below so h and h' are evaluated at the same argument.
  h_values <- divergence_h_pair(linear_terms, divergence)$h

  weights <- p0_active * h_values

  return(weights)
}

#' Weighted second-moment scatter of the extended design (internal normalization)
#'
#' @description
#' Computes \eqn{S_2 = B_{ext}^\top \mathrm{diag}(p_0 \odot h'(u)) B_{ext}} with
#' \eqn{u = B_{ext}\lambda}, the single \eqn{(m+1) \times (m+1)} object from which
#' every Jacobian entry is read off (see the "Matrix form" section in this file's
#' header for the derivation). `divergence_h_pair()` is called ONCE on the full
#' length-\eqn{n_a} linear-term vector, replacing the pre-2026-08-11 pattern of
#' one scalar call per \eqn{(i, j, k)} triple.
#'
#' @param lambda Numeric vector of length m+1
#' @param p0_active Numeric vector, active propensity scores
#' @param b_ext_active Matrix, extended design `cbind(1, vectors_active)`
#' @param divergence Character, divergence type
#' @return Symmetric (m+1) x (m+1) numeric matrix
#' @keywords internal
compute_moment_scatter_internal <- function(lambda, p0_active, b_ext_active, divergence) {

  # u = B_ext lambda: all n_active linear terms in one BLAS call
  linear_terms <- as.vector(b_ext_active %*% lambda)

  # Canonical h'(x) (see divergence_h_functions.R), vectorized over the whole
  # linear-term vector. Same clipped argument as compute_q_delta_weights_internal()
  # uses for h, so J remains a valid Jacobian of the residual actually computed.
  h_prime <- divergence_h_pair(linear_terms, divergence)$h_prime

  weights <- as.vector(p0_active * h_prime)

  crossprod(b_ext_active, weights * b_ext_active)
}

#' Compute analytical Jacobian for internal normalization (m+1)×(m+1)
#'
#' @description
#' Matrix-form assembly: the full Jacobian is the weighted second-moment scatter
#' `compute_moment_scatter_internal()` with its rows permuted from
#' intercept-first to marginals-first order. See the "Matrix form" section in this
#' file's header for why this is equivalent to the paper's sum over \eqn{2^m}
#' atoms.
#'
#' @param lambda Numeric vector of length m+1
#' @param p0_active Numeric vector, active propensity scores
#' @param vectors_active Matrix, active treatment combinations (n_active x m)
#' @param divergence Character, divergence type
#' @param b_ext_active Optional pre-built extended design `cbind(1, vectors_active)`
#' @param moment_scatter Optional pre-computed \eqn{S_2}. Supplied by
#'   `newton_qdelta()` for `divergence = "chi2"`, whose \eqn{h' \equiv 1/2} makes
#'   \eqn{S_2} independent of `lambda` and therefore hoistable out of the Newton
#'   loop. Must be `NULL` for every other divergence.
#' @param row_perm Optional pre-computed row permutation `c(2:(m+1), 1)`
#' @return (m+1) x (m+1) numeric matrix
#' @keywords internal
compute_jacobian_analytical_internal <- function(lambda, p0_active, vectors_active, divergence,
                                                b_ext_active = NULL, moment_scatter = NULL,
                                                row_perm = NULL) {

  m <- length(lambda) - 1

  if (is.null(moment_scatter)) {
    if (is.null(b_ext_active)) b_ext_active <- cbind(1, vectors_active)
    moment_scatter <- compute_moment_scatter_internal(lambda, p0_active,
                                                     b_ext_active, divergence)
  }

  if (is.null(row_perm)) row_perm <- c(seq_len(m) + 1L, 1L)

  # Rows 1..m are the marginal constraints (row multiplier a_k[i] = B_ext[, i+1]);
  # row m+1 is normalization (row multiplier 1 = B_ext[, 1]).
  moment_scatter[row_perm, , drop = FALSE]
}

#' Compute a single marginal-constraint Jacobian entry for internal normalization
#'
#' @description
#' Single-cell accessor for \eqn{\partial E_Q[A_i] / \partial \lambda_j}, kept for
#' diagnostic use and API stability. The Newton loop does NOT call this: it
#' assembles the whole Jacobian at once via
#' `compute_jacobian_analytical_internal()`, since computing all \eqn{(m+1)^2}
#' entries cell-by-cell re-evaluates the same \eqn{2^m} linear terms
#' \eqn{(m+1)^2} times (the pre-2026-08-11 bottleneck). This entry equals
#' `compute_moment_scatter_internal(...)[i + 1, j]`.
#'
#' @param i Integer, marginal constraint index in 1..m
#' @param j Integer, lambda index in 1..(m+1) (1 = intercept)
#' @param lambda Numeric vector of length m+1
#' @param p0_active Numeric vector, active propensity scores
#' @param vectors_active Matrix, active treatment combinations
#' @param divergence Character, divergence type
#' @return Numeric scalar
#' @keywords internal
compute_marginal_jacobian_internal <- function(i, j, lambda, p0_active, vectors_active, divergence) {

  b_ext_active <- cbind(1, vectors_active)

  # Canonical h' (see divergence_h_functions.R) -- one vectorized call over all
  # active combinations, evaluated at the SAME clipped linear terms the weights
  # function uses.
  h_prime <- divergence_h_pair(as.vector(b_ext_active %*% lambda), divergence)$h_prime

  # ∂E_Q[A_i]/∂λ_j = Σ_k p_0(a_k|X) * a_k[i] * h'(u_k) * B_ext[k, j]
  sum(as.vector(p0_active) * vectors_active[, i] * as.vector(h_prime) * b_ext_active[, j])
}

#' Compute the normalization-constraint Jacobian row entry for internal normalization
#'
#' @description
#' Single-cell accessor for \eqn{\partial (\sum_k q_k) / \partial \lambda_j}; see
#' `compute_marginal_jacobian_internal()` for why the Newton loop does not use it.
#' This entry equals `compute_moment_scatter_internal(...)[1, j]`.
#'
#' @param j Integer, lambda index in 1..(m+1) (1 = intercept)
#' @param lambda Numeric vector of length m+1
#' @param p0_active Numeric vector, active propensity scores
#' @param vectors_active Matrix, active treatment combinations
#' @param divergence Character, divergence type
#' @return Numeric scalar
#' @keywords internal
compute_normalization_jacobian_internal <- function(j, lambda, p0_active, vectors_active, divergence) {

  b_ext_active <- cbind(1, vectors_active)

  h_prime <- divergence_h_pair(as.vector(b_ext_active %*% lambda), divergence)$h_prime

  # ∂(Σ q_delta)/∂λ_j = Σ_k p_0(a_k|X) * h'(u_k) * B_ext[k, j]
  sum(as.vector(p0_active) * as.vector(h_prime) * b_ext_active[, j])
}

#' Check constraint satisfaction for internal normalization
#'
#' @param lambda Numeric vector of length m+1
#' @param p0_active Numeric vector, active propensity scores (unused; kept for
#'   signature stability)
#' @param vectors_active Matrix, active treatment combinations
#' @param divergence Character, divergence type
#' @param b_ext_active Optional pre-built extended design `cbind(1, vectors_active)`
#' @return `TRUE` if every active linear term is strictly inside the divergence's
#'   domain (with the numerical buffer), `FALSE` otherwise. Divergences with no
#'   domain restriction (`kl`) always return `TRUE`.
#' @keywords internal
check_constraint_satisfaction_internal <- function(lambda, p0_active, vectors_active, divergence,
                                                 b_ext_active = NULL) {

  if (!(divergence %in% c("reverse_kl", "chi2", "hellinger"))) {
    return(TRUE)
  }

  if (is.null(b_ext_active)) b_ext_active <- cbind(1, vectors_active)

  linear_terms <- as.vector(b_ext_active %*% lambda)

  # all(x < c) is the vectorized complement of the scalar loop's
  # "return FALSE on the first x >= c".
  if (divergence == "reverse_kl") {
    all(linear_terms < -DOMAIN_VALIDATION_STRICT_TOL)
  } else if (divergence == "chi2") {
    # For chi2, domain is linear_term > -2 (so h(x) = x/2 + 1 > 0)
    all(linear_terms > CHI2_BOUNDARY_LAMBDA0 + CHI2_DOMAIN_BUFFER)
  } else {
    all(linear_terms < 1 - DOMAIN_VALIDATION_STRICT_TOL)
  }
}

#' Line search for internal normalization with Armijo condition
#'
#' Enhanced line search with:
#' - More steps (20 instead of 10)
#' - Armijo condition for sufficient decrease
#' - Proper objective function evaluation
#'
#' @param lambda_old Current lambda
#' @param lambda_new Proposed lambda (from Newton step)
#' @param p0_active Active propensity scores
#' @param vectors_active Active treatment combinations
#' @param target_marginals Target marginal probabilities
#' @param divergence Divergence type
#' @param max_steps Maximum line search steps (default 20)
#' @param c1 Armijo condition parameter (default 1e-4)
#' @param b_ext_active Optional pre-built extended design `cbind(1, vectors_active)`
#' @return Accepted lambda (either candidate or old)
line_search_newton_internal <- function(lambda_old, lambda_new, p0_active,
                                      vectors_active, target_marginals, divergence,
                                      max_steps = LINE_SEARCH_MAX_STEPS,
                                      c1 = LINE_SEARCH_ARMIJO_C1,
                                      b_ext_active = NULL) {

  if (is.null(b_ext_active)) b_ext_active <- cbind(1, vectors_active)

  # Current objective (constraint norm squared)
  constraint_old <- compute_constraint_vector_internal(
    lambda_old, p0_active, vectors_active, target_marginals, divergence
  )
  objective_old <- sum(constraint_old^2)

  direction <- lambda_new - lambda_old

  # For Chi2, compute max safe step size analytically
  # Domain: linear_term > -2 for all combinations
  # linear_term_j(lambda_old + alpha * direction) = lt_old_j + alpha * lt_dir_j
  # Need: lt_old_j + alpha * lt_dir_j > -2
  # For lt_dir_j < 0: alpha < (lt_old_j + 2) / (-lt_dir_j)
  alpha_max <- 1.0
  if (divergence == "chi2") {
    lt_old <- as.vector(b_ext_active %*% lambda_old)
    lt_dir <- as.vector(b_ext_active %*% direction)
    descending <- lt_dir < 0
    if (any(descending)) {
      # 0.1% buffer from boundary
      safe_alpha <- (lt_old[descending] + 2) / (-lt_dir[descending]) * 0.999
      alpha_max <- min(alpha_max, min(safe_alpha))
    }
  }

  # Try full step first (alpha_max), then backtrack
  alphas <- c(alpha_max, alpha_max / 2^(1:max_steps))

  for (alpha in alphas) {
    lambda_candidate <- lambda_old + alpha * direction

    # Check domain feasibility
    if (!check_constraint_satisfaction_internal(
      lambda_candidate, p0_active, vectors_active, divergence,
      b_ext_active = b_ext_active
    )) {
      next
    }

    # Armijo condition: sufficient decrease
    constraint_candidate <- compute_constraint_vector_internal(
      lambda_candidate, p0_active, vectors_active, target_marginals, divergence
    )
    objective_candidate <- sum(constraint_candidate^2)

    if (objective_candidate <= objective_old * (1 - c1 * alpha)) {
      return(lambda_candidate)
    }
  }

  return(lambda_old)  # No improvement found
}

#' Numerical Jacobian for internal normalization (fallback)
compute_jacobian_numerical_internal <- function(lambda, p0_active, vectors_active,
                                              target_marginals, divergence, h = 1e-6) {

  m <- length(lambda) - 1
  jacobian <- matrix(0, nrow = m + 1, ncol = m + 1)

  for (j in 1:(m + 1)) {
    lambda_plus <- lambda
    lambda_plus[j] <- lambda_plus[j] + h

    lambda_minus <- lambda
    lambda_minus[j] <- lambda_minus[j] - h

    constraints_plus <- compute_constraint_vector_internal(lambda_plus, p0_active, vectors_active,
                                                         target_marginals, divergence)
    constraints_minus <- compute_constraint_vector_internal(lambda_minus, p0_active, vectors_active,
                                                          target_marginals, divergence)

    jacobian[, j] <- (constraints_plus - constraints_minus) / (2 * h)
  }

  return(jacobian)
}

# =============================================================================
# EXTERNAL NORMALIZATION FUNCTIONS (m-dimensional system)
# =============================================================================

#' Compute constraint vector for m-dimensional system
compute_constraint_vector_m <- function(lambda_active, p0_active, vectors_active,
                                       target_marginals, divergence, lambda0 = 0) {

  m <- length(lambda_active)

  # Compute q_delta weights with fixed λ₀
  lambda <- c(lambda0, lambda_active)
  q_delta_weights <- compute_q_delta_weights(lambda, p0_active,
                                           vectors_active, divergence)

  # Normalize
  q_delta_normalized <- q_delta_weights / sum(q_delta_weights)

  # Compute constraint vector (m constraints only)
  constraints <- numeric(m)

  # Marginal constraints (vectorized)
  expected_A <- as.numeric(q_delta_normalized %*% vectors_active)
  constraints <- expected_A - target_marginals

  return(constraints)
}

#' Compute q_delta weights
compute_q_delta_weights <- function(lambda, p0_active, vectors_active, divergence) {

  m <- length(lambda) - 1
  n_active <- length(p0_active)

  # Vectorized: compute linear terms for all combinations at once
  linear_terms <- lambda[1] + vectors_active %*% lambda[2:(m+1)]

  # Apply h-function based on divergence (vectorized)
  # Canonical h-function (see divergence_h_functions.R).
  h_values <- divergence_h_pair(linear_terms, divergence)$h

  weights <- p0_active * h_values

  return(weights)
}

#' Compute analytical Jacobian for m-dimensional system
#'
#' @description
#' Matrix form of the external-normalization (`fix_lambda0 = TRUE`) Jacobian.
#' With \eqn{u = \lambda_0 + B\lambda_a}, \eqn{w = p_0 \odot h'(u)},
#' \eqn{Z = \sum_k q_k} and \eqn{\bar q = q / Z}:
#'
#' \deqn{J = \frac{B^\top \mathrm{diag}(w) B}{Z}
#'          - \frac{(\bar q^\top B)^\top (w^\top B)}{Z},}
#'
#' which is exactly the per-cell formula
#' \eqn{J_{ij} = \mathrm{sum1}_{ij}/Z - E_Q[A_i]\,\mathrm{sum2}_j/Z} of the
#' pre-2026-08-11 triple loop, with \eqn{h'} evaluated once over the whole active
#' vector instead of once per \eqn{(i, j, k)} triple.
#'
#' @param lambda_active Numeric vector of length m
#' @param p0_active Numeric vector, active propensity scores
#' @param vectors_active Matrix, active treatment combinations (n_active x m)
#' @param divergence Character, divergence type
#' @param lambda0 Numeric scalar, fixed intercept
#' @return m x m numeric matrix
#' @keywords internal
compute_jacobian_analytical_m <- function(lambda_active, p0_active, vectors_active, divergence, lambda0) {
  lambda <- c(lambda0, lambda_active)

  # Compute q_delta weights and normalization
  q_delta_weights <- compute_q_delta_weights(lambda, p0_active, vectors_active, divergence)
  Z <- sum(q_delta_weights)
  q_delta_normalized <- as.vector(q_delta_weights) / Z

  # Canonical h' (see divergence_h_functions.R), one vectorized call
  linear_terms <- as.vector(lambda0 + vectors_active %*% lambda_active)
  weights <- as.vector(p0_active) * as.vector(divergence_h_pair(linear_terms, divergence)$h_prime)

  sum1 <- crossprod(vectors_active, weights * vectors_active)         # m x m
  sum2 <- as.vector(crossprod(vectors_active, weights))               # length m
  expected_A <- as.vector(crossprod(q_delta_normalized, vectors_active))  # length m

  sum1 / Z - outer(expected_A, sum2) / Z
}

#' Compute numerical Jacobian for m-dimensional system
compute_jacobian_numerical_m <- function(lambda_active, p0_active, vectors_active, target_marginals, divergence, lambda0, h = 1e-6) {
  m <- length(lambda_active)
  jacobian <- matrix(0, nrow = m, ncol = m)
  for (j in 1:m) {
    lambda_plus <- lambda_active; lambda_plus[j] <- lambda_plus[j] + h
    lambda_minus <- lambda_active; lambda_minus[j] <- lambda_minus[j] - h
    f_plus <- compute_constraint_vector_m(lambda_plus, p0_active, vectors_active, target_marginals, divergence, lambda0)
    f_minus <- compute_constraint_vector_m(lambda_minus, p0_active, vectors_active, target_marginals, divergence, lambda0)
    jacobian[, j] <- (f_plus - f_minus) / (2 * h)
  }
  return(jacobian)
}

#' Get initial lambda values
get_initial_lambda <- function(delta, divergence) {

  m <- length(delta)

  if (divergence == "kl") {
    lambda <- c(0, delta * 0.1)
  } else if (divergence == "reverse_kl") {
    lambda <- c(-1, -abs(delta) * 0.5)
  } else if (divergence == "chi2") {
    lambda <- c(0, rep(0, m))
  } else if (divergence == "hellinger") {
    lambda <- c(-0.5, delta * 0.1)
  } else {
    lambda <- c(0, rep(0, m))
  }

  return(lambda)
}

#' Check constraint satisfaction for m-dimensional system
#'
#' @param lambda_active Numeric vector of length m
#' @param p0_active Numeric vector, active propensity scores (unused; kept for
#'   signature stability)
#' @param vectors_active Matrix, active treatment combinations
#' @param divergence Character, divergence type
#' @param lambda0 Numeric scalar, fixed intercept
#' @return Logical, `TRUE` if all active linear terms are inside the domain
#' @keywords internal
check_constraint_satisfaction_m <- function(lambda_active, p0_active, vectors_active, divergence, lambda0 = 0) {
  if (!(divergence %in% c("reverse_kl", "chi2", "hellinger"))) {
    return(TRUE)
  }
  linear_terms <- as.vector(lambda0 + vectors_active %*% lambda_active)
  if (divergence == "reverse_kl") {
    all(linear_terms < -NEAR_ZERO_THRESHOLD)
  } else if (divergence == "chi2") {
    all(linear_terms > CHI2_BOUNDARY_LAMBDA0)
  } else {
    all(linear_terms < 1)
  }
}

#' Line search for m-dimensional system
line_search_newton_m <- function(lambda_old, lambda_new, p0_active,
                                vectors_active, target_marginals, divergence, lambda0,
                                max_steps = 10) {

  for (step in 1:max_steps) {
    alpha <- 1 / (2^step)
    lambda_candidate <- (1 - alpha) * lambda_old + alpha * lambda_new

    if (check_constraint_satisfaction_m(lambda_candidate, p0_active,
                                      vectors_active, divergence, lambda0)) {
      return(lambda_candidate)
    }
  }

  # If line search fails, return old lambda
  return(lambda_old)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Closed-form Q_delta for a single eligible treatment (m = 1)
#'
#' With a single binary treatment, Q_delta is a two-point distribution over
#' \{a = 0, a = 1\}. The marginal constraint \eqn{E_Q[A] = m^*} (the target
#' marginal) together with normalization \eqn{q_0 + q_1 = 1} pins it down
#' exactly: \eqn{q_1 = m^*} and \eqn{q_0 = 1 - m^*}. This solution does not
#' depend on the divergence, so no Newton iteration is required.
#'
#' Downstream EIF code re-derives the h-values from \code{lambda}, so we also
#' recover the internal-normalization lambda by inverting the divergence's
#' h-function at the required weight ratios \eqn{w_a = q_a / p_0(a)}:
#' \eqn{\lambda_0 = h^{-1}(w_0)} and \eqn{\lambda_1 = h^{-1}(w_1) - \lambda_0}.
#' Because both weight ratios are strictly positive, every inverse stays in the
#' valid domain of its h-function.
#'
#' @param active_idx Integer vector, indices (into the length-2 combination
#'   vector) of combinations with non-negligible propensity. Both are active
#'   here: \code{handle_zero_propensities()} already errored out otherwise
#'   (rank of the m = 1 design requires both combinations).
#' @param p0_active Numeric vector, active propensity scores (length 2).
#' @param treatment_combinations Matrix, active treatment combinations
#'   (2 x 1, rows \code{0} and \code{1}).
#' @param target_marginals Numeric scalar (length 1), target marginal
#'   \eqn{Q_\delta(A = 1)}.
#' @param divergence Character, divergence type.
#' @param n_combinations Integer, length of the full weight vector (2 for m = 1).
#' @return List matching \code{newton_qdelta()}'s return contract.
solve_qdelta_single_treatment <- function(active_idx, p0_active, treatment_combinations,
                                           target_marginals, divergence, n_combinations) {

  # Row a = 1 in the active combination matrix (the treated point).
  treated_row <- which(treatment_combinations[, 1] == 1)
  control_row <- which(treatment_combinations[, 1] == 0)

  # Two-point Q_delta: q_1 = target marginal, q_0 = 1 - target marginal.
  q_active <- numeric(length(p0_active))
  q_active[treated_row] <- target_marginals[1]
  q_active[control_row] <- 1 - target_marginals[1]

  # Weight ratios w_a = q_a / p_0(a); strictly positive since both active.
  w <- q_active / p0_active

  # Invert the divergence h-function to recover lambda = (lambda_0, lambda_1)
  # under internal normalization, matching compute_q_delta_weights_internal().
  h_inverse <- function(weight) {
    if (divergence == "kl") {
      log(weight)                    # h(x) = exp(x)
    } else if (divergence == "reverse_kl") {
      -1 / weight                    # h(x) = -1/x  (x < 0, weight > 0)
    } else if (divergence == "chi2") {
      2 * (weight - 1)               # h(x) = x/2 + 1
    } else if (divergence == "hellinger") {
      1 - 1 / sqrt(weight)           # h(x) = 1/(1-x)^2  (x < 1)
    } else {
      stop("Unsupported divergence in closed-form single-treatment solver: ", divergence)
    }
  }

  x_control <- h_inverse(w[control_row])   # linear term at a = 0  => lambda_0
  x_treated <- h_inverse(w[treated_row])   # linear term at a = 1  => lambda_0 + lambda_1
  lambda <- c(x_control, x_treated - x_control)

  # Map active weights back to the full-length (2^m = 2) vector.
  q_delta_weights_full <- numeric(n_combinations)
  q_delta_weights_full[active_idx] <- q_active

  list(
    lambda = as.vector(lambda),
    q_delta_weights = q_delta_weights_full,
    iterations = 0L,
    converged = TRUE,
    constraint_norm = 0,
    adaptive_diagnostics = NULL
  )
}

#' Handle zero and near-zero propensity scores
#'
#' @param propensity_scores_i Vector of propensity scores
#' @param divergence Type of divergence
#' @param m Number of treatments (constraints)
#' @param combination_matrix Optional `2^m x m` matrix of all binary treatment
#'   combinations, as returned by [generate_binary_vectors()]. `NULL` (default)
#'   generates it internally; `newton_qdelta()` threads down the matrix it already
#'   built so this function does not regenerate it per observation.
#' @return Vector of propensity scores
#' @keywords internal
handle_zero_propensities <- function(propensity_scores_i, divergence, m,
                                    combination_matrix = NULL) {

  n_combinations <- length(propensity_scores_i)
  min_threshold <- PROPENSITY_ACTIVE_THRESHOLD

  # Count non-zero combinations
  active_indices <- which(propensity_scores_i > min_threshold)
  n_nonzeros <- length(active_indices)

  # If all propensities are zero/near-zero, this is impossible
  if (n_nonzeros == 0) {
    stop("All propensity scores are zero or near-zero (< ", min_threshold, ") - cannot compute Q_delta.\n",
         "This indicates a severe data issue:\n",
         "  1. Check your propensity score estimation - it may have failed\n",
         "  2. Verify your data contains valid treatment assignments\n",
         "  3. Check for numerical underflow in probability calculations\n",
         "Maximum propensity score: ", sprintf("%.2e", max(propensity_scores_i)))
  }

  # Check for under-identification: rank of active vectors matrix
  # We need the support to span the full space of m+1 constraints
  if (is.null(combination_matrix)) {
    combination_matrix <- generate_binary_vectors(m)
  }
  active_vectors <- combination_matrix[active_indices, , drop = FALSE]
  design_matrix <- cbind(1, active_vectors)

  # Use qr rank check
  rank <- qr(design_matrix)$rank

  if (rank < m + 1) {
    stop(sprintf(
      "Under-identification: Support rank (%d) < number of constraints (%d).\n",
      rank, m + 1),
      "The target is under-identified for this observation.\n",
      "Possible causes:\n",
      "  1. Some treatments are redundant (perfectly correlated)\n",
      "  2. Some treatments are constant (no variation)\n",
      "  3. Insufficient treatment variation in this observation's neighborhood\n",
      "Solutions:\n",
      "  1. Reduce the number of treatments (m)\n",
      "  2. Check your data for treatment variation\n",
      "  3. Remove redundant treatments\n",
      "Active combinations: ", n_nonzeros, " out of ", n_combinations
    )
  }

  # For reverse_kl, we still want to ensure a small safe minimum for the active set
  # to avoid numerical issues with the h-function derivatives (1/linear_term^2)
  if (divergence == "reverse_kl") {
    min_safe <- MIN_SAFE_VALUE
    if (any(propensity_scores_i[active_indices] < min_safe)) {
      propensity_scores_i[active_indices][propensity_scores_i[active_indices] < min_safe] <- min_safe
      # Re-normalize ONLY the active set
      sum_active <- sum(propensity_scores_i[active_indices])
      propensity_scores_i[active_indices] <- propensity_scores_i[active_indices] / sum_active
    }
  }

  return(propensity_scores_i)
}

#' Compute target marginals for given delta
#'
#' Computes target marginal probabilities using the odds scaling formula:
#' Q_delta(A_k = 1) = delta_k * P_0(A_k = 1) / (1 + (delta_k - 1) * P_0(A_k = 1))
#'
#' @param propensity_scores_i Numeric vector, propensity scores for observation i
#' @param treatment_combinations Matrix, all treatment combinations. Optional if `m` is provided.
#' @param delta Numeric vector, odds scaling factors
#' @param m Integer, number of treatments. Required if `treatment_combinations` is NULL. Can be inferred from `length(delta)` as convenience fallback.
#' @return Numeric vector, target marginal probabilities
compute_target_marginals <- function(propensity_scores_i, treatment_combinations = NULL, delta, m = NULL) {

  # Handle treatment_combinations parameter
  if (is.null(treatment_combinations)) {
    if (is.null(m)) {
      # Convenience: infer from delta if available
      if (!is.null(delta)) {
        m <- length(delta)
      } else {
        stop("Either 'treatment_combinations' or 'm' must be provided")
      }
    }
    treatment_combinations <- generate_binary_vectors(m)
  } else {
    m <- ncol(treatment_combinations)
  }

  # Current marginals P_0(A_k = 1 | X_i)
  p0_marginals <- numeric(m)
  for (k in 1:m) {
    p0_marginals[k] <- sum(propensity_scores_i * treatment_combinations[, k])
  }

  # Target marginals using odds scaling formula (from theory, main.tex line 81)
  # Q_delta(A_k = 1) = delta_k * P_0(A_k = 1) / (1 + (delta_k - 1) * P_0(A_k = 1))
  target_marginals <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)

  return(target_marginals)
}


#' Get initial lambda for m-dimensional system
#'
#' @param delta Numeric vector, odds scaling factors (length m)
#' @param p0_active Numeric vector, active propensity scores (unused; kept for
#'   signature stability)
#' @param vectors_active Matrix, active treatment combinations
#' @param divergence Character, divergence type
#' @return Numeric vector of length m
#' @keywords internal
get_initial_lambda_active <- function(delta, p0_active, vectors_active, divergence) {
  m <- length(delta)
  if (divergence == "kl") {
    lambda_active <- delta * 0.1
  } else if (divergence == "reverse_kl") {
    lambda_active <- -abs(delta) * 0.01
    test_linear_terms <- as.vector(vectors_active %*% lambda_active)
    if (any(test_linear_terms >= 0)) {
      lambda_active <- -abs(delta) * 0.5
    }
  } else if (divergence == "chi2") {
    lambda_active <- rep(0, m)
  } else if (divergence == "hellinger") {
    lambda_active <- delta * 0.05
  } else {
    lambda_active <- rep(0, m)
  }
  return(lambda_active)
}
