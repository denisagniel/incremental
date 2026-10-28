#' @file fast_qdelta_solver_legacy.R
#' @title Frozen pre-vectorization Q_delta solver (equivalence reference ONLY)
#' @description
#' Verbatim snapshot of the scalar-loop Q_delta Newton solver as it stood
#' immediately before the 2026-08-11 vectorization
#' (`quality_reports/specs/2026-08-10_qdelta-solver-performance-fix.md`, MUST
#' item "Numerical equivalence, old vs. new": *"Keep the current implementation
#' accessible ... for direct in-session comparison -- do not rely on git
#' history/checkout for this."*).
#'
#' "Current implementation" for this snapshot means the **post-bugfix** code of
#' 2026-08-10, i.e. it already includes both prerequisite fixes:
#'   1. `divergence_h_pair()` unification across all h/h' call sites
#'      (`R/divergence_h_functions.R`);
#'   2. the `prev_constraint_norm` update moved to the END of the Newton loop
#'      body, so the two "stalled" checks compare against the previous
#'      iteration's norm (see `tests/testthat/test-qdelta-newton-stall-bugfix.R`).
#'
#' Every function here is suffixed `_legacy` and is called by nothing in the
#' package except `tests/testthat/test-qdelta-vectorized-equivalence.R`. It is
#' NOT exported, NOT a fallback, and NOT maintained: it exists so the vectorized
#' rewrite can be diffed numerically against the exact code that produced the
#' methods paper's existing application results. Functions that the rewrite did
#' NOT change (`compute_constraint_vector_internal()`,
#' `compute_q_delta_weights_internal()`, `compute_jacobian_numerical_internal()`,
#' `compute_target_marginals()`, `solve_qdelta_single_treatment()`,
#' `compute_constraint_vector_m()`, `compute_q_delta_weights()`,
#' `compute_jacobian_numerical_m()`, `get_initial_lambda_active()`) are shared
#' with `fast_qdelta_solver.R` rather than duplicated, so the comparison isolates
#' exactly the changed code.
#'
#' Do not "fix" anything in this file. If it needs to change, the rewrite is
#' wrong.
#' @author Denis Agniel
#' @keywords internal
#' @name fast_qdelta_solver_legacy
NULL

#' Legacy (scalar-loop) Newton-Raphson solver for per-observation Q_delta
#'
#' @inheritParams newton_qdelta
#' @return List with lambda solution, convergence info, and q_delta weights
#' @keywords internal
newton_qdelta_legacy <- function(propensity_scores_i, delta, divergence,
                         max_iter = 200, tol = 1e-10, use_analytical_jacobian = TRUE,
                         fix_lambda0 = FALSE, adaptive_constraints = FALSE,
                         max_constraint_drops = 2, condition_threshold = 1e10,
                         verbose_adaptive = FALSE) {

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

  # Initialize adaptive constraint tracking
  constraints_dropped_history <- list()
  total_adaptive_drops <- 0
  adaptive_diagnostics <- list()

  # Handle zero/near-zero propensity scores
  # Now throws error if n_nonzeros < m + 1
  propensity_scores_i <- handle_zero_propensities_legacy(propensity_scores_i, divergence, m)

  # Identify active subset (combinations with p0 > threshold)
  active_idx <- which(propensity_scores_i > PROPENSITY_ACTIVE_THRESHOLD)
  p0_active <- propensity_scores_i[active_idx]

  # Generate all treatment combinations and subset to active
  treatment_combinations_all <- generate_binary_vectors(m)
  treatment_combinations <- treatment_combinations_all[active_idx, , drop = FALSE]

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
    lambda_active <- get_initial_lambda_active_legacy(delta, p0_active, treatment_combinations, divergence)

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
        jacobian <- compute_jacobian_analytical_m_legacy(lambda_active, p0_active,
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
          if (!check_constraint_satisfaction_m_legacy(lambda_new, p0_active,
                                             treatment_combinations, divergence, lambda0)) {
            # If constraints violated, use line search
            lambda_active <- line_search_newton_m_legacy(lambda_active, lambda_new, p0_active,
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
    lambda <- get_initial_lambda_internal_legacy(delta, p0_active, treatment_combinations, divergence)

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
          jacobian <- compute_jacobian_analytical_internal_legacy(lambda, p0_active,
                                                         treatment_combinations, divergence)
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
            if (check_constraint_satisfaction_internal_legacy(lambda_candidate, p0_active,
                                                     treatment_combinations, divergence)) {
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
            check_constraint_satisfaction_internal_legacy(lambda_new, p0_active,
                                                  treatment_combinations, divergence)
          }, error = function(e) {
            FALSE  # If check fails, assume constraints violated
          })

          if (isTRUE(constraint_check)) {
            lambda <- lambda_new
          } else {
            # Line search
            lambda <- line_search_newton_internal_legacy(lambda, lambda_new, p0_active,
                                                treatment_combinations, target_marginals, divergence)
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
        # Check if at domain boundary
        boundary_linear <- sapply(1:nrow(treatment_combinations), function(j) {
          lambda[1] + sum(treatment_combinations[j, ] * lambda[2:(m+1)])
        })
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
    test_linear <- sapply(1:nrow(treatment_combinations), function(j) {
      lambda[1] + sum(treatment_combinations[j, ] * lambda[2:(m+1)])
    })
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
        test_linear <- sapply(1:nrow(treatment_combinations), function(j) {
          lambda[1] + sum(treatment_combinations[j, ] * lambda[2:(m+1)])
        })
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

#' Legacy initialization for internal normalization (scalar sapply loops)
#' @keywords internal
get_initial_lambda_internal_legacy <- function(delta, p0_active, vectors_active, divergence) {

  m <- length(delta)

  if (divergence == "kl") {
    # For KL: start with small values, λ₀ can be 0
    lambda <- c(0, delta * 0.1)

  } else if (divergence == "reverse_kl") {
    # More conservative initialization for reverse_kl
    lambda <- c(-1.0, -abs(delta) * 0.2)

    # Verify all linear terms are safely negative
    test_linear_terms <- sapply(1:nrow(vectors_active), function(j) {
      a_j <- vectors_active[j, ]
      lambda[1] + sum(a_j * lambda[2:(m+1)])
    })

    # Adjust to ensure max is at most -0.1
    if (max(test_linear_terms) >= -0.1) {
      lambda[1] <- lambda[1] + min(test_linear_terms) - 0.2
    }

  } else if (divergence == "chi2") {
    # For chi2: use conservative initialization
    lambda_k <- delta * 0.1
    lambda_0 <- -0.5

    lambda <- c(lambda_0, lambda_k)

    test_linear <- sapply(1:nrow(vectors_active), function(j) {
      lambda[1] + sum(vectors_active[j, ] * lambda[2:(m+1)])
    })

    if (min(test_linear) <= -1.5) {
      # Too close to boundary, scale everything down
      lambda <- lambda * 0.3
      lambda[1] <- -0.3
    }

  } else if (divergence == "hellinger") {
    # For hellinger: need all linear terms < 1
    lambda <- c(-0.1, delta * 0.05)

    # Verify constraints
    test_linear_terms <- sapply(1:nrow(vectors_active), function(j) {
      a_j <- vectors_active[j, ]
      lambda[1] + sum(a_j * lambda[2:(m+1)])
    })

    if (any(test_linear_terms >= 1 - DOMAIN_VALIDATION_STRICT_TOL)) {
      lambda[1] <- -0.5  # More negative
    }
  }

  return(lambda)
}

#' Legacy analytical Jacobian for internal normalization, (m+1)x(m+1) via (i,j) loops
#' @keywords internal
compute_jacobian_analytical_internal_legacy <- function(lambda, p0_active, vectors_active, divergence) {

  m <- length(lambda) - 1

  # Initialize Jacobian matrix (m+1) × (m+1)
  jacobian <- matrix(0, nrow = m + 1, ncol = m + 1)

  # Compute q_delta weights
  q_delta_weights <- compute_q_delta_weights_internal(lambda, p0_active,
                                                    vectors_active, divergence)

  # Compute Jacobian entries
  for (i in 1:(m + 1)) {
    for (j in 1:(m + 1)) {

      if (i <= m) {
        # Marginal constraint derivatives: ∂E_Q[A_i]/∂λ_j
        jacobian[i, j] <- compute_marginal_jacobian_internal_legacy(i, j, lambda, p0_active,
                                                          vectors_active, divergence)
      } else {
        # Normalization constraint derivatives: ∂(Σ q_delta)/∂λ_j
        jacobian[i, j] <- compute_normalization_jacobian_internal_legacy(j, lambda, p0_active,
                                                               vectors_active, divergence)
      }
    }
  }

  return(jacobian)
}

#' Legacy marginal-constraint Jacobian entry (single (i,j) cell, scalar loop over k)
#' @keywords internal
compute_marginal_jacobian_internal_legacy <- function(i, j, lambda, p0_active, vectors_active, divergence) {

  m <- length(lambda) - 1
  n_active <- length(p0_active)

  # ∂E_Q[A_i]/∂λ_j = Σ_k p_0(a_k|X) * a_k[i] * ∂h_k/∂λ_j
  sum_term <- 0

  for (k in 1:n_active) {
    a_k <- vectors_active[k, ]
    linear_term <- lambda[1] + sum(a_k * lambda[2:(m+1)])

    h_prime <- divergence_h_pair(linear_term, divergence)$h_prime
    dh_dlambda <- h_prime * (if (j == 1) 1 else a_k[j - 1])

    sum_term <- sum_term + p0_active[k] * a_k[i] * dh_dlambda
  }

  return(sum_term)
}

#' Legacy normalization-constraint Jacobian entry (single j cell, scalar loop over k)
#' @keywords internal
compute_normalization_jacobian_internal_legacy <- function(j, lambda, p0_active, vectors_active, divergence) {

  m <- length(lambda) - 1
  n_active <- length(p0_active)

  # ∂(Σ q_delta)/∂λ_j = Σ_k p_0(a_k|X) * ∂h_k/∂λ_j
  sum_term <- 0

  for (k in 1:n_active) {
    a_k <- vectors_active[k, ]
    linear_term <- lambda[1] + sum(a_k * lambda[2:(m+1)])

    h_prime <- divergence_h_pair(linear_term, divergence)$h_prime
    dh_dlambda <- h_prime * (if (j == 1) 1 else a_k[j - 1])

    sum_term <- sum_term + p0_active[k] * dh_dlambda
  }

  return(sum_term)
}

#' Legacy domain check for internal normalization (scalar loop, early return)
#' @keywords internal
check_constraint_satisfaction_internal_legacy <- function(lambda, p0_active, vectors_active, divergence) {

  if (!(divergence %in% c("reverse_kl", "chi2", "hellinger"))) {
    return(TRUE)
  }

  m <- length(lambda) - 1
  n_active <- length(p0_active)

  for (j in 1:n_active) {
    a_j <- vectors_active[j, ]
    linear_term <- lambda[1] + sum(a_j * lambda[2:(m+1)])

    if (divergence == "reverse_kl") {
      if (linear_term >= -DOMAIN_VALIDATION_STRICT_TOL) return(FALSE)
    } else if (divergence == "chi2") {
      # For chi2, domain is linear_term > -2 (so h(x) = x/2 + 1 > 0)
      if (linear_term <= CHI2_BOUNDARY_LAMBDA0 + CHI2_DOMAIN_BUFFER) return(FALSE)
    } else if (divergence == "hellinger") {
      if (linear_term >= 1 - DOMAIN_VALIDATION_STRICT_TOL) return(FALSE)
    }
  }

  return(TRUE)
}

#' Legacy Armijo line search for internal normalization (scalar chi2 alpha_max loop)
#' @keywords internal
line_search_newton_internal_legacy <- function(lambda_old, lambda_new, p0_active,
                                      vectors_active, target_marginals, divergence,
                                      max_steps = LINE_SEARCH_MAX_STEPS,
                                      c1 = LINE_SEARCH_ARMIJO_C1) {

  # Current objective (constraint norm squared)
  constraint_old <- compute_constraint_vector_internal(
    lambda_old, p0_active, vectors_active, target_marginals, divergence
  )
  objective_old <- sum(constraint_old^2)

  direction <- lambda_new - lambda_old

  alpha_max <- 1.0
  if (divergence == "chi2") {
    m <- length(lambda_old) - 1
    for (j in 1:nrow(vectors_active)) {
      a_j <- vectors_active[j, ]
      lt_old <- lambda_old[1] + sum(a_j * lambda_old[2:(m+1)])
      lt_dir <- direction[1] + sum(a_j * direction[2:(m+1)])
      if (lt_dir < 0) {
        # Need lt_old + alpha * lt_dir > -2
        safe_alpha <- (lt_old + 2) / (-lt_dir) * 0.999  # 0.1% buffer from boundary
        alpha_max <- min(alpha_max, safe_alpha)
      }
    }
  }

  # Try full step first (alpha_max), then backtrack
  alphas <- c(alpha_max, alpha_max / 2^(1:max_steps))

  for (alpha in alphas) {
    lambda_candidate <- lambda_old + alpha * direction

    # Check domain feasibility
    if (!check_constraint_satisfaction_internal_legacy(
      lambda_candidate, p0_active, vectors_active, divergence
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

#' Legacy analytical Jacobian for the m-dimensional (external normalization) system
#' @keywords internal
compute_jacobian_analytical_m_legacy <- function(lambda_active, p0_active, vectors_active, divergence, lambda0) {
  m <- length(lambda_active)
  n_active <- length(p0_active)
  lambda <- c(lambda0, lambda_active)

  jacobian <- matrix(0, nrow = m, ncol = m)

  # Compute q_delta weights and normalization
  q_delta_weights <- compute_q_delta_weights(lambda, p0_active, vectors_active, divergence)
  Z <- sum(q_delta_weights)
  q_delta_normalized <- q_delta_weights / Z

  for (i in 1:m) {
    for (j in 1:m) {
      sum1 <- 0
      sum2 <- 0
      for (k in 1:n_active) {
        a_k <- vectors_active[k, ]
        linear_term <- lambda[1] + sum(a_k * lambda[2:(m+1)])

        h_prime <- divergence_h_pair(linear_term, divergence)$h_prime
        dh_dlambda <- h_prime * a_k[j]

        sum1 <- sum1 + p0_active[k] * a_k[i] * dh_dlambda
        sum2 <- sum2 + p0_active[k] * dh_dlambda
      }

      E_Ai <- sum(q_delta_normalized * vectors_active[, i])
      jacobian[i, j] <- sum1 / Z - E_Ai * sum2 / Z
    }
  }
  return(jacobian)
}

#' Legacy domain check for the m-dimensional system (scalar loop)
#' @keywords internal
check_constraint_satisfaction_m_legacy <- function(lambda_active, p0_active, vectors_active, divergence, lambda0 = 0) {
  if (!(divergence %in% c("reverse_kl", "chi2", "hellinger"))) {
    return(TRUE)
  }
  m <- length(lambda_active)
  n_active <- length(p0_active)
  for (j in 1:n_active) {
    a_j <- vectors_active[j, ]
    linear_term <- lambda0 + sum(a_j * lambda_active)
    if (divergence == "reverse_kl") {
      if (linear_term >= -NEAR_ZERO_THRESHOLD) return(FALSE)
    } else if (divergence == "chi2") {
      if (linear_term <= CHI2_BOUNDARY_LAMBDA0) return(FALSE)
    } else if (divergence == "hellinger") {
      if (linear_term >= 1) return(FALSE)
    }
  }
  return(TRUE)
}

#' Legacy line search for the m-dimensional system
#' @keywords internal
line_search_newton_m_legacy <- function(lambda_old, lambda_new, p0_active,
                                vectors_active, target_marginals, divergence, lambda0,
                                max_steps = 10) {

  for (step in 1:max_steps) {
    alpha <- 1 / (2^step)
    lambda_candidate <- (1 - alpha) * lambda_old + alpha * lambda_new

    if (check_constraint_satisfaction_m_legacy(lambda_candidate, p0_active,
                                      vectors_active, divergence, lambda0)) {
      return(lambda_candidate)
    }
  }

  # If line search fails, return old lambda
  return(lambda_old)
}

#' Legacy initialization for the m-dimensional system (scalar sapply loop)
#' @keywords internal
get_initial_lambda_active_legacy <- function(delta, p0_active, vectors_active, divergence) {
  m <- length(delta)
  if (divergence == "kl") {
    lambda_active <- delta * 0.1
  } else if (divergence == "reverse_kl") {
    lambda_active <- -abs(delta) * 0.01
    test_linear_terms <- sapply(1:nrow(vectors_active), function(j) {
      a_j <- vectors_active[j, ]
      0 + sum(a_j * lambda_active)
    })
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

#' Legacy zero/near-zero propensity handling (regenerates the combination matrix internally)
#' @keywords internal
handle_zero_propensities_legacy <- function(propensity_scores_i, divergence, m) {

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
    treatment_combinations_all <- generate_binary_vectors(m)
  active_vectors <- treatment_combinations_all[active_indices, , drop = FALSE]
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

  # For reverse_kl, ensure a small safe minimum for the active set
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
