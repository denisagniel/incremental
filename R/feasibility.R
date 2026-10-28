#' @file feasibility.R
#' @title Delta Feasibility Optimization
#' @description
#' Functions for finding and optimizing feasible delta values. These functions
#' help identify maximum feasible deltas, project infeasible deltas onto the
#' feasible set, and characterize the feasible region.
#' @author Denis Agniel

#' Find Maximum Feasible Uniform Delta
#'
#' @description
#' Finds the largest scalar d such that delta = (d, d, ..., d) is feasible
#' for all observations. Uses binary search for efficiency.
#'
#' @param propensity_scores Matrix (n × 2^m) of propensity scores
#' @param divergence Character, divergence type ("kl", "reverse_kl", "chi2", "hellinger")
#' @param tol Numeric, tolerance for binary search convergence (default 1e-3)
#' @param max_iter Integer, maximum binary search iterations (default 50)
#' @param upper_bound Numeric, initial upper bound for search (default 10.0)
#' @param verbose Logical, print progress messages (default FALSE)
#'
#' @details
#' Performs binary search over the scalar d to find the maximum value such that
#' delta = rep(d, m) is feasible for all n observations. This is the most common
#' use case - finding the maximum uniform increase across all treatments.
#'
#' The search starts with lower bound d = 1.0 (no change) and upper bound specified
#' by \code{upper_bound}. If the upper bound is feasible, it doubles until an
#' infeasible value is found.
#'
#' @return List containing:
#' \describe{
#'   \item{max_delta}{Numeric scalar, maximum feasible uniform delta value}
#'   \item{delta_vector}{Numeric vector (length m), the maximum feasible delta as a vector}
#'   \item{n_feasible}{Integer, number of observations feasible at maximum (should be n)}
#'   \item{iterations}{Integer, number of binary search iterations}
#'   \item{converged}{Logical, whether binary search converged within tolerance}
#' }
#'
#' @examples
#' \dontrun{
#' # Estimate propensity scores
#' ps_est <- estimate_propensity_scores(data, method = "glm")
#'
#' # Find maximum uniform delta
#' result <- find_max_uniform_delta(ps_est$propensity_scores)
#' cat("Maximum uniform delta:", result$max_delta, "\n")
#'
#' # Use with safety margin
#' safe_delta <- rep(result$max_delta * 0.9, m)
#' }
#'
#' @export
find_max_uniform_delta <- function(propensity_scores,
                                   divergence = "kl",
                                   tol = 1e-3,
                                   max_iter = 50,
                                   upper_bound = 10.0,
                                   verbose = FALSE) {

  # Validate inputs
  if (!is.matrix(propensity_scores)) {
    propensity_scores <- as.matrix(propensity_scores)
  }

  n <- nrow(propensity_scores)
  n_combinations <- ncol(propensity_scores)
  m <- log2(n_combinations)

  if (m != round(m)) {
    stop("Number of columns in propensity_scores must be a power of 2 (2^m)")
  }
  m <- as.integer(m)

  if (tol <= 0) stop("tol must be positive")
  if (max_iter <= 0) stop("max_iter must be positive")

  # Initialize bounds
  left <- 1.0  # delta = 1 means no change, always feasible
  right <- upper_bound

  # Check if upper bound is feasible - if so, double until infeasible
  iter <- 0
  while (iter < 10) {  # Max 10 doublings
    delta_test <- rep(right, m)
    feasibility <- analyze_delta_feasibility_all_obs(
      propensity_scores, delta_test, tol = 1e-8,
      early_stop = TRUE  # Stop on first infeasible (faster)
    )

    if (feasibility$n_feasible < n) {
      # Found infeasible upper bound
      break
    }

    # Upper bound still feasible, double it
    left <- right
    right <- right * 2

    if (verbose) {
      cat("Upper bound", left, "feasible, trying", right, "\n")
    }

    iter <- iter + 1
  }

  if (iter == 10) {
    warning("Upper bound search did not find infeasible value after 10 doublings. ",
            "Maximum delta may be > ", right)
  }

  # Binary search
  iter <- 0
  converged <- FALSE

  if (verbose) {
    cat("Binary search: initial bounds [", left, ",", right, "]\n")
  }

  while ((right - left) > tol && iter < max_iter) {
    mid <- (left + right) / 2
    delta_test <- rep(mid, m)

    # Test feasibility (stop on first infeasible for speed)
    feasibility <- analyze_delta_feasibility_all_obs(
      propensity_scores, delta_test, tol = 1e-8,
      early_stop = TRUE  # Stop on first infeasible (faster)
    )

    if (feasibility$n_feasible == n) {
      # All checked observations were feasible, try larger delta
      left <- mid
    } else {
      # Found at least one infeasible observation, try smaller delta
      right <- mid
    }

    iter <- iter + 1

    if (verbose && iter %% 5 == 0) {
      cat("Iteration", iter, ": bounds [", round(left, 4), ",", round(right, 4), "]\n")
    }
  }

  converged <- (right - left) <= tol

  if (verbose) {
    cat("Binary search completed in", iter, "iterations\n")
    cat("Maximum uniform delta:", left, "\n")
  }

  # Final verification
  final_delta <- rep(left, m)
  final_feasibility <- analyze_delta_feasibility_all_obs(
    propensity_scores, final_delta, tol = 1e-8,
    early_stop = FALSE  # Need complete verification
  )

  return(list(
    max_delta = left,
    delta_vector = final_delta,
    n_feasible = final_feasibility$n_feasible,
    iterations = iter,
    converged = converged
  ))
}


#' Find Maximum Feasible Delta Per Treatment
#'
#' @description
#' Finds the maximum feasible delta for each treatment individually, holding
#' all other treatments at a baseline value (typically 1.0).
#'
#' @param propensity_scores Matrix (n × 2^m) of propensity scores
#' @param baseline_delta Numeric vector (length m), baseline delta values for other treatments (default rep(1, m))
#' @param divergence Character, divergence type ("kl", "reverse_kl", "chi2", "hellinger")
#' @param tol Numeric, tolerance for binary search convergence (default 1e-3)
#' @param max_iter Integer, maximum binary search iterations per treatment (default 50)
#' @param upper_bound Numeric, initial upper bound for search (default 10.0)
#' @param verbose Logical, print progress messages (default FALSE)
#'
#' @details
#' For each treatment k, performs binary search to find the maximum delta[k]
#' such that the resulting delta vector is feasible, with all other treatments
#' set to their baseline values.
#'
#' This is useful for understanding per-treatment constraints and identifying
#' which treatments are most constrained by the support.
#'
#' @return List containing:
#' \describe{
#'   \item{max_deltas}{Numeric vector (length m), maximum feasible delta for each treatment}
#'   \item{baseline_delta}{Numeric vector (length m), the baseline used}
#'   \item{iterations}{Integer vector (length m), iterations per treatment}
#'   \item{converged}{Logical vector (length m), convergence status per treatment}
#' }
#'
#' @examples
#' \dontrun{
#' # Find maximum per-treatment deltas
#' result <- find_max_per_treatment_delta(ps_est$propensity_scores)
#'
#' # Show results
#' for (k in 1:length(result$max_deltas)) {
#'   cat("Treatment", k, ": max delta =", round(result$max_deltas[k], 3), "\n")
#' }
#'
#' # Identify most constrained treatment
#' most_constrained <- which.min(result$max_deltas)
#' cat("Most constrained: treatment", most_constrained, "\n")
#' }
#'
#' @export
find_max_per_treatment_delta <- function(propensity_scores,
                                         baseline_delta = NULL,
                                         divergence = "kl",
                                         tol = 1e-3,
                                         max_iter = 50,
                                         upper_bound = 10.0,
                                         verbose = FALSE) {

  # Validate inputs
  if (!is.matrix(propensity_scores)) {
    propensity_scores <- as.matrix(propensity_scores)
  }

  n <- nrow(propensity_scores)
  n_combinations <- ncol(propensity_scores)
  m <- log2(n_combinations)

  if (m != round(m)) {
    stop("Number of columns in propensity_scores must be a power of 2 (2^m)")
  }
  m <- as.integer(m)

  # Set baseline
  if (is.null(baseline_delta)) {
    baseline_delta <- rep(1.0, m)
  }

  if (length(baseline_delta) != m) {
    stop("baseline_delta must have length m = ", m)
  }

  # Storage for results
  max_deltas <- numeric(m)
  iterations_per_treatment <- integer(m)
  converged_per_treatment <- logical(m)

  # Search for each treatment
  for (k in 1:m) {
    if (verbose) {
      cat("\n--- Treatment", k, "---\n")
    }

    # Binary search for this treatment
    left <- baseline_delta[k]  # Start from baseline
    right <- upper_bound

    # Find infeasible upper bound
    iter_expand <- 0
    while (iter_expand < 10) {
      delta_test <- baseline_delta
      delta_test[k] <- right

      feasibility <- analyze_delta_feasibility_all_obs(
        propensity_scores, delta_test, tol = 1e-8
      )

      if (feasibility$n_feasible < n) {
        break
      }

      left <- right
      right <- right * 2
      iter_expand <- iter_expand + 1
    }

    # Binary search
    iter <- 0
    converged <- FALSE

    while ((right - left) > tol && iter < max_iter) {
      mid <- (left + right) / 2

      delta_test <- baseline_delta
      delta_test[k] <- mid

      feasibility <- analyze_delta_feasibility_all_obs(
        propensity_scores, delta_test, tol = 1e-8
      )

      if (feasibility$n_feasible == n) {
        left <- mid
      } else {
        right <- mid
      }

      iter <- iter + 1
    }

    converged <- (right - left) <= tol

    max_deltas[k] <- left
    iterations_per_treatment[k] <- iter
    converged_per_treatment[k] <- converged

    if (verbose) {
      cat("Treatment", k, ": max =", round(left, 4),
          "(", iter, "iterations,",
          if(converged) "converged" else "NOT converged", ")\n")
    }
  }

  return(list(
    max_deltas = max_deltas,
    baseline_delta = baseline_delta,
    iterations = iterations_per_treatment,
    converged = converged_per_treatment
  ))
}


#' Find Maximum Feasible Delta in a Given Direction
#'
#' @description
#' Finds the largest scalar α such that delta = α * direction is feasible
#' for all observations. This generalizes uniform delta (direction = rep(1, m))
#' and allows finding maximum deltas along any direction in R^m.
#'
#' @param propensity_scores Matrix (n × 2^m) of propensity scores
#' @param direction Numeric vector (length m), direction to search along
#' @param divergence Character, divergence type ("kl", "reverse_kl", "chi2", "hellinger")
#' @param tol Numeric, tolerance for binary search convergence (default 1e-3)
#' @param max_iter Integer, maximum binary search iterations (default 50)
#' @param upper_bound Numeric, initial upper bound for scalar search (default 10.0)
#' @param normalize_direction Logical, normalize direction to unit length (default FALSE)
#' @param verbose Logical, print progress messages (default FALSE)
#'
#' @details
#' This function performs a 1D search along a ray: delta = α * direction,
#' starting from α = 0 (if all direction components ≥ 0) or α where
#' delta[k] = 1 for the most restrictive k.
#'
#' The direction vector can be any values, but must satisfy:
#' - If direction[k] < 0: delta[k] = 1 + α * direction[k] ≥ 1 for feasible α
#' - All components of (α * direction) must keep delta[k] ≥ 1
#'
#' Common use cases:
#' - Uniform: direction = rep(1, m)
#' - Single treatment: direction = c(0, ..., 1, ..., 0) (1 in position k)
#' - Trade-offs: direction = c(1, -0.5, 0, ...) (increase A1, decrease A2)
#'
#' @return List containing:
#' \describe{
#'   \item{max_alpha}{Numeric scalar, maximum feasible scaling factor}
#'   \item{max_delta}{Numeric vector (length m), maximum feasible delta = alpha * direction}
#'   \item{direction}{Numeric vector (length m), the direction used (possibly normalized)}
#'   \item{n_feasible}{Integer, number of observations feasible at maximum}
#'   \item{iterations}{Integer, number of binary search iterations}
#'   \item{converged}{Logical, whether binary search converged}
#' }
#'
#' @examples
#' \dontrun{
#' # Uniform direction (same as find_max_uniform_delta)
#' result <- find_max_delta_in_direction(
#'   ps_est$propensity_scores,
#'   direction = rep(1, m)
#' )
#'
#' # Single treatment k=2
#' direction <- rep(0, m)
#' direction[2] <- 1
#' result <- find_max_delta_in_direction(
#'   ps_est$propensity_scores,
#'   direction = direction
#' )
#'
#' # Trade-off: increase A1, decrease A2
#' result <- find_max_delta_in_direction(
#'   ps_est$propensity_scores,
#'   direction = c(1, -0.5, rep(0, m-2))
#' )
#' }
#'
#' @export
find_max_delta_in_direction <- function(propensity_scores,
                                        direction,
                                        divergence = "kl",
                                        tol = 1e-3,
                                        max_iter = 50,
                                        upper_bound = 10.0,
                                        normalize_direction = FALSE,
                                        verbose = FALSE) {

  # Validate inputs
  if (!is.matrix(propensity_scores)) {
    propensity_scores <- as.matrix(propensity_scores)
  }

  n <- nrow(propensity_scores)
  n_combinations <- ncol(propensity_scores)
  m <- log2(n_combinations)

  if (m != round(m)) {
    stop("Number of columns in propensity_scores must be a power of 2 (2^m)")
  }
  m <- as.integer(m)

  if (length(direction) != m) {
    stop("direction must have length m = ", m)
  }

  if (all(direction == 0)) {
    stop("direction cannot be zero vector")
  }

  # Normalize if requested
  if (normalize_direction) {
    direction <- direction / sqrt(sum(direction^2))
  }

  # Find lower bound: minimum alpha such that all delta[k] >= 1
  # delta[k] = 1 + alpha * direction[k] >= 1
  # If direction[k] > 0: no constraint from below
  # If direction[k] < 0: alpha >= (1 - 1) / direction[k] = 0
  # If direction[k] = 0: no constraint
  #
  # Actually, we want delta[k] to start from 1:
  # delta[k] = 1 + alpha * direction[k]
  # So alpha = 0 gives delta = rep(1, m)

  # For negative directions, find minimum alpha
  negative_dirs <- direction < 0
  if (any(negative_dirs)) {
    # delta[k] >= 1 means 1 + alpha * direction[k] >= 1
    # For direction[k] < 0: alpha <= 0 / direction[k] = 0
    # So we can't go negative in alpha for decreasing directions
    # Set lower bound to 0
    left <- 0.0
  } else {
    # All directions non-negative, start from alpha = 0
    left <- 0.0
  }

  right <- upper_bound

  # Helper function to compute delta from alpha
  compute_delta <- function(alpha) {
    rep(1.0, m) + alpha * direction
  }

  # Find infeasible upper bound
  iter_expand <- 0
  while (iter_expand < 10) {
    delta_test <- compute_delta(right)

    # Check that all delta[k] >= 1 (required constraint)
    if (any(delta_test < 1.0)) {
      # This alpha gives delta < 1 for some treatment, not allowed
      right <- right / 2
      iter_expand <- iter_expand + 1
      next
    }

    feasibility <- analyze_delta_feasibility_all_obs(
      propensity_scores, delta_test, tol = 1e-8,
      early_stop = TRUE  # Stop on first infeasible (faster)
    )

    if (feasibility$n_feasible < n) {
      break
    }

    left <- right
    right <- right * 2
    iter_expand <- iter_expand + 1
  }

  # Binary search
  iter <- 0
  converged <- FALSE

  if (verbose) {
    cat("Binary search along direction:", paste(round(direction, 3), collapse = ", "), "\n")
    cat("Initial bounds: alpha ∈ [", left, ",", right, "]\n")
  }

  while ((right - left) > tol && iter < max_iter) {
    mid <- (left + right) / 2
    delta_test <- compute_delta(mid)

    # Check delta >= 1 constraint
    if (any(delta_test < 1.0 - 1e-10)) {
      # Infeasible due to delta < 1
      right <- mid
      iter <- iter + 1
      next
    }

    feasibility <- analyze_delta_feasibility_all_obs(
      propensity_scores, delta_test, tol = 1e-8,
      early_stop = TRUE  # Stop on first infeasible (faster)
    )

    if (feasibility$n_feasible == n) {
      left <- mid
    } else {
      right <- mid
    }

    iter <- iter + 1

    if (verbose && iter %% 5 == 0) {
      cat("Iteration", iter, ": alpha ∈ [", round(left, 4), ",", round(right, 4), "]\n")
    }
  }

  converged <- (right - left) <= tol

  # Final result
  max_alpha <- left
  max_delta <- compute_delta(max_alpha)

  # Final verification
  final_feasibility <- analyze_delta_feasibility_all_obs(
    propensity_scores, max_delta, tol = 1e-8,
    early_stop = FALSE  # Need complete verification
  )

  if (verbose) {
    cat("Search completed in", iter, "iterations\n")
    cat("Maximum alpha:", max_alpha, "\n")
    cat("Maximum delta:", paste(round(max_delta, 3), collapse = ", "), "\n")
  }

  return(list(
    max_alpha = max_alpha,
    max_delta = max_delta,
    direction = direction,
    n_feasible = final_feasibility$n_feasible,
    iterations = iter,
    converged = converged
  ))
}


#' Compute Feasible Bounds (Axis-Aligned Bounding Box)
#'
#' @description
#' Computes axis-aligned bounding box for the feasible region by finding
#' minimum and maximum feasible delta for each treatment coordinate.
#'
#' @param propensity_scores Matrix (n × 2^m) of propensity scores
#' @param divergence Character, divergence type ("kl", "reverse_kl", "chi2", "hellinger")
#' @param lower_bound Numeric, lower bound for search (default 1.0)
#' @param upper_bound Numeric, upper bound for search (default 10.0)
#' @param baseline_delta Numeric vector (length m), baseline for other coordinates (default rep(1, m))
#' @param tol Numeric, tolerance for binary search (default 1e-3)
#' @param verbose Logical, print progress (default FALSE)
#'
#' @details
#' For each treatment k, finds the range [L_k, U_k] of feasible values while
#' holding other treatments at baseline. The feasible region is contained in
#' the hyperrectangle: ∏_{k=1}^m [L_k, U_k].
#'
#' This provides a conservative outer approximation of the feasible region.
#'
#' @return List containing:
#' \describe{
#'   \item{lower_bounds}{Numeric vector (length m), lower bounds per treatment}
#'   \item{upper_bounds}{Numeric vector (length m), upper bounds per treatment}
#'   \item{uniform_lower}{Numeric scalar, max of lower bounds (uniform constraint)}
#'   \item{uniform_upper}{Numeric scalar, min of upper bounds (uniform constraint)}
#' }
#'
#' @examples
#' \dontrun{
#' bounds <- compute_feasible_bounds(ps_est$propensity_scores)
#'
#' cat("Feasible region (axis-aligned box):\n")
#' for (k in 1:length(bounds$lower_bounds)) {
#'   cat("  delta[", k, "] ∈ [",
#'       round(bounds$lower_bounds[k], 2), ", ",
#'       round(bounds$upper_bounds[k], 2), "]\n", sep = "")
#' }
#' }
#'
#' @export
compute_feasible_bounds <- function(propensity_scores,
                                   divergence = "kl",
                                   lower_bound = 1.0,
                                   upper_bound = 10.0,
                                   baseline_delta = NULL,
                                   tol = 1e-3,
                                   verbose = FALSE) {

  # Get dimensions
  n_combinations <- ncol(propensity_scores)
  m <- as.integer(log2(n_combinations))

  if (is.null(baseline_delta)) {
    baseline_delta <- rep(1.0, m)
  }

  # Use find_max_per_treatment_delta for upper bounds
  max_result <- find_max_per_treatment_delta(
    propensity_scores,
    baseline_delta = baseline_delta,
    divergence = divergence,
    tol = tol,
    upper_bound = upper_bound,
    verbose = verbose
  )

  upper_bounds <- max_result$max_deltas

  # For lower bounds, we typically have delta >= 1 (no decrease)
  # But we can check if smaller values are feasible
  lower_bounds <- rep(lower_bound, m)

  # Compute uniform bounds (intersection)
  uniform_lower <- max(lower_bounds)
  uniform_upper <- min(upper_bounds)

  return(list(
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    uniform_lower = uniform_lower,
    uniform_upper = uniform_upper,
    baseline_delta = baseline_delta
  ))
}

# =============================================================================
# FAST SAMPLE-BASED METHODS FOR LARGE DATASETS
# =============================================================================

#' Find Maximum Uniform Delta (Fast, Sample-Based)
#'
#' Fast approximation for large datasets using random sampling for binary search.
#' Much faster than full dataset search with minimal accuracy loss.
#'
#' @param propensity_scores Matrix of propensity scores (n × 2^m)
#' @param sample_size Number of observations to sample for search (default: 10000)
#' @param verify_full Logical; verify result on full dataset after search (default: TRUE)
#' @param tol Numeric, convergence tolerance for binary search
#' @param max_iter Maximum iterations for binary search  
#' @param upper_bound Starting upper bound for search
#' @param verbose Logical; print progress messages
#'
#' @details
#' Uses random sampling to dramatically speed up binary search:
#' 1. Sample n_sample observations
#' 2. Run binary search on sample
#' 3. Optionally verify on full dataset with 5% safety margin
#'
#' For 381K observations, reduces search from hours to minutes.
#'
#' @return List containing max_delta, sample_based=TRUE, and verification info
#' @export
find_max_uniform_delta_fast <- function(propensity_scores,
                                        sample_size = 10000,
                                        verify_full = TRUE,
                                        tol = 1e-3,
                                        max_iter = 50,
                                        upper_bound = 10.0,
                                        verbose = FALSE) {

  n <- nrow(propensity_scores)
  m <- log2(ncol(propensity_scores))
  m <- as.integer(m)

  # Sample observations for search
  actual_sample_size <- min(sample_size, n)
  sample_idx <- sample(n, actual_sample_size, replace = FALSE)
  ps_sample <- propensity_scores[sample_idx, , drop = FALSE]

  if (verbose) {
    cat("Sample-based search: using", actual_sample_size, "/", n, "observations\n\n")
  }

  # Binary search on sample (reuse existing function logic)
  left <- 1.0
  right <- upper_bound

  # Find infeasible upper bound
  iter <- 0
  while (iter < 10) {
    delta_test <- rep(right, m)
    feasibility <- analyze_delta_feasibility_all_obs(
      ps_sample, delta_test, tol = 1e-8, early_stop = TRUE
    )

    if (feasibility$n_feasible < actual_sample_size) {
      break
    }

    left <- right
    right <- right * 2
    iter <- iter + 1
  }

  # Binary search
  iter <- 0
  while ((right - left) > tol && iter < max_iter) {
    mid <- (left + right) / 2
    delta_test <- rep(mid, m)

    feasibility <- analyze_delta_feasibility_all_obs(
      ps_sample, delta_test, tol = 1e-8, early_stop = TRUE
    )

    if (feasibility$n_feasible == actual_sample_size) {
      left <- mid
    } else {
      right <- mid
    }

    iter <- iter + 1

    if (verbose && iter %% 5 == 0) {
      cat("Iteration", iter, ": bounds [", round(left, 4), ",", round(right, 4), "]\n")
    }
  }

  max_delta_sample <- left

  if (verbose) {
    cat("\nSample-based max delta:", round(max_delta_sample, 4), "\n")
  }

  # Verify on full dataset if requested
  if (verify_full && n > actual_sample_size) {
    if (verbose) {
      cat("Verifying on full dataset...\n")
    }

    # Apply 5% safety margin
    max_delta_safe <- max_delta_sample * 0.95
    delta_verify <- rep(max_delta_safe, m)

    feasibility_full <- analyze_delta_feasibility_all_obs(
      propensity_scores, delta_verify, tol = 1e-8, early_stop = TRUE
    )

    if (verbose) {
      cat("With 5% safety margin (", round(max_delta_safe, 4), "):\n", sep = "")
      cat("  Feasible:", feasibility_full$n_feasible, "/", n, "observations\n")
    }

    return(list(
      max_delta = max_delta_safe,
      delta_vector = rep(max_delta_safe, m),
      sample_based = TRUE,
      sample_size = actual_sample_size,
      n_total = n,
      verified_full = TRUE,
      n_feasible_full = feasibility_full$n_feasible,
      iterations = iter
    ))

  } else {
    # No full verification - return sample estimate with warning
    if (n > actual_sample_size) {
      warning("Sample-based estimate not verified on full dataset. ",
              "Set verify_full=TRUE for verification.")
    }

    return(list(
      max_delta = max_delta_sample,
      delta_vector = rep(max_delta_sample, m),
      sample_based = TRUE,
      sample_size = actual_sample_size,
      n_total = n,
      verified_full = FALSE,
      iterations = iter
    ))
  }
}


#' Quick Grid Search for Feasible Delta
#'
#' Fast feasibility check using grid search on sample.
#' Tests a sequence of delta values and returns the largest feasible one.
#'
#' @param propensity_scores Matrix of propensity scores
#' @param delta_grid Numeric vector of delta values to test (default: seq(1.1, 2.0, 0.1))
#' @param sample_size Number of observations to sample (default: 10000)
#' @param verbose Logical
#'
#' @details
#' Extremely fast method for large datasets:
#' 1. Sample observations
#' 2. Test each delta in grid sequentially
#' 3. Stop at first infeasible (with early stopping)
#' 4. Return largest feasible with 5% safety margin
#'
#' Trade-off: Less precise than binary search, but much faster.
#' Suitable for quick diagnostics.
#'
#' @return List with max feasible delta
#' @export
quick_feasible_delta <- function(propensity_scores,
                                 delta_grid = seq(1.1, 2.0, by = 0.1),
                                 sample_size = 10000,
                                 verbose = FALSE) {

  n <- nrow(propensity_scores)
  m <- log2(ncol(propensity_scores))
  m <- as.integer(m)

  # Sample
  actual_sample_size <- min(sample_size, n)
  sample_idx <- sample(n, actual_sample_size, replace = FALSE)
  ps_sample <- propensity_scores[sample_idx, , drop = FALSE]

  if (verbose) {
    cat("Quick grid search: sample", actual_sample_size, "/", n, "\n")
    cat("Testing deltas:", paste(round(delta_grid, 2), collapse = ", "), "\n\n")
  }

  # Test each delta
  max_feasible <- 1.0

  for (delta_val in delta_grid) {
    delta_test <- rep(delta_val, m)

    feasibility <- analyze_delta_feasibility_all_obs(
      ps_sample, delta_test, tol = 1e-8, early_stop = TRUE
    )

    if (feasibility$n_feasible == actual_sample_size) {
      max_feasible <- delta_val
      if (verbose) {
        cat("δ =", delta_val, "✓\n")
      }
    } else {
      if (verbose) {
        cat("δ =", delta_val, "✗ (stopped)\n")
      }
      break  # Found first infeasible
    }
  }

  # Apply safety margin
  max_feasible_safe <- max_feasible * 0.95

  if (verbose) {
    cat("\nMax feasible (sample):", round(max_feasible, 2), "\n")
    cat("With 5% safety margin:", round(max_feasible_safe, 2), "\n")
  }

  return(list(
    max_delta = max_feasible_safe,
    delta_vector = rep(max_feasible_safe, m),
    method = "grid_search",
    sample_size = actual_sample_size,
    n_total = n
  ))
}
