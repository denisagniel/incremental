#' Generate all binary vectors of length m
#'
#' Creates a matrix containing all 2^m possible binary combinations
#' for m binary treatments.
#'
#' @param m Integer, number of binary treatments
#' @return Matrix with 2^m rows and m columns, each row is a binary vector
#' @export
#' @examples
#' generate_binary_vectors(2)
#' generate_binary_vectors(3)
generate_binary_vectors <- local({
  # Memoised by m: the matrix is a plain value (copy-on-modify), so sharing is safe, and several
  # per-observation loops previously rebuilt it for every observation.
  cache <- new.env(parent = emptyenv())
  max_m <- 20L  # 2^20 x 20 doubles is 168 MB (2^24 would be 3 GB and stay cached for the whole session)

  function(m) {
    if (!is.numeric(m) || length(m) != 1L || is.na(m) || m <= 0 || m != as.integer(m)) {
      stop("m must be a positive integer")
    }
    if (m > max_m) {
      stop("m > ", max_m, " is not supported by the dense 2^m enumeration")
    }

    key <- as.character(as.integer(m))
    if (is.null(cache[[key]])) {
      # Doubling construction; column j is bit (j - 1), so column 1 varies fastest, the same
      # order as as.integer(intToBits(i)) for row i + 1.
      vectors <- matrix(0, nrow = 1L, ncol = 0L)
      for (j in seq_len(m)) {
        vectors <- rbind(cbind(vectors, 0), cbind(vectors, 1))
      }
      cache[[key]] <- vectors
    }
    cache[[key]]
  }
})


#' Compute marginal probabilities from joint distribution
#'
#' Extracts marginal probabilities P(A_k = 1) for each treatment k
#' from a joint probability distribution over all treatment combinations.
#'
#' @param joint_dist Numeric vector of length 2^m, joint probabilities
#' @param vectors Matrix of binary vectors (from generate_binary_vectors). Optional if `m` is provided.
#' @param m Integer, number of treatments. Required if `vectors` is NULL.
#' @return Numeric vector of length m, marginal probabilities
#' @export
#' @examples
#' # Using vectors parameter
#' vectors <- generate_binary_vectors(2)
#' joint_dist <- c(0.3, 0.2, 0.1, 0.4)
#' compute_marginals(joint_dist, vectors)
#' 
#' # Using m parameter (vectors will be generated internally)
#' compute_marginals(joint_dist, m = 2)
compute_marginals <- function(joint_dist, vectors = NULL, m = NULL) {
  # Handle vectors parameter
  if (is.null(vectors)) {
    if (is.null(m)) {
      stop("Either 'vectors' or 'm' must be provided")
    }
    vectors <- generate_binary_vectors(m)
  } else {
    m <- ncol(vectors)
  }
  
  # Validate dimensions
  if (length(joint_dist) != nrow(vectors)) {
    stop("joint_dist length must equal number of rows in vectors (expected ", nrow(vectors), ", got ", length(joint_dist), ")")
  }
  
  # Vectorized: compute all marginals at once using matrix multiplication
  # marginals[k] = sum(joint_dist * vectors[, k])
  marginals <- as.numeric(joint_dist %*% vectors)
  
  return(marginals)
}

#' Check if delta values are feasible given support constraints
#'
#' Enhanced version that uses marginal polytope computation to rigorously
#' validate that target marginal probabilities implied by delta lie within
#' the marginal polytope of the support of P_0.
#'
#' @param delta Numeric vector of length m, odds scaling factors
#' @param p0_marginals Numeric vector of length m, observed marginal probabilities
#' @param support_matrix Logical vector of length 2^m, indicating which treatment
#'   combinations have positive probability under P_0
#' @param tol Numeric, tolerance for feasibility check
#' @return List containing:
#'   - feasible: logical, whether delta is feasible
#'   - target_marginals: computed target marginal probabilities
#'   - polytope_vertices: vertices of the marginal polytope
#'   - distance_to_boundary: distance to polytope boundary
#'   - reason: character, reason for infeasibility if applicable
#' @export
#' @examples
#' # Example with m=2, full support
#' delta <- c(1.5, 2.0)
#' p0_marginals <- c(0.3, 0.4)
#' support_matrix <- rep(TRUE, 4)
#' result <- check_delta_feasibility(delta, p0_marginals, support_matrix)
#' 
#' # Example with missing support
#' support_partial <- c(TRUE, TRUE, TRUE, FALSE)
#' result2 <- check_delta_feasibility(delta, p0_marginals, support_partial)
check_delta_feasibility <- function(delta, p0_marginals, support_matrix, tol = 1e-8) {
  
  if (length(delta) != length(p0_marginals)) {
    stop("delta and p0_marginals must have same length")
  }
  
  m <- length(delta)
  
  # Compute target marginals using marginal odds ratio formula
  target_marginals <- delta * p0_marginals / (1 + (delta - 1) * p0_marginals)
  
  # Basic validity check: target marginals must be in [0,1]
  if (any(target_marginals < 0) || any(target_marginals > 1)) {
    return(list(
      feasible = FALSE,
      target_marginals = target_marginals,
      polytope_vertices = NULL,
      distance_to_boundary = 0,
      reason = "Target marginals outside [0,1] range"
    ))
  }
  
  # Check if any target marginal is exactly 0 or 1 (boundary cases)
  if (any(target_marginals == 0) || any(target_marginals == 1)) {
    # For boundary cases, check if they're compatible with support
    for (k in 1:m) {
      if (target_marginals[k] == 0) {
        # If target marginal is 0, no supported combination can have A_k = 1
        vectors <- generate_binary_vectors(m)
        for (i in 1:nrow(vectors)) {
          if (support_matrix[i] && vectors[i, k] == 1) {
            return(list(
              feasible = FALSE,
              target_marginals = target_marginals,
              polytope_vertices = NULL,
              distance_to_boundary = 0,
              reason = paste("Target marginal for A", k, "= 0 but supported combinations have A", k, "= 1")
            ))
          }
        }
      } else if (target_marginals[k] == 1) {
        # If target marginal is 1, all supported combinations must have A_k = 1
        vectors <- generate_binary_vectors(m)
        for (i in 1:nrow(vectors)) {
          if (support_matrix[i] && vectors[i, k] == 0) {
            return(list(
              feasible = FALSE,
              target_marginals = target_marginals,
              polytope_vertices = NULL,
              distance_to_boundary = 0,
              reason = paste("Target marginal for A", k, "= 1 but supported combinations have A", k, "= 0")
            ))
          }
        }
      }
    }
  }
  
  # Full support: the marginal polytope is the unit cube, and every target marginal is strictly
  # inside (0, 1) here (values outside [0, 1] and boundary values were handled above), so the target
  # is feasible and no LP over 2^m vertices is needed. The LP's weights (and the "distance" taken
  # from them) were an arbitrary basic solution; the distance reported here is the margin to the
  # cube's boundary, which is what the field is documented to mean.
  # A target within `tol` of the boundary is NOT taken here: it falls through to the LP, as before.
  if (all(support_matrix) && min(target_marginals, 1 - target_marginals) > tol) {
    return(list(
      feasible = TRUE,
      target_marginals = target_marginals,
      polytope_vertices = compute_marginal_polytope(support_matrix, m),
      distance_to_boundary = min(target_marginals, 1 - target_marginals),
      reason = "Feasible"
    ))
  }

  # Compute marginal polytope vertices
  polytope_result <- tryCatch({
    compute_marginal_polytope(support_matrix, m)
  }, error = function(e) {
    return(list(
      feasible = FALSE,
      target_marginals = target_marginals,
      polytope_vertices = NULL,
      distance_to_boundary = 0,
      reason = paste("Error computing marginal polytope:", e$message)
    ))
  })
  
  # Check if polytope computation failed
  if (is.list(polytope_result) && !is.null(polytope_result$feasible) && !polytope_result$feasible) {
    return(polytope_result)
  }
  
  polytope_vertices <- polytope_result
  
  # Check if target marginals lie within the polytope
  polytope_check <- check_point_in_marginal_polytope(target_marginals, polytope_vertices, tol)
  
  return(list(
    feasible = polytope_check$feasible,
    target_marginals = target_marginals,
    polytope_vertices = polytope_vertices,
    distance_to_boundary = polytope_check$distance_to_boundary,
    reason = if (polytope_check$feasible) "Feasible" else "Target marginals outside marginal polytope"
  ))
}

#' Find active constraint indices
#'
#' Identifies which marginal constraints are active (not at boundary)
#' for a given covariate value.
#'
#' @param p0_marginals Numeric vector, marginal probabilities P_0(A_k = 1 | X = x)
#' @param tol Numeric, tolerance for determining if probability is at boundary
#' @return Integer vector of active constraint indices
#' @export
find_active_constraints <- function(p0_marginals, tol = 1e-8) {
  active <- which(p0_marginals > tol & p0_marginals < (1 - tol))
  return(active)
}

#' Compute marginal polytope vertices from support
#'
#' Enumerates the vertices of the marginal polytope M(supp(P_0)) for a given
#' support pattern. The marginal polytope is the set of all possible marginal
#' probability vectors that can be realized by distributions supported only
#' on the given support.
#'
#' @param support_vector Logical vector of length 2^m, indicating which treatment
#'   combinations have positive probability under P_0
#' @param m Integer, number of treatments
#' @return Matrix with m columns and k rows, where each row is a vertex
#'   of the marginal polytope
#' @export
#' @examples
#' # Full support case (m=2)
#' support_full <- rep(TRUE, 4)
#' vertices_full <- compute_marginal_polytope(support_full, 2)
#' 
#' # Missing one combination (m=2)
#' support_partial <- c(TRUE, TRUE, TRUE, FALSE)
#' vertices_partial <- compute_marginal_polytope(support_partial, 2)
compute_marginal_polytope <- function(support_vector, m) {
  
  if (length(support_vector) != 2^m) {
    stop("support_vector length must equal 2^m")
  }
  
  # Generate all treatment combinations
  vectors <- generate_binary_vectors(m)

  # Find supported combinations
  supported_indices <- which(support_vector)
  if (length(supported_indices) == 0) {
    stop("No supported combinations found")
  }

  # Rows of the 2^m design are distinct, so each supported combination is its own vertex
  # (all mass on that combination puts the marginals at the combination itself): no de-duplication
  # is needed. Full support returns the cached design itself, so per-observation results share it.
  if (length(supported_indices) == nrow(vectors)) {
    return(vectors)
  }

  vectors[supported_indices, , drop = FALSE]
}

#' Check if target marginals lie within marginal polytope
#'
#' Uses linear programming to determine if a target marginal vector lies
#' within the marginal polytope defined by the given vertices.
#'
#' @param target_marginals Numeric vector of length m, target marginal probabilities
#' @param polytope_vertices Matrix with m columns, each row is a vertex
#' @param tol Numeric, tolerance for feasibility check
#' @return List containing:
#'   - feasible: logical, whether point is in polytope
#'   - in_interior: logical, whether point is in relative interior
#'   - weights: numeric vector, convex combination weights (if feasible)
#'   - distance_to_boundary: numeric, distance to polytope boundary
#'
#' @details
#' This function uses the lpSolve package for linear programming.
#' lpSolve is a required dependency for accurate polytope membership checks.
#'
#' @export
#' @examples
#' # Check if point is in polytope
#' vertices <- matrix(c(0,0, 1,0, 0,1, 1,1), ncol=2, byrow=TRUE)
#' target <- c(0.5, 0.5)
#' result <- check_point_in_marginal_polytope(target, vertices)
check_point_in_marginal_polytope <- function(target_marginals, polytope_vertices, tol = 1e-8) {
  
  m <- length(target_marginals)
  n_vertices <- nrow(polytope_vertices)
  
  if (ncol(polytope_vertices) != m) {
    stop("polytope_vertices must have m columns")
  }
  
  # Set up linear program: find weights w_i >= 0, sum(w_i) = 1
  # such that sum(w_i * vertex_i) = target_marginals
  
  # Objective: minimize sum of weights (to find if feasible)
  obj <- rep(1, n_vertices)
  
  # Constraints: 
  # 1. sum(w_i) = 1 (equality constraint)
  # 2. sum(w_i * vertex_i[j]) = target_marginals[j] for each j (m equality constraints)
  # 3. w_i >= 0 for all i (n_vertices inequality constraints)
  
  # Constraint matrix: (1 + m) rows, n_vertices columns
  const_mat <- matrix(0, nrow = 1 + m, ncol = n_vertices)
  
  # Row 1: sum(w_i) = 1
  const_mat[1, ] <- 1
  
  # Rows 2 to (m+1): sum(w_i * vertex_i[j]) = target_marginals[j]
  for (j in 1:m) {
    const_mat[1 + j, ] <- polytope_vertices[, j]
  }
  
  # Right-hand side
  rhs <- c(1, target_marginals)
  
  # Direction (all equality constraints)
  dir <- rep("==", 1 + m)
  
  # Use lpSolve for linear programming (required dependency)
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("lpSolve package is required for polytope checks. Install with: install.packages('lpSolve')")
  }

  # Solve LP
  lp_result <- lpSolve::lp("min", obj, const_mat, dir, rhs)
  
  if (lp_result$status == 0) {
    # Feasible
    weights <- lp_result$solution
    feasible <- TRUE
    
    # Check if in relative interior (all weights > tol)
    in_interior <- all(weights > tol)
    
    # Distance to boundary (minimum positive weight)
    positive_weights <- weights[weights > tol]
    distance_to_boundary <- if (length(positive_weights) > 0) min(positive_weights) else 0
    
  } else {
    # Infeasible
    weights <- rep(0, n_vertices)
    feasible <- FALSE
    in_interior <- FALSE
    distance_to_boundary <- 0
  }
  
  return(list(
    feasible = feasible,
    in_interior = in_interior,
    weights = weights,
    distance_to_boundary = distance_to_boundary
  ))
}

#' Simplified polytope check without LP solver
#'
#' Fallback method when lpSolve is not available.
#'
#' @param target_marginals Numeric vector
#' @param polytope_vertices Matrix
#' @param tol Numeric tolerance
#' @return List with feasibility information
check_point_in_marginal_polytope_simple <- function(target_marginals, polytope_vertices, tol) {
  
  # Simple check: see if target is close to any vertex
  n_vertices <- nrow(polytope_vertices)
  distances <- numeric(n_vertices)
  
  for (i in 1:n_vertices) {
    distances[i] <- sqrt(sum((target_marginals - polytope_vertices[i, ])^2))
  }
  
  min_distance <- min(distances)
  closest_vertex <- which.min(distances)
  
  # If very close to a vertex, consider feasible
  if (min_distance < tol) {
    weights <- rep(0, n_vertices)
    weights[closest_vertex] <- 1
    return(list(
      feasible = TRUE,
      in_interior = FALSE,
      weights = weights,
      distance_to_boundary = 1
    ))
  }
  
  # Otherwise, check if target is in convex hull using simple heuristics
  # This is a very simplified check - not mathematically rigorous
  feasible <- min_distance < 0.1  # Arbitrary threshold
  
  return(list(
    feasible = feasible,
    in_interior = FALSE,
    weights = rep(0, n_vertices),
    distance_to_boundary = min_distance
  ))
}

#' Analyze delta feasibility across all observations
#'
#' Checks feasibility of delta values for each observation and provides
#' aggregate statistics about feasibility across the dataset.
#'
#' @param propensity_scores Matrix of propensity scores (n × 2^m)
#' @param delta Numeric vector of length m, odds scaling factors
#' @param tol Numeric, tolerance for feasibility check
#' @param early_stop Logical; if TRUE (default), stop checking as soon as the first
#'   infeasible observation is found and return immediately with partial results.
#'   If FALSE, checks all n observations and returns full report (slower but complete).
#'   Early stopping is recommended for large datasets (n > 100K) to avoid wasting time
#'   checking all observations when delta is clearly infeasible.
#' @return List containing:
#'   - n_feasible: number of feasible observations
#'   - n_total: total number of observations
#'   - feasibility_rate: proportion of feasible observations
#'   - infeasible_indices: indices of infeasible observations
#'   - polytope_info: list of polytope information for each observation
#'   - most_restrictive_obs: index of observation with tightest polytope
#'   - early_stopped: logical flag indicating if early stopping was triggered (only present if TRUE)
#'   - first_infeasible: index of first infeasible observation (only present if early_stopped=TRUE)
#' @export
#' @examples
#' # Generate example data
#' n <- 100; m <- 2
#' propensity_scores <- matrix(runif(n * 2^m), nrow = n, ncol = 2^m)
#' propensity_scores <- propensity_scores / rowSums(propensity_scores)
#' delta <- c(1.5, 2.0)
#' result <- analyze_delta_feasibility_all_obs(propensity_scores, delta)
analyze_delta_feasibility_all_obs <- function(propensity_scores, delta, tol = 1e-8,
                                              early_stop = TRUE) {
  
  n <- nrow(propensity_scores)
  m <- length(delta)
  
  if (ncol(propensity_scores) != 2^m) {
    stop("propensity_scores must have 2^m columns")
  }
  
  # Generate treatment combinations
  vectors <- generate_binary_vectors(m)
  
  # Check feasibility for each observation
  feasibility_results <- vector("list", n)
  feasible_indices <- logical(n)

  # Progress bar for complete checks on large datasets
  show_progress <- !early_stop && n > 1000
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = n, style = 3, width = 50)
  }

  for (i in 1:n) {
    # Update progress bar
    if (show_progress && i %% 100 == 0) {
      setTxtProgressBar(pb, i)
    }

    # Compute marginals for this observation
    p0_marginals <- compute_marginals(propensity_scores[i, ], vectors)

    # Extract support (combinations with positive probability)
    support_matrix <- propensity_scores[i, ] > tol

    # Check feasibility
    feasibility_result <- check_delta_feasibility(delta, p0_marginals, support_matrix, tol)
    feasibility_results[[i]] <- feasibility_result
    feasible_indices[i] <- feasibility_result$feasible

    # Early stopping if infeasible found
    if (early_stop && !feasibility_result$feasible) {
      # Build partial result and return immediately
      return(list(
        n_feasible = i - 1,  # All before this one were feasible
        n_infeasible = 1,    # At least this one is infeasible
        n_total = n,
        feasibility_rate = (i - 1) / n,
        feasible_per_obs = feasible_indices[1:i],  # Consistent field name
        feasible_indices = which(feasible_indices[1:i]),
        infeasible_indices = i,  # First infeasible observation
        first_infeasible = i,     # Flag for early stop
        feasibility_results = feasibility_results[1:i],
        polytope_info = feasibility_result$polytope_info,  # From infeasible obs
        early_stopped = TRUE   # Flag indicating early stop
      ))
    }
  }

  # Close progress bar if it was shown
  if (show_progress) {
    setTxtProgressBar(pb, n)  # Complete the bar
    close(pb)
    cat("\n")  # New line after progress bar
  }

  # If we get here, all observations were checked (no early stop triggered)
  # Find most restrictive observation (smallest polytope volume)
  polytope_volumes <- numeric(n)
  for (i in 1:n) {
    if (feasibility_results[[i]]$feasible) {
      # Approximate volume by number of vertices (simplified)
      polytope_volumes[i] <- nrow(feasibility_results[[i]]$polytope_vertices)
    } else {
      polytope_volumes[i] <- Inf
    }
  }

  most_restrictive_obs <- which.min(polytope_volumes)

  return(list(
    n_feasible = sum(feasible_indices),
    n_infeasible = sum(!feasible_indices),
    n_total = n,
    feasibility_rate = sum(feasible_indices) / n,
    feasible_per_obs = feasible_indices,  # Add this field for consistency
    infeasible_indices = which(!feasible_indices),
    polytope_info = feasibility_results,
    most_restrictive_obs = most_restrictive_obs
  ))
}

#' Suggest feasible delta values
#'
#' Finds delta values that are feasible for all observations by solving
#' a constrained optimization problem to minimize distance from requested delta
#' while ensuring feasibility across all observations.
#'
#' @param propensity_scores Matrix of propensity scores (n × 2^m)
#' @param requested_delta Numeric vector of length m, requested delta values
#' @param tol Numeric, tolerance for feasibility check
#' @param max_iter Integer, maximum iterations for optimization
#' @return List containing:
#'   - suggested_delta: feasible delta values closest to requested
#'   - feasibility_margin: minimum distance to polytope boundary
#'   - n_obs_feasible: number of observations for which suggested delta is feasible
#'   - optimization_info: details about the optimization process
#' @export
#' @examples
#' # Generate example data
#' n <- 50; m <- 2
#' propensity_scores <- matrix(runif(n * 2^m), nrow = n, ncol = 2^m)
#' propensity_scores <- propensity_scores / rowSums(propensity_scores)
#' requested_delta <- c(3.0, 4.0)  # Potentially infeasible
#' result <- suggest_feasible_delta(propensity_scores, requested_delta)
suggest_feasible_delta <- function(propensity_scores, requested_delta, tol = 1e-8, max_iter = 100) {
  
  n <- nrow(propensity_scores)
  m <- length(requested_delta)
  
  if (ncol(propensity_scores) != 2^m) {
    stop("propensity_scores must have 2^m columns")
  }
  
  # Generate treatment combinations
  vectors <- generate_binary_vectors(m)
  
  # Objective function: minimize ||delta - requested_delta||^2
  objective <- function(delta) {
    sum((delta - requested_delta)^2)
  }
  
  # Constraint function: check feasibility for all observations
  constraint <- function(delta) {
    feasible_count <- 0
    for (i in 1:n) {
      p0_marginals <- compute_marginals(propensity_scores[i, ], vectors)
      support_matrix <- propensity_scores[i, ] > tol
      feasibility_result <- check_delta_feasibility(delta, p0_marginals, support_matrix, tol)
      if (feasibility_result$feasible) {
        feasible_count <- feasible_count + 1
      }
    }
    return(feasible_count - n)  # Should be >= 0 for all feasible
  }
  
  # Use constrained optimization
  # Start from requested delta
  initial_delta <- requested_delta
  
  # Set bounds: delta must be positive
  lower_bounds <- rep(1e-6, m)
  upper_bounds <- rep(10, m)  # Reasonable upper bound
  
  # Try to find feasible delta using optimization
  tryCatch({
    # Use constrOptim for constrained optimization
    opt_result <- constrOptim(
      theta = initial_delta,
      f = objective,
      grad = function(delta) 2 * (delta - requested_delta),
      ui = diag(m),  # delta >= 0 constraints
      ci = rep(0, m),  # delta >= 0
      method = "BFGS",
      control = list(maxit = max_iter)
    )
    
    suggested_delta <- opt_result$par

    # Check if the suggested delta is actually feasible (need complete check)
    feasibility_check <- analyze_delta_feasibility_all_obs(propensity_scores, suggested_delta, tol,
                                                           early_stop = FALSE)
    
    if (feasibility_check$n_feasible == n) {
      # Successfully found feasible delta
      return(list(
        suggested_delta = suggested_delta,
        feasibility_margin = min(sapply(feasibility_check$polytope_info, function(x) x$distance_to_boundary)),
        n_obs_feasible = n,
        optimization_info = list(
          converged = opt_result$convergence == 0,
          iterations = opt_result$iterations,
          value = opt_result$value
        )
      ))
    }
  }, error = function(e) {
    # Optimization failed, try alternative approach
  })
  
  # Fallback: try to find a feasible delta by scaling down the requested delta
  scale_factors <- seq(0.1, 1.0, by = 0.1)
  
  for (scale in scale_factors) {
    candidate_delta <- 1 + scale * (requested_delta - 1)  # Scale from 1 towards requested

    # Check feasibility (need complete check to verify all observations)
    feasibility_check <- analyze_delta_feasibility_all_obs(propensity_scores, candidate_delta, tol,
                                                           early_stop = FALSE)
    
    if (feasibility_check$n_feasible == n) {
      return(list(
        suggested_delta = candidate_delta,
        feasibility_margin = min(sapply(feasibility_check$polytope_info, function(x) x$distance_to_boundary)),
        n_obs_feasible = n,
        optimization_info = list(
          converged = TRUE,
          iterations = 0,
          value = sum((candidate_delta - requested_delta)^2),
          method = "scaling"
        )
      ))
    }
  }
  
  # If no feasible delta found, return the original with warning
  warning("No feasible delta found. Returning requested delta with feasibility issues.")
  
  return(list(
    suggested_delta = requested_delta,
    feasibility_margin = 0,
    n_obs_feasible = 0,
    optimization_info = list(
      converged = FALSE,
      iterations = 0,
      value = Inf,
      method = "none"
    )
  ))
}

#' Validate that required columns are present in the data
#'
#' Asserts that all named treatment, covariate, outcome, and (optionally)
#' eligibility columns exist in `data`, failing loudly with a message that lists
#' every missing column. This is the single validation point used by the
#' estimation entry points, which require explicit column specifications rather
#' than auto-detecting columns by name pattern.
#'
#' Validation checks for the *presence* of the required columns only; it does not
#' forbid extra columns. The data legitimately carries id, missingness, date, and
#' eligibility columns (and may carry treatments outside the analysis set); the
#' estimators use exactly the named columns and ignore the rest.
#'
#' @param data Data frame to validate.
#' @param treatment_cols Character vector of treatment column names.
#' @param covariate_cols Character vector of covariate column names.
#' @param outcome_col Character scalar, outcome column name. Default "Y".
#' @param eligibility_cols Character vector of eligibility column names, or NULL
#'   when eligibility is not used.
#' @return Invisibly `TRUE`. Called for its side effect (an error on failure).
#' @keywords internal
validate_columns <- function(data,
                             treatment_cols,
                             covariate_cols,
                             outcome_col = "Y",
                             eligibility_cols = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  # Each column spec must be a non-empty character vector (catches NULL, NA,
  # numeric mistakes, and the "forgot to pass it" case). eligibility may be NULL.
  check_spec <- function(x, nm, allow_null = FALSE) {
    if (allow_null && is.null(x)) return(invisible(NULL))
    if (!is.character(x) || length(x) == 0 || anyNA(x)) {
      stop("`", nm, "` must be a non-empty character vector of column names.",
           call. = FALSE)
    }
    invisible(NULL)
  }
  check_spec(treatment_cols, "treatment_cols")
  check_spec(covariate_cols, "covariate_cols")
  check_spec(outcome_col, "outcome_col")
  check_spec(eligibility_cols, "eligibility_cols", allow_null = TRUE)

  requested <- c(treatment_cols, covariate_cols, outcome_col, eligibility_cols)
  missing_cols <- setdiff(requested, names(data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".\n",
         "Available columns: ", paste(names(data), collapse = ", "), ".",
         call. = FALSE)
  }

  invisible(TRUE)
}
