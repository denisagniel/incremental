#' Construct Wald-type confidence intervals using EIF variance
#'
#' Computes confidence intervals for the counterfactual mean under incremental
#' distribution using the asymptotic variance from the efficient influence function.
#'
#' @param estimate Numeric, point estimate
#' @param eif_values Numeric vector, efficient influence function values
#' @param confidence_level Numeric, confidence level (default 0.95)
#' @param method Character, method for variance estimation ("asymptotic", "bootstrap")
#' @return List containing:
#'   - estimate: point estimate
#'   - se: standard error
#'   - ci_lower: lower confidence bound
#'   - ci_upper: upper confidence bound
#'   - confidence_level: confidence level used
#' @export
#' @examples
#' # Example usage (typically called internally by estimate_qdelta_mean)
#' eif_vals <- rnorm(100, 0, 1)
#' ci_result <- inference_qdelta_mean(estimate = 0.5, eif_values = eif_vals)
#' print(ci_result)
inference_qdelta_mean <- function(estimate, eif_values, confidence_level = 0.95,
                                method = "asymptotic") {

  # Input validation
  if (!is.numeric(eif_values) || !is.finite(estimate)) {
    stop("eif_values must be numeric and estimate must be finite", call. = FALSE)
  }
  if (any(!is.finite(eif_values))) {
    stop("eif_values contains NA, NaN, or Inf values", call. = FALSE)
  }
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be in (0, 1)", call. = FALSE)
  }

  n <- length(eif_values)
  if (n == 0) {
    stop("eif_values cannot be empty", call. = FALSE)
  }
  
  if (method == "asymptotic") {
    # Asymptotic variance from EIF
    eif_var <- mean(eif_values^2)
    if (eif_var < .Machine$double.eps) {
      warning("EIF variance is effectively zero; standard error set to 0. ",
              "Confidence intervals may be degenerate.", call. = FALSE)
      se <- 0
    } else {
      se <- sqrt(eif_var / n)
    }

  } else if (method == "bootstrap") {
    # Bootstrap variance (would need bootstrap samples)
    stop("Bootstrap method requires bootstrap samples - use bootstrap_qdelta_mean instead")
    
  } else {
    stop("method must be 'asymptotic'")
  }
  
  # Compute confidence interval
  alpha <- 1 - confidence_level
  z_alpha <- qnorm(1 - alpha/2)
  ci_lower <- estimate - z_alpha * se
  ci_upper <- estimate + z_alpha * se
  
  return(list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    confidence_level = confidence_level,
    method = method
  ))
}

#' Bootstrap inference for counterfactual mean under incremental distribution
#'
#' Performs bootstrap inference by resampling the data and re-estimating
#' the counterfactual mean multiple times.
#'
#' @param data tibble, observed data
#' @param delta Numeric vector, odds scaling factors
#' @param divergence Character, divergence type
#' @param estimator Character, estimation method
#' @param n_bootstrap Integer, number of bootstrap samples
#' @param confidence_level Numeric, confidence level
#' @param ... Additional arguments passed to estimate_qdelta_mean
#' @return List containing:
#'   - bootstrap_estimates: vector of bootstrap estimates
#'   - estimate: original point estimate
#'   - se: bootstrap standard error
#'   - ci_lower: lower confidence bound
#'   - ci_upper: upper confidence bound
#'   - confidence_level: confidence level used
#'
#' @details
#' Random number generation: The caller is responsible for setting the random
#' seed before calling this function (e.g., set.seed(123)) to ensure reproducibility.
#'
#' @export
#' @examples
#' # Generate example data
#' sim_data <- simulate_vector_treatment_data(n = 200, m = 2, p = 2, seed = 123)
#'
#' # Bootstrap inference (caller controls RNG)
#' set.seed(456)
#' boot_result <- bootstrap_qdelta_mean(sim_data$data, delta = c(1.5, 2.0),
#'                                   treatment_cols = c("A1", "A2"),
#'                                   covariate_cols = c("X1", "X2"),
#'                                   n_bootstrap = 100)
#' print(boot_result$se)
#' print(boot_result$ci_lower)
#' print(boot_result$ci_upper)
bootstrap_qdelta_mean <- function(data, delta, divergence = "kl", estimator = "one_step",
                                n_bootstrap = 1000, confidence_level = 0.95,
                                ...) {
  
  # Get original estimate
  original_result <- estimate_qdelta_mean(data, delta, divergence = divergence,
                                        estimator = estimator, ...)
  original_estimate <- original_result$estimate
  
  # Bootstrap estimates
  bootstrap_estimates <- numeric(n_bootstrap)
  
  # Track convergence across bootstrap samples
  boot_convergence_rates <- numeric(n_bootstrap)

  for (b in 1:n_bootstrap) {
    # Bootstrap sample
    boot_idx <- sample(1:nrow(data), nrow(data), replace = TRUE)
    boot_data <- data[boot_idx, ]

    # Re-estimate (simplified - in practice would need to re-estimate nuisance parameters)
    tryCatch({
      boot_result <- estimate_qdelta_mean(boot_data, delta, divergence = divergence,
                                        estimator = estimator, ...)
      bootstrap_estimates[b] <- boot_result$estimate

      # Track Q_delta convergence for this bootstrap sample
      n_converged <- sum(boot_result$convergence_info$converged)
      n_obs <- nrow(boot_data)
      boot_convergence_rates[b] <- n_converged / n_obs

      # Warn if many observations failed to converge
      if (boot_convergence_rates[b] < 0.9) {
        warning("Bootstrap sample ", b, ": Q_delta optimization converged for only ",
                round(boot_convergence_rates[b] * 100, 1), "% of observations. ",
                "Results may be unreliable.", call. = FALSE)
      }

    }, error = function(e) {
      # If bootstrap fails, use NA (don't create misleading results)
      warning("Bootstrap sample ", b, " failed: ", conditionMessage(e))
      bootstrap_estimates[b] <<- NA_real_
      boot_convergence_rates[b] <<- NA_real_
    })
  }
  
  # Compute bootstrap statistics
  se <- sd(bootstrap_estimates, na.rm = TRUE)
  alpha <- 1 - confidence_level
  ci_lower <- quantile(bootstrap_estimates, alpha/2, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_estimates, 1 - alpha/2, na.rm = TRUE)
  
  return(list(
    bootstrap_estimates = bootstrap_estimates,
    estimate = original_estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    confidence_level = confidence_level,
    n_bootstrap = n_bootstrap
  ))
}

#' Hypothesis testing for counterfactual mean under incremental distribution
#'
#' Performs hypothesis tests for the counterfactual mean parameter,
#' including tests for zero effect and comparisons between different delta values.
#'
#' @param estimate Numeric, point estimate
#' @param se Numeric, standard error
#' @param null_value Numeric, null hypothesis value (default 0)
#' @param alternative Character, alternative hypothesis ("two.sided", "greater", "less")
#' @return List containing:
#'   - test_statistic: z-statistic
#'   - p_value: p-value
#'   - null_value: null hypothesis value
#'   - alternative: alternative hypothesis
#' @export
#' @examples
#' # Test H0: psi_delta = 0
#' test_result <- hypothesis_test_qdelta_mean(estimate = 0.5, se = 0.1)
#' print(test_result$p_value)
hypothesis_test_qdelta_mean <- function(estimate, se, null_value = 0, 
                                      alternative = "two.sided") {
  
  # Compute test statistic
  test_statistic <- (estimate - null_value) / se
  
  # Compute p-value
  if (alternative == "two.sided") {
    p_value <- 2 * (1 - pnorm(abs(test_statistic)))
  } else if (alternative == "greater") {
    p_value <- 1 - pnorm(test_statistic)
  } else if (alternative == "less") {
    p_value <- pnorm(test_statistic)
  } else {
    stop("alternative must be one of: 'two.sided', 'greater', 'less'")
  }
  
  return(list(
    test_statistic = test_statistic,
    p_value = p_value,
    null_value = null_value,
    alternative = alternative
  ))
}

#' Sensitivity analysis for counterfactual mean under incremental distribution
#'
#' Performs sensitivity analysis by varying the delta values and examining
#' how the estimates change.
#'
#' @param data tibble, observed data
#' @param delta_grid List of numeric vectors, grid of delta values to test
#' @param divergence Character, divergence type
#' @param estimator Character, estimation method
#' @param confidence_level Numeric, confidence level
#' @param ... Additional arguments passed to estimate_qdelta_mean
#' @return List containing:
#'   - estimates: vector of estimates for each delta combination
#'   - confidence_intervals: matrix of confidence intervals
#'   - delta_grid: grid of delta values used
#'   - summary: summary statistics
#' @export
#' @examples
#' # Sensitivity analysis with different delta values
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 2, p = 2, seed = 123)
#' 
#' # Create delta grid
#' delta_grid <- list(
#'   c(1.2, 1.2), c(1.5, 1.5), c(2.0, 2.0),
#'   c(1.2, 1.5), c(1.5, 2.0), c(2.0, 1.2)
#' )
#' 
#' sens_result <- sensitivity_analysis(sim_data$data, delta_grid,
#'                                   treatment_cols = c("A1", "A2"),
#'                                   covariate_cols = c("X1", "X2"))
#' print(sens_result$summary)
sensitivity_analysis <- function(data, delta_grid, divergence = "kl", 
                                estimator = "one_step", confidence_level = 0.95, ...) {
  
  n_delta <- length(delta_grid)
  estimates <- numeric(n_delta)
  ci_lower <- numeric(n_delta)
  ci_upper <- numeric(n_delta)
  p_values <- numeric(n_delta)
  
  for (i in 1:n_delta) {
    delta_i <- delta_grid[[i]]
    
    tryCatch({
      result <- estimate_qdelta_mean(data, delta_i, divergence = divergence,
                                   estimator = estimator, ...)
      estimates[i] <- result$estimate
      ci_lower[i] <- result$ci_lower
      ci_upper[i] <- result$ci_upper
      p_values[i] <- result$p_value
    }, error = function(e) {
      warning("sensitivity_analysis(): delta_grid[[", i, "]] = c(",
              paste(delta_i, collapse = ", "), ") failed: ", conditionMessage(e),
              call. = FALSE)
      estimates[i] <<- NA
      ci_lower[i] <<- NA
      ci_upper[i] <<- NA
      p_values[i] <<- NA
    })
  }
  
  # Create summary
  summary_df <- data.frame(
    delta_1 = sapply(delta_grid, function(x) x[1]),
    delta_2 = sapply(delta_grid, function(x) x[2]),
    estimate = estimates,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_values,
    significant = p_values < (1 - confidence_level)
  )
  
  return(list(
    estimates = estimates,
    confidence_intervals = cbind(ci_lower, ci_upper),
    delta_grid = delta_grid,
    summary = summary_df,
    confidence_level = confidence_level
  ))
}
