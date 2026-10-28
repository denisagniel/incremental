#' Helpers for nuisance-cache tests
#'
#' Synthetic multi-pattern eligibility data small enough for fast tests but
#' structurally like the application: m treatments, several eligibility
#' patterns, treatment only possible where eligible, and several synthetic
#' "cost types" that differ ONLY in the outcome column (same cohort, same
#' covariates, same treatments) -- which is exactly what makes sharing
#' propensity fits across cost types valid.
#'
#' As in the application (`application/analysis/00_config.R`:
#' `outcome_name = "Y"`), each cost type is a data frame whose outcome column is
#' named "Y"; `cache_test_dataset()` builds it.

#' @param n Number of observations.
#' @param m Number of treatments.
#' @param p Number of covariates.
#' @param patterns Character vector of eligibility pattern strings (m chars each).
#' @param probs Pattern prevalences.
#' @param seed Seed for data generation.
#' @param n_cost_types Number of synthetic cost types.
#' @return List with `base` (data frame without outcome), `outcomes` (list of
#'   numeric vectors), `covariate_cols`, `eligibility_cols`, `cost_types`.
make_cache_test_data <- function(n = 600, m = 6, p = 3,
                                 patterns = c("111000", "110100", "101010"),
                                 probs = c(0.4, 0.35, 0.25),
                                 seed = 20260811L,
                                 n_cost_types = 2) {
  stopifnot(all(nchar(patterns) == m), length(probs) == length(patterns))

  withr::with_seed(seed, {
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- paste0("X", seq_len(p))

    pat_idx <- sample(seq_along(patterns), n, replace = TRUE, prob = probs)
    elig <- do.call(rbind, lapply(pat_idx, function(k) {
      as.integer(strsplit(patterns[k], "")[[1]])
    }))
    colnames(elig) <- paste0("elig_A", seq_len(m))

    A <- matrix(0L, nrow = n, ncol = m)
    for (k in seq_len(m)) {
      # Moderate, covariate-dependent propensities kept well away from 0/1 so
      # small deltas stay feasible for every observation.
      pk <- plogis(0.2 * X[, 1] - 0.15 * X[, min(2, p)] + 0.1 * k - 0.2)
      A[, k] <- ifelse(elig[, k] == 1, rbinom(n, 1, pk), 0L)
    }
    colnames(A) <- paste0("A", seq_len(m))

    lin <- 1 + 0.6 * X[, 1] + 0.3 * X[, min(2, p)] +
      as.numeric(A %*% seq(0.5, 1.0, length.out = m))
    outcomes <- lapply(seq_len(n_cost_types), function(j) {
      as.numeric(lin * (1 + 0.4 * (j - 1)) + rnorm(n, sd = 0.5 + 0.2 * (j - 1)))
    })
    cost_types <- paste0("ct", seq_len(n_cost_types))
    names(outcomes) <- cost_types

    list(
      base = cbind(as.data.frame(A), as.data.frame(elig), as.data.frame(X)),
      outcomes = outcomes,
      covariate_cols = paste0("X", seq_len(p)),
      eligibility_cols = paste0("elig_A", seq_len(m)),
      cost_types = cost_types
    )
  })
}

#' Build the analysis data frame for one synthetic cost type (outcome named "Y").
#'
#' @param d Output of make_cache_test_data().
#' @param cost_type Name of the cost type (element of `d$cost_types`).
#' @return Data frame with a "Y" column.
cache_test_dataset <- function(d, cost_type = d$cost_types[1]) {
  out <- d$base
  out$Y <- d$outcomes[[cost_type]]
  out
}

#' Standard cache-test estimation call (kept in one place so tests stay short).
#'
#' @param d Output of make_cache_test_data().
#' @param delta Numeric vector of length m.
#' @param cost_type Which synthetic cost type to estimate.
#' @param cache Optional nuisance cache.
#' @param tag Cache tag (defaults to cost_type when a cache is supplied).
#' @param estimator "one_step" (or "plugin"; see estimate_qdelta_mean()).
#' @param ... Passed to estimate_qdelta_mean_stratified().
run_cache_test_estimate <- function(d, delta, cost_type = d$cost_types[1],
                                    cache = NULL, tag = cost_type,
                                    estimator = "one_step", ...) {
  suppressWarnings(estimate_qdelta_mean_stratified(
    data = cache_test_dataset(d, cost_type),
    delta = delta,
    covariate_cols = d$covariate_cols,
    eligibility_cols = d$eligibility_cols,
    outcome_col = "Y",
    divergence = "kl",
    min_stratum_size = 100,
    min_realized_coverage = 0.8,
    max_patterns = 10,
    ps_method = "sequential",
    binary_method = "glm",
    outcome_method = "glm",
    estimator = estimator,
    cross_fit = TRUE,
    n_folds = 3,
    verbose = FALSE,
    nuisance_cache = cache,
    nuisance_tag = if (is.null(cache)) NULL else tag,
    ...
  ))
}
