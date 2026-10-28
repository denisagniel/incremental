#' Minimal self-contained DGP for oracle-estimation tests
#'
#' Reproduces "DGP-1: Linear-Linear (Simple)" from the paper's simulation
#' study: linear logistic propensities, linear outcome model, no
#' interactions. The full 15-DGP simulation library lives in the research
#' monorepo (sims/dgps/), not in this package, so this helper duplicates
#' just the DGP-1 generative logic to keep the package test suite
#' self-contained.
#'
#' @param n Sample size.
#' @param m Number of binary treatments.
#' @param p Number of covariates.
#' @param seed Random seed.
#' @return List with `data` (tibble), `true_propensities` (n x 2^m matrix),
#'   and `true_outcomes` (n x 2^m matrix), matching the format expected by
#'   `estimate_qdelta_mean(oracle_propensities = , oracle_outcomes = )`.
generate_dgp1_oracle_data <- function(n = 100, m = 2, p = 2, seed = 123) {
  set.seed(seed)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)

  linear_pred <- 0.5 * X[, 1] + if (p >= 2) 0.3 * X[, 2] else 0
  prob <- stats::plogis(linear_pred)

  A <- matrix(0, nrow = n, ncol = m)
  colnames(A) <- paste0("A", 1:m)
  for (k in 1:m) A[, k] <- rbinom(n, 1, prob)

  treatment_effects <- c(0.5, 0.3)[1:m]
  covariate_effects <- c(0.3, 0.4)[1:p]

  Y <- 1.0 + X %*% covariate_effects + A %*% treatment_effects + rnorm(n, 0, 1.0)

  data_list <- c(
    setNames(lapply(1:p, function(i) X[, i]), paste0("X", 1:p)),
    setNames(lapply(1:m, function(i) A[, i]), paste0("A", 1:m)),
    list(Y = as.numeric(Y))
  )
  data <- tibble::as_tibble(data_list)

  all_combos <- expand.grid(lapply(1:m, function(k) c(0, 1)))
  colnames(all_combos) <- paste0("A", 1:m)
  combo_names <- apply(all_combos, 1, function(x) paste0(paste0("A", 1:m, "=", x), collapse = "_"))

  true_propensities <- matrix(NA_real_, nrow = n, ncol = nrow(all_combos))
  true_outcomes <- matrix(NA_real_, nrow = n, ncol = nrow(all_combos))
  for (i in seq_len(nrow(all_combos))) {
    combo <- as.numeric(all_combos[i, ])
    prob_combo <- rep(1, n)
    for (k in 1:m) prob_combo <- prob_combo * ifelse(combo[k] == 1, prob, 1 - prob)
    true_propensities[, i] <- prob_combo
    true_outcomes[, i] <- 1.0 + X %*% covariate_effects + sum(combo * treatment_effects)
  }
  colnames(true_propensities) <- combo_names
  colnames(true_outcomes) <- combo_names

  list(data = data, true_propensities = true_propensities, true_outcomes = true_outcomes)
}
