# Regression tests for the pooled stratified SE (.combine_patterns).
# Weighting the within-pattern influence value by w_k alone understated the SE by
# roughly w_k^2 per pattern (e.g. ~$74 vs ~$340 for total cost at delta = 1).
# When phi_i = Y_i (as at delta = 1) the pooled estimator is the plain mean of Y,
# so the SE must equal the population-variance (mean(x^2), not var()'s n-1
# correction -- see R/inference.R's EIF SE convention) analog of sd(Y) / sqrt(n):
# sqrt(mean((Y - mean(Y))^2) / n).

test_that("pooled SE equals the population-variance analog of sd(y)/sqrt(n) when phi = y (full coverage)", {
  set.seed(20261026)
  n <- 2000
  k <- sample(1:4, n, replace = TRUE, prob = c(.5, .3, .15, .05))
  y <- rnorm(n, mean = c(0, 5, -3, 10)[k], sd = c(1, 4, 2, 8)[k])
  est <- vapply(1:4, function(j) mean(y[k == j]), numeric(1))
  prev <- tabulate(k, 4) / n
  eif <- y - est[k]  # within-pattern centered and unscaled, as the fits return it

  out <- .combine_patterns(est, prev, obs_eif = eif, obs_pattern_pos = k, n_full = n)

  expect_equal(out$estimate, mean(y), tolerance = 1e-12)
  expect_equal(out$se, sqrt(mean((y - mean(y))^2) / n), tolerance = 1e-10)
})

test_that("uncentered per-pattern influence values give the same SE", {
  set.seed(20261027)
  n <- 1500
  k <- sample(1:3, n, replace = TRUE)
  y <- rnorm(n, c(0, 2, 7)[k], c(1, 3, 5)[k])
  est <- vapply(1:3, function(j) mean(y[k == j]), numeric(1))
  prev <- tabulate(k, 3) / n
  a <- .combine_patterns(est, prev, y - est[k], k, n)
  b <- .combine_patterns(est, prev, y, k, n)  # uncentered phi = y
  expect_equal(a$se, b$se, tolerance = 1e-10)
})

test_that("subset of patterns: SE is that of the mean over included observations", {
  set.seed(20261028)
  n <- 2000
  k <- sample(1:3, n, replace = TRUE)
  y <- rnorm(n, c(0, 2, 7)[k], c(1, 3, 5)[k])
  est <- vapply(1:3, function(j) mean(y[k == j]), numeric(1))
  keep <- 1:2
  pos <- ifelse(k %in% keep, k, NA_integer_)

  out <- .combine_patterns(est[keep], (tabulate(k, 3) / n)[keep],
                           obs_eif = ifelse(is.na(pos), 0, y - est[k]),
                           obs_pattern_pos = pos, n_full = n)

  yS <- y[k %in% keep]
  nS <- length(yS)
  c_hat <- nS / n
  # IF_i = (y_i - ybar_S) / c on S, 0 elsewhere; mean(IF) = 0.
  # Population-variance convention (mean(x^2), not var()'s n-1 correction --
  # see R/inference.R's EIF SE convention).
  expected <- sqrt(sum((yS - mean(yS))^2) / (c_hat^2 * n) / n)
  expect_equal(out$estimate, mean(yS), tolerance = 1e-12)
  expect_equal(out$se, expected, tolerance = 1e-10)
  expect_equal(out$coverage, c_hat, tolerance = 1e-12)
})

test_that("fixed custom weights drop the between-pattern term", {
  set.seed(20261029)
  n <- 3000
  k <- sample(1:2, n, replace = TRUE, prob = c(.6, .4))
  y <- rnorm(n, c(0, 50)[k], c(1, 2)[k])  # large between-pattern spread
  est <- vapply(1:2, function(j) mean(y[k == j]), numeric(1))
  n_k <- tabulate(k, 2)
  w <- c(.5, .5)
  out <- .combine_patterns(est, w, y - est[k], k, n, prevalence_weights = FALSE)
  expected <- sqrt(sum(w^2 * vapply(1:2, function(j) var(y[k == j]), 1) / n_k))
  expect_equal(out$se, expected, tolerance = 1e-3)
})

# --- stratified_contrast() ---------------------------------------------------
.make_contrib <- function(y, k, est_mult = 1) {
  K <- max(k); n <- length(y)
  est <- vapply(seq_len(K), function(j) mean(est_mult * y[k == j]), numeric(1))
  list(
    patterns = data.frame(
      pattern_id = seq_len(K), n_obs = tabulate(k, K),
      prevalence = tabulate(k, K) / n, estimate = est, se = NA_real_
    ),
    obs_pattern_id = k,
    obs_eif = est_mult * y - est[k],
    n_full = n,
    confidence_level = 0.95
  )
}

test_that("contrast of identical fits is exactly zero with zero SE", {
  set.seed(20261030)
  k <- sample(1:3, 800, replace = TRUE); y <- rnorm(800, c(0, 2, 5)[k])
  a <- .make_contrib(y, k)
  out <- stratified_contrast(a, a)
  expect_equal(out$estimate, 0, tolerance = 1e-12)
  expect_equal(out$se, 0, tolerance = 1e-12)
})

test_that("contrast of phi = 2y vs phi = y has the population-variance analog of SE sd(y)/sqrt(n)", {
  set.seed(20261031)
  n <- 1500
  k <- sample(1:4, n, replace = TRUE, prob = c(.4, .3, .2, .1))
  y <- rnorm(n, c(0, 3, -2, 8)[k], c(1, 2, 3, 4)[k])
  out <- stratified_contrast(.make_contrib(y, k, 1), .make_contrib(y, k, 2))
  expect_equal(out$estimate, mean(y), tolerance = 1e-12)  # mean(2y) - mean(y)
  # Population-variance convention (mean(x^2), not var()'s n-1 correction --
  # see R/inference.R's EIF SE convention).
  expect_equal(out$se, sqrt(mean((y - mean(y))^2) / n), tolerance = 1e-10)
})

test_that("contrast refuses fits on different observations or patterns", {
  set.seed(20261032)
  k <- sample(1:2, 300, replace = TRUE); y <- rnorm(300)
  a <- .make_contrib(y, k); b <- .make_contrib(y, k)
  b$obs_pattern_id <- rev(b$obs_pattern_id)
  expect_error(stratified_contrast(a, b), "same observations")
})
