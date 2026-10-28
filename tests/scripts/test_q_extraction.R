library(incremental)

set.seed(2026)
n <- 10  # Small sample
m <- 8
p <- 10

data_test <- data.frame(matrix(rnorm(n * p), nrow = n, ncol = p))
colnames(data_test) <- paste0("X", 1:p)

for (j in 1:m) {
  relevant_covs <- sample(1:p, size = min(2, p), replace = FALSE)
  linear_predictor <- rowSums(data_test[, relevant_covs, drop = FALSE]) / sqrt(2)
  prob_j <- plogis(linear_predictor - 1)
  data_test[[paste0("A", j)]] <- rbinom(n, 1, prob_j)
}
data_test$Y <- rnorm(n, mean = 0.5 + rowSums(data_test[, paste0("A", 1:m)]) * 0.1, sd = 1)

ps_est <- estimate_propensity_scores(
  data = data_test,
  method = "binary_per_combination",
  binary_method = "xgboost",
  cross_fit = FALSE,
  seed = 2026
)

delta_test <- rep(1.2, m)

# Compute Q_delta for all observations
qdelta_results <- vector("list", n)
for (i in 1:n) {
  qdelta_results[[i]] <- compute_qdelta_per_obs(
    ps_est$propensity_scores[i, ],
    delta_test,
    divergence = "kl",
    tol = 1e-8
  )
}

# Check first few results
cat("First observation q_delta:\n")
cat("  Length:", length(qdelta_results[[1]]$q_delta), "\n")
cat("  Sum:", sum(qdelta_results[[1]]$q_delta), "\n")
cat("  First 10 values:", paste(round(head(qdelta_results[[1]]$q_delta, 10), 4), collapse = ", "), "\n\n")

# Extract Q matrix (like in eif.R line 136)
Q <- as.matrix(do.call(rbind, lapply(qdelta_results, function(x) x$q_delta)))

cat("Q matrix after extraction:\n")
cat("  Dimensions:", nrow(Q), "×", ncol(Q), "\n")
cat("  Row sums (first 5):", paste(round(head(rowSums(Q), 5), 4), collapse = ", "), "\n")
cat("  Any NA?", any(is.na(Q)), "\n")
cat("  Any Inf?", any(is.infinite(Q)), "\n")
cat("  Min value:", min(Q), "\n")
cat("  Max value:", max(Q), "\n")
