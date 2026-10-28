# Fidelity test for the Omega eigenvalue diagnostic
# (application/analysis/12_omega_eigenvalue_diagnostic.R).
#
# The diagnostic reconstructs the standard estimator's per-observation Omega
# using compute_omega_reduced() and reports its eigenvalue spectrum. This test
# asserts that the reconstructed Omega is IDENTICAL to the matrix
# calculate_eif_lambda() forms and inverts (R/eif.R:771-775), so the diagnostic
# characterizes the estimator's actual conditioning -- not an approximation.

test_that("compute_omega_reduced reproduces calculate_eif_lambda's Omega (KL)", {
  set.seed(20260706)
  m <- 3L
  a_bar <- cbind(1, generate_binary_vectors(m))  # 2^m x (m+1)

  for (trial in 1:10) {
    pi_i <- runif(2^m); pi_i <- as.vector(pi_i / sum(pi_i))
    delta <- rep(1.5, m)
    qd <- compute_qdelta_per_obs(pi_i, delta, "kl", tol = 1e-8)

    # Estimator's construction (calculate_eif_lambda, KL: h'(x) = exp(x))
    lin <- as.vector(a_bar %*% qd$lambda)
    weights <- pi_i * exp(lin)
    omega_eif <- t(a_bar) %*% diag(weights) %*% a_bar

    # Diagnostic's construction
    omega_diag <- compute_omega_reduced(pi_i, qd$lambda, a_bar, "kl")

    expect_equal(omega_diag, omega_eif, tolerance = 1e-10)
  }
})

test_that("active-constraint subsetting yields a positive-definite Omega (interior)", {
  set.seed(1)
  m <- 3L
  a_bar <- cbind(1, generate_binary_vectors(m))
  pi_i <- runif(2^m); pi_i <- as.vector(pi_i / sum(pi_i))  # interior
  qd <- compute_qdelta_per_obs(pi_i, rep(1.5, m), "kl", tol = 1e-8)

  omega_full <- compute_omega_reduced(pi_i, qd$lambda, a_bar, "kl")
  E_A <- compute_marginals(pi_i, a_bar[, -1, drop = FALSE], m = m)
  active <- c(1L, find_active_constraints(E_A, tol = 1e-8) + 1L)
  omega_active <- omega_full[active, active, drop = FALSE]

  ev <- eigen(omega_active, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))                 # positive definite under interiority
  expect_true(is.finite(max(ev) / min(ev)))  # bounded condition number
})
