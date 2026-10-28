library(testthat)
library(incremental)

# Regression tests for the canonical divergence_h_pair() primitive
# (R/divergence_h_functions.R), added 2026-08-10 when six independent,
# inconsistent h(x)/h'(x) implementations across fast_qdelta_solver.R,
# eif.R, eif_active.R, and estimation.R were unified into one. See
# [LEARN 2026-08-10] in R/divergence_h_functions.R and MEMORY.md.

test_that("divergence_h_pair matches analytic h/h' away from boundaries", {
  x <- c(-3, -1, -0.5, 0.2, 0.9)

  # KL: h(x) = exp(x), h'(x) = exp(x)
  hp <- divergence_h_pair(x, "kl")
  expect_equal(hp$h, exp(x), tolerance = 1e-12)
  expect_equal(hp$h_prime, exp(x), tolerance = 1e-12)
  expect_equal(hp$n_clipped, 0L)

  # reverse_kl: h(x) = -1/x, h'(x) = 1/x^2 (domain x < 0)
  xn <- c(-5, -2, -0.5, -0.01)
  hp <- divergence_h_pair(xn, "reverse_kl")
  expect_equal(hp$h, -1 / xn, tolerance = 1e-12)
  expect_equal(hp$h_prime, 1 / xn^2, tolerance = 1e-12)
  expect_equal(hp$n_clipped, 0L)

  # chi2: h(x) = x/2 + 1, h'(x) = 1/2 (constant), never clipped
  hp <- divergence_h_pair(x, "chi2")
  expect_equal(hp$h, x / 2 + 1, tolerance = 1e-12)
  expect_equal(hp$h_prime, rep(0.5, length(x)))
  expect_equal(hp$n_clipped, 0L)

  # hellinger: h(x) = 1/(1-x)^2, h'(x) = 2/(1-x)^3 (domain x < 1)
  xh <- c(-2, 0, 0.5, 0.9)
  hp <- divergence_h_pair(xh, "hellinger")
  expect_equal(hp$h, 1 / (1 - xh)^2, tolerance = 1e-12)
  expect_equal(hp$h_prime, 2 / (1 - xh)^3, tolerance = 1e-12)
  expect_equal(hp$n_clipped, 0L)
})

test_that("divergence_h_pair saturates (does not blow up or go negative) at domain boundaries", {
  # reverse_kl: as x -> 0^-, h and h' -> +Inf. Clipped h/h' must be large,
  # finite, and positive -- not the old H_VALUE_SMALL = 1e-6 (wrong
  # direction) and not the inline convention's unguarded negative/-Inf.
  hp <- divergence_h_pair(0, "reverse_kl")          # exactly at the boundary
  expect_true(is.finite(hp$h) && hp$h > 1e6)
  expect_true(is.finite(hp$h_prime) && hp$h_prime > 1e6)
  expect_gt(hp$n_clipped, 0)

  hp <- divergence_h_pair(5, "reverse_kl")          # outside the domain entirely
  expect_true(is.finite(hp$h) && hp$h > 0)          # not negative, not -Inf
  expect_true(is.finite(hp$h_prime) && hp$h_prime > 0)

  # hellinger: as x -> 1^-, h and h' -> +Inf.
  hp <- divergence_h_pair(1, "hellinger")
  expect_true(is.finite(hp$h) && hp$h > 1e6)
  expect_true(is.finite(hp$h_prime) && hp$h_prime > 1e6)
  expect_gt(hp$n_clipped, 0)

  hp <- divergence_h_pair(1.5, "hellinger")         # outside the domain entirely
  expect_true(is.finite(hp$h) && hp$h > 0)
  expect_true(is.finite(hp$h_prime) && hp$h_prime > 0)

  # kl: overflow guard only, well past exp() overflow range.
  hp <- divergence_h_pair(1000, "kl")
  expect_true(is.finite(hp$h))
  expect_gt(hp$n_clipped, 0)
})

test_that("h and h' are evaluated at the SAME clipped x (fixes the Newton-Jacobian defect)", {
  # Before the fix, the solver's Jacobian used the unclipped linear term
  # while the moment sum used the clipped one -- so, near a boundary, h and
  # h' corresponded to different arguments. Here both must come from
  # identical clipping.
  for (div in c("reverse_kl", "hellinger")) {
    x_boundary <- if (div == "reverse_kl") 0.5 else 1.5   # outside domain
    hp <- divergence_h_pair(x_boundary, div)
    # Reconstruct h' from h analytically and confirm they're consistent with
    # a SINGLE clipped argument (not e.g. h at clipped x but h' at raw x).
    if (div == "reverse_kl") {
      xs <- -1 / hp$h            # implied clipped x from h
      expect_equal(hp$h_prime, 1 / xs^2, tolerance = 1e-10)
    } else {
      d <- 1 / sqrt(hp$h)        # implied (1 - clipped x) from h
      expect_equal(hp$h_prime, 2 / d^3, tolerance = 1e-10)
    }
  }
})

test_that("divergence_h_pair works elementwise on vectors and matrices identically", {
  x_vec <- c(-2, -0.5, 0.3)
  x_mat <- matrix(x_vec, nrow = 1)
  for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
    hp_vec <- divergence_h_pair(x_vec, div)
    hp_mat <- divergence_h_pair(x_mat, div)
    expect_equal(as.numeric(hp_mat$h), hp_vec$h, tolerance = 1e-12)
    expect_equal(as.numeric(hp_mat$h_prime), hp_vec$h_prime, tolerance = 1e-12)
  }
})

test_that("solver weights and compute_H_vectorized_general agree cell-for-cell, including near boundary", {
  # Regression test for the cross-site inconsistency: before 2026-08-10, the
  # solver's compute_q_delta_weights_internal() and eif.R's
  # compute_H_vectorized_general() could disagree by orders of magnitude
  # near reverse_kl/hellinger boundaries because they used different
  # override conventions. Both now delegate to divergence_h_pair(), so they
  # must agree exactly (not just approximately) for the same lambda.
  m <- 2
  vectors <- generate_binary_vectors(m)
  a_extended <- cbind(1, vectors)

  boundary_lambdas <- list(
    reverse_kl = c(0, -0.5, -0.5),    # linear terms hit 0 for a=(0,0)
    hellinger  = c(1, -0.5, -0.5)     # linear terms hit 1 for a=(0,0)
  )

  for (div in names(boundary_lambdas)) {
    lambda <- boundary_lambdas[[div]]
    p0_active <- rep(1 / nrow(vectors), nrow(vectors))

    solver_weights <- compute_q_delta_weights_internal(lambda, p0_active, vectors, div)
    solver_h <- as.numeric(solver_weights) / p0_active

    H_result <- compute_H_vectorized_general(matrix(lambda, nrow = 1), a_extended, div)
    eif_h <- as.numeric(H_result$H_mat)

    expect_equal(solver_h, eif_h, tolerance = 1e-10,
                 info = paste("divergence:", div))
  }
})

test_that("analytical Jacobian matches numerical Jacobian in the high-curvature (but valid, unclipped) region", {
  # High-curvature sanity check: close enough to the reverse_kl/hellinger
  # boundary that h/h' are large, but well outside the 1e-12 numerical
  # safety buffer, so clipping never engages and a finite-difference
  # Jacobian (step 1e-6) is well-posed. This is the regime Newton actually
  # visits for hard-to-satisfy constraints; it must still match after
  # unifying the h-function implementations.
  m <- 2
  vectors <- generate_binary_vectors(m)
  p0_active <- rep(1 / nrow(vectors), nrow(vectors))
  target_marginals <- c(0.5, 0.5)

  for (div in c("reverse_kl", "hellinger")) {
    lambda <- if (div == "reverse_kl") c(-0.01, -0.005, -0.005) else c(0.99, -0.005, -0.005)

    analytical <- compute_jacobian_analytical_internal(lambda, p0_active, vectors, div)
    numerical <- compute_jacobian_numerical_internal(lambda, p0_active, vectors,
                                                       target_marginals, div)

    expect_equal(analytical, numerical, tolerance = 1e-3,
                 info = paste("divergence:", div))
  }
})

test_that("h and h' saturate together rather than disagreeing deep inside the invalid domain", {
  # Deep inside the invalid domain (far past the 1e-12 buffer, e.g. exactly
  # at or beyond the true boundary), the clipped argument xs is pinned to a
  # constant, so a finite-difference Jacobian correctly reports ~0 local
  # sensitivity (the composite h(clip(x)) is locally flat there) while
  # divergence_h_pair()'s h_prime reports h'(xs) (large, not 0) -- these are
  # DIFFERENT questions (derivative of h at the clipped point vs. derivative
  # of the clipped composite w.r.t. the original x) and are not expected to
  # agree. What must hold, and is the actual regression guard here, is that
  # h and h' stay mutually consistent (both large and positive, evaluated at
  # the identical clipped argument) rather than the old bug where they could
  # point in opposite directions (e.g. h -> 1e-6 while a sibling call site's
  # h' -> 1e12). Per Oracle's review this region is unreachable for any
  # accepted Newton iterate (check_constraint_satisfaction_* rejects
  # candidates long before reaching it), so this is a safety-net test, not a
  # claim about Newton's actual trajectory.
  for (div in c("reverse_kl", "hellinger")) {
    x_deep <- if (div == "reverse_kl") 0 else 1
    hp <- divergence_h_pair(x_deep, div)
    expect_true(is.finite(hp$h) && hp$h > 1e6)
    expect_true(is.finite(hp$h_prime) && hp$h_prime > 1e6)
    expect_gt(hp$n_clipped, 0)
  }
})


test_that("compute_h_function_value delegates to the canonical primitive and still fails loudly out of domain", {
  expect_equal(compute_h_function_value(c(1, 0), c(0.5, 0.3, 0.7), "kl"),
               exp(0.5 + 0.3), tolerance = 1e-12)
  expect_error(compute_h_function_value(c(1, 0), c(1, 1, 1), "reverse_kl"), "Invalid domain")
  expect_error(compute_h_function_value(c(1, 0), c(1, 1, 1), "hellinger"), "Invalid domain")
})
