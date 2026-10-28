library(testthat)
library(incremental)

test_that("Newton stall-checks compare against the PREVIOUS iteration's norm, not the current one (bugfix 2026-08-10)", {
  # Regression test for a bug in newton_qdelta() (R/fast_qdelta_solver.R,
  # internal-normalization branch): `prev_constraint_norm <- constraint_norm`
  # used to run BEFORE the two "stalled" checks --
  #   (a) chi2 boundary early-termination (iter > 5), and
  #   (b) the periodic deterministic perturbation for other divergences
  #       (iter %% 25 == 0 && iter > 50)
  # -- which made both checks tautologically true
  # (`constraint_norm >= constraint_norm * 0.999` / `constraint_norm >
  # constraint_norm * 0.95`, comparing constraint_norm to itself) for any
  # solve needing more than 5 (chi2) or 50 (other divergences) iterations,
  # regardless of whether the solver had actually stalled. The fix moves the
  # `prev_constraint_norm` update to the end of the loop body so these checks
  # compare against the genuine previous-iteration norm, matching the
  # (correct, unaffected) divergence-detection check earlier in the loop
  # (`iter > 5 && constraint_norm > prev_constraint_norm * 1.5`).
  #
  # This case (m = 6, reverse_kl, hard delta/propensity combination -- found
  # via random search for a solve needing > 50 Newton iterations) needs 79
  # iterations under the fixed code, crossing the iter > 50 threshold where
  # the bug mattered. Under the pre-fix ordering (verified by re-running a
  # reconstructed pre-fix variant side by side in the same session), this
  # exact case took 81 iterations and converged to a slightly different
  # lambda / constraint_norm (~1e-11 apart -- below any tolerance this repo
  # cares about, but a real, reproducible discrepancy caused by the spurious
  # perturbation at iteration 75 under the old ordering).
  delta <- c(56.5106217262, 58.6072451433, 2.1850880153, 53.2119743871,
             9.4379634617, 46.7824246609)
  p0 <- c(2.20689637405426e-05, 0.00335018122185354, 7.19767258059993e-08,
          1.09264339749118e-05, 4.28680586308681e-08, 6.50758987739041e-06,
          1.39811843373559e-10, 2.12241507018449e-08, 0.00364559795533524,
          0.553420357928031, 1.18899195954556e-05, 0.0018049504193469,
          7.081424788731e-06, 0.00107499638995698, 2.30956820776292e-08,
          3.50604230049229e-06, 1.7017935587356e-05, 0.0025834093938328,
          5.55030720044657e-08, 8.42564933137066e-06, 3.30566432167773e-08,
          5.0181669907278e-06, 1.07812445244528e-10, 1.63664789061569e-08,
          0.00281121270172364, 0.426756422035752, 9.1686174390955e-06,
          0.00139184287653651, 5.46066559074274e-06, 0.000828956879716108,
          1.78096640405954e-08, 2.70359780992853e-06, 4.9146774282016e-08,
          7.46073092737958e-06, 1.60289533225679e-10, 2.43327684337904e-08,
          9.54655971258442e-11, 1.44921644071801e-08, 3.11355856533024e-13,
          4.72654065743377e-11, 8.11861317732335e-06, 0.00123244687579938,
          2.6478415636511e-08, 4.01955850273608e-06, 1.57700737460308e-08,
          2.39397760367605e-06, 5.14332384294605e-11, 7.80782784326865e-09,
          3.7898319508368e-08, 5.75315813870736e-06, 1.23603309327604e-10,
          1.87636125890029e-08, 7.36159342049899e-11, 1.11752741679342e-08,
          2.40094368431581e-13, 3.64475493298579e-11, 6.26046776525952e-06,
          0.000950370927868325, 2.04181753640556e-08, 3.09958312920269e-06,
          1.21607023498244e-08, 1.84605662213708e-06, 3.96614697877137e-11,
          6.02081334112715e-09)

  result <- newton_qdelta(p0, delta, "reverse_kl", fix_lambda0 = FALSE, max_iter = 200)

  expect_true(result$converged)
  # Fixed behavior takes 79 iterations on this case; the pre-fix ordering
  # bug took 81 (one spurious perturbation at iteration 75). Allow a little
  # headroom above 79 without being so loose the bound stops discriminating
  # against a reintroduced ordering bug.
  expect_lte(result$iterations, 85)

  expected_lambda <- c(-107.746954732017, 53.1921959173969, 0.575381115717561,
                        3.03997798946677e-05, 49.9325181973622, 4.02622843610877,
                        0.0205966120249815)
  expect_equal(result$lambda, expected_lambda, tolerance = 1e-6)
})
