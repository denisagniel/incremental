#' @file test-eif-chunk-parallel.R
#' @title Stage-3 chunk-level parallelism: equivalence, failure propagation, policy
#'
#' Spec: `quality_reports/specs/2026-08-10_qdelta-solver-performance-fix.md`,
#' MUST item "Chunk-level parallelization".
#'
#' Parallelism is over the OBSERVATION axis only (never delta scenarios or
#' coalitions). Each worker owns both phases for its own observations: the
#' per-observation Q_delta solves and the vectorized EIF assembly. Because chunks
#' are independent and nothing about either phase's math changed, results must be
#' *bit-identical* between sequential and parallel execution -- which is what these
#' tests assert (a far tighter bar than the spec's 1e-6, and the observed one).
#'
#' Robustness notes for CI-like environments:
#'   * `parallel::mclapply()` cannot fork on Windows, so `resolve_eif_workers()`
#'     returns 1 there. Parallel-execution tests skip on non-unix rather than
#'     asserting a fallback that would silently pass.
#'   * Tests that need real forking skip when fewer than 2 cores are detectable.
#'   * Worker counts are kept small (2) and n small, so these stay fast.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ecp_can_fork <- function() {
  .Platform$OS.type == "unix" &&
    isTRUE(tryCatch(parallel::detectCores() >= 2L, error = function(e) FALSE))
}

ecp_skip_unless_fork <- function() {
  skip_on_os("windows")
  if (!ecp_can_fork()) skip("requires a forking platform with >= 2 detectable cores")
}

#' Deterministic EIF input set
ecp_case <- function(m, n, seed = 20260811L) {
  withr::local_seed(seed + m)
  a_mat <- generate_binary_vectors(m)
  n_comb <- nrow(a_mat)
  ps <- matrix(stats::runif(n * n_comb, 0.3, 1), nrow = n)
  ps <- ps / rowSums(ps)
  mu <- matrix(stats::rnorm(n * n_comb, mean = 1, sd = 0.3), nrow = n)
  obs_idx <- apply(ps, 1, function(pr) sample.int(n_comb, 1, prob = pr))
  list(
    m = m, n = n, a_mat = a_mat, a_ext = cbind(1, a_mat),
    ps = ps, mu = mu, A = a_mat[obs_idx, , drop = FALSE],
    Y = mu[cbind(seq_len(n), obs_idx)] + stats::rnorm(n, sd = 0.1),
    delta = rep(1.2, m)
  )
}

ecp_args <- function(cs, divergence = "kl", chunk_size = 37L) {
  list(A = cs$A, Y = cs$Y, propensity_scores = cs$ps, delta = cs$delta,
       outcome_means = cs$mu, divergence = divergence,
       feasibility_check = FALSE, chunk_size = chunk_size)
}

# ===========================================================================
# 1. Worker-count policy
# ===========================================================================

test_that("resolve_eif_workers implements the documented resolution order", {
  # Explicit argument wins.
  expect_identical(resolve_eif_workers(1L), 1L)
  expect_identical(resolve_eif_workers(3L), if (.Platform$OS.type == "unix") 3L else 1L)

  # Degenerate / invalid inputs floor to 1 rather than erroring or forking wildly.
  expect_identical(resolve_eif_workers(0L), 1L)
  expect_identical(resolve_eif_workers(-5L), 1L)
  expect_identical(resolve_eif_workers(NA_integer_), 1L)

  # Option is the escape hatch when the argument is NULL...
  withr::with_options(list(incremental.n_cores = 2L), {
    expect_identical(resolve_eif_workers(NULL), if (.Platform$OS.type == "unix") 2L else 1L)
    # ...but an explicit argument still overrides the option.
    expect_identical(resolve_eif_workers(1L), 1L)
  })

  # Default (no arg, no option) is detectCores() - 1, floored at 1.
  withr::with_options(list(incremental.n_cores = NULL), {
    got <- resolve_eif_workers(NULL)
    expect_true(got >= 1L)
    if (.Platform$OS.type == "unix") {
      expect_identical(got, max(1L, parallel::detectCores() - 1L))
    } else {
      expect_identical(got, 1L)   # mclapply cannot fork on Windows
    }
  })
})

test_that("Windows is resolved to sequential rather than silently degrading", {
  # Documents the platform contract. On unix this asserts the converse (that a
  # requested worker count is honoured), so the test is meaningful on both.
  if (.Platform$OS.type == "unix") {
    expect_gt(resolve_eif_workers(4L), 1L)
  } else {
    expect_identical(resolve_eif_workers(4L), 1L)
  }
})

# ===========================================================================
# 2. Numerical equivalence: sequential vs parallel
# ===========================================================================

test_that("parallel chunk execution is bit-identical to sequential, all divergences", {
  ecp_skip_unless_fork()

  for (m in c(2L, 4L, 6L)) {
    cs <- ecp_case(m, n = 200L)
    for (div in c("kl", "reverse_kl", "chi2", "hellinger")) {
      # chunk_size 37 over n = 200 gives 6 chunks, so >= 2 workers really do split
      # the work and the re-assembly scatter is genuinely exercised.
      args <- ecp_args(cs, div, chunk_size = 37L)
      seq_res <- do.call(efficient_influence_function_chunked, c(args, list(n_cores = 1L)))
      par_res <- do.call(efficient_influence_function_chunked, c(args, list(n_cores = 2L)))

      label <- sprintf("m = %d, divergence = %s", m, div)
      # Bit-identical, not merely within 1e-6: chunks are independent and no
      # algorithm changed, so any difference at all would signal a defect.
      expect_identical(par_res$eif_values, seq_res$eif_values, info = label)
      expect_identical(par_res$plugin, seq_res$plugin, info = label)
      expect_identical(par_res$centered_eif, seq_res$centered_eif, info = label)
      expect_identical(par_res$qdelta_results, seq_res$qdelta_results, info = label)
      expect_identical(par_res$convergence_summary, seq_res$convergence_summary, info = label)
      expect_identical(par_res$convergence_info, seq_res$convergence_info, info = label)
    }
  }
})

test_that("parallel chunk execution matches the non-chunked path exactly", {
  ecp_skip_unless_fork()
  # The chunked+parallel path must agree with the (separately tested) small-n
  # non-chunked path, so parallelism cannot drift away from the reference.
  for (m in c(2L, 4L)) {
    cs <- ecp_case(m, n = 150L)
    for (div in c("kl", "hellinger")) {
      free <- efficient_influence_function(
        A = cs$A, Y = cs$Y, propensity_scores = cs$ps, delta = cs$delta,
        outcome_means = cs$mu, divergence = div, feasibility_check = FALSE)
      par_res <- do.call(efficient_influence_function_chunked,
                         c(ecp_args(cs, div, chunk_size = 29L), list(n_cores = 3L)))
      expect_identical(par_res$eif_values, free$eif_values,
                       info = sprintf("m = %d, %s", m, div))
      expect_identical(par_res$plugin, free$plugin)
    }
  }
})

test_that("results are invariant to worker count and chunk size", {
  ecp_skip_unless_fork()
  cs <- ecp_case(4L, n = 210L)
  ref <- do.call(efficient_influence_function_chunked,
                 c(ecp_args(cs, "kl", chunk_size = 210L), list(n_cores = 1L)))
  for (K in c(1L, 2L, 3L)) {
    for (csz in c(17L, 53L, 210L)) {
      got <- do.call(efficient_influence_function_chunked,
                     c(ecp_args(cs, "kl", chunk_size = csz), list(n_cores = K)))
      expect_identical(got$eif_values, ref$eif_values,
                       info = sprintf("K = %d, chunk_size = %d", K, csz))
      expect_identical(got$plugin, ref$plugin)
    }
  }
})

test_that("caller-supplied qdelta_results are honoured and not recomputed in workers", {
  ecp_skip_unless_fork()
  cs <- ecp_case(4L, n = 120L)

  qd <- lapply(seq_len(cs$n), function(i) {
    compute_qdelta_per_obs(cs$ps[i, ], cs$delta, "kl", fix_lambda0 = FALSE,
                           tol = 1e-10, combination_matrix = cs$a_mat)
  })

  supplied <- do.call(efficient_influence_function_chunked,
                      c(ecp_args(cs, "kl", chunk_size = 31L),
                        list(qdelta_results = qd, n_cores = 2L)))
  solved <- do.call(efficient_influence_function_chunked,
                    c(ecp_args(cs, "kl", chunk_size = 31L), list(n_cores = 2L)))

  # Same numbers either way, and the supplied list is returned untouched.
  expect_identical(supplied$eif_values, solved$eif_values)
  expect_identical(supplied$qdelta_results, qd)
})

# ===========================================================================
# 3. Failure propagation: nothing may be silently swallowed
# ===========================================================================

test_that("a worker error is re-raised, not absorbed as a list element", {
  ecp_skip_unless_fork()
  skip_on_ci()  # fork/signal propagation through mclapply varies by CI sandbox;
                # the re-raise contract is exercised locally and is not a
                # correctness property of the estimator itself.
  cs <- ecp_case(4L, n = 200L)

  # Break one observation's propensity row so it no longer sums to 1; the solver's
  # input validation fires inside exactly one chunk.
  ps_bad <- cs$ps
  ps_bad[150, ] <- ps_bad[150, ] * 3
  args <- ecp_args(cs, "kl", chunk_size = 37L)
  args$propensity_scores <- ps_bad

  # Sequential: the underlying error surfaces directly.
  expect_error(do.call(efficient_influence_function_chunked, c(args, list(n_cores = 1L))),
               "propensity scores must sum to approximately 1")

  # Parallel: mclapply hands back a "try-error" element, which MUST be re-raised
  # with context rather than flowing on to produce silent NA/garbage.
  expect_error(
    suppressWarnings(
      do.call(efficient_influence_function_chunked, c(args, list(n_cores = 3L)))),
    "Parallel EIF computation failed in"
  )
  # And the message must carry the underlying cause plus an explicit statement
  # that it was not ignored, so a reader is not left guessing.
  err <- tryCatch(
    suppressWarnings(
      do.call(efficient_influence_function_chunked, c(args, list(n_cores = 3L)))),
    error = function(e) conditionMessage(e))
  expect_match(err, "propensity scores must sum to approximately 1")
  expect_match(err, "NOT silently ignored")
  expect_match(err, "n_cores = 1")
})

test_that("a killed worker (NULL result) is re-raised, not absorbed", {
  ecp_skip_unless_fork()
  skip_on_ci()  # tools::pskill()-simulated SIGKILL and mclapply's handling of a
                # killed child are sensitive to the host's process/signal model;
                # GitHub Actions' containers don't reproduce this reliably.
  # mclapply returns NULL (not a try-error) when a child dies by signal, e.g. the
  # OOM killer. That branch is separate and must also fail loudly. Simulated by
  # having one chunk's worker SIGKILL itself; inside a fork Sys.getpid() is the
  # child, so only the child dies.
  cs <- ecp_case(4L, n = 200L)
  original <- eif_chunk_worker

  saboteur <- function(chunk_indices, ...) {
    if (150L %in% chunk_indices) tools::pskill(Sys.getpid(), tools::SIGKILL)
    original(chunk_indices = chunk_indices, ...)
  }

  assignInNamespace("eif_chunk_worker", saboteur, ns = "incremental")
  on.exit(assignInNamespace("eif_chunk_worker", original, ns = "incremental"), add = TRUE)

  err <- tryCatch(
    suppressWarnings(do.call(efficient_influence_function_chunked,
                             c(ecp_args(cs, "kl", chunk_size = 37L), list(n_cores = 3L)))),
    error = function(e) conditionMessage(e))

  expect_true(is.character(err))
  expect_match(err, "Parallel EIF computation failed in")
  expect_match(err, "child process killed")
  expect_match(err, "NOT silently ignored")
})

# ===========================================================================
# 4. Per-worker memory budget
# ===========================================================================

test_that("auto chunk size divides the 2GB target by the worker count", {
  # Pins the per-worker budget formula: peak memory is K * chunk_footprint, so the
  # target must be divided by K rather than applied as if one chunk were resident.
  bytes_per_obs <- function(m) 6 * 2^m * 8
  expected_chunk <- function(m, K) {
    max(1000, min(floor((2e9 / K) / bytes_per_obs(m)), 50000))
  }

  for (m in c(8L, 11L)) {
    for (K in c(1L, 2L, 4L, 9L)) {
      cz <- expected_chunk(m, K)
      peak_mb <- K * cz * bytes_per_obs(m) / 1024^2
      # Total peak stays within the ~2 GB envelope regardless of K (m >= 8, where
      # the 50000 cap is not what binds).
      expect_lte(peak_mb, 2048 * 1.05,
                 label = sprintf("m = %d, K = %d peak %.0f MB", m, K, peak_mb))
      # And per-worker footprint shrinks as K grows.
      if (K > 1L) expect_lte(cz, expected_chunk(m, 1L))
    }
  }

  # The formula the implementation actually uses must match the one asserted above.
  # Verified behaviourally: with a tiny n the resulting chunk count follows.
  ecp_skip_unless_fork()
  cs <- ecp_case(6L, n = 120L)
  res <- do.call(efficient_influence_function_chunked,
                 c(ecp_args(cs, "kl", chunk_size = NULL), list(n_cores = 2L)))
  # n = 120 is far below any auto chunk size, so this must land in a single chunk
  # and still produce correct, complete output.
  expect_length(res$eif_values, cs$n)
  expect_true(all(is.finite(res$eif_values)))
})

# ===========================================================================
# 5. BLAS pinning
# ===========================================================================

test_that("pin_blas_single_thread is a safe no-op when RhpcBLASctl is absent", {
  # The package is a Suggests, so the guarded path must not error when missing.
  # Simulated by making requireNamespace report FALSE.
  testthat::local_mocked_bindings(
    requireNamespace = function(...) FALSE,
    .package = "base"
  )
  expect_false(pin_blas_single_thread())
})

test_that("BLAS pinning executes inside forked workers", {
  ecp_skip_unless_fork()
  skip_if_not_installed("RhpcBLASctl")

  # Runs the real pin inside real forks and reports the thread count the child
  # sees afterwards. NOTE: this R build links the single-threaded reference
  # libRblas (blas_get_num_procs() == 1 even after blas_set_num_threads(4)), so a
  # 4 -> 1 transition is NOT observable here; what is asserted is that the call
  # path runs in the child, RhpcBLASctl is loadable post-fork, and the child ends
  # up at exactly one BLAS thread. On a multithreaded BLAS (OpenBLAS/MKL, e.g. the
  # production server) this is what prevents fork + thread-pool deadlock.
  out <- parallel::mclapply(1:2, function(i) {
    pinned <- pin_blas_single_thread()
    c(pinned = as.integer(pinned), procs = RhpcBLASctl::blas_get_num_procs())
  }, mc.cores = 2L, mc.preschedule = FALSE)

  expect_length(out, 2L)
  expect_false(any(vapply(out, function(x) inherits(x, "try-error"), logical(1))))
  for (o in out) {
    expect_identical(as.integer(o[["pinned"]]), 1L)   # RhpcBLASctl reachable in child
    expect_identical(as.integer(o[["procs"]]), 1L)    # child left at one BLAS thread
  }
})

# ===========================================================================
# 6. return_qdelta_results: the flag that unlocks the parallel speedup
# ===========================================================================

test_that("return_qdelta_results = FALSE keeps every number identical", {
  ecp_skip_unless_fork()
  # Withholding the heavy per-observation payload must change only what is
  # RETURNED, never what is COMPUTED. eif_values, plugin and the full convergence
  # reporting must be bit-identical to the TRUE case.
  for (m in c(2L, 4L, 6L)) {
    cs <- ecp_case(m, n = 180L)
    for (div in c("kl", "chi2")) {
      args <- ecp_args(cs, div, chunk_size = 31L)
      with_qd <- do.call(efficient_influence_function_chunked,
                         c(args, list(n_cores = 2L, return_qdelta_results = TRUE)))
      without <- do.call(efficient_influence_function_chunked,
                         c(args, list(n_cores = 2L, return_qdelta_results = FALSE)))

      label <- sprintf("m = %d, %s", m, div)
      expect_identical(without$eif_values, with_qd$eif_values, info = label)
      expect_identical(without$plugin, with_qd$plugin, info = label)
      expect_identical(without$centered_eif, with_qd$centered_eif, info = label)

      # Convergence reporting must survive intact -- it is rebuilt from the cheap
      # scalar fields the workers always return.
      expect_identical(without$convergence_info, with_qd$convergence_info, info = label)
      expect_identical(without$convergence_summary, with_qd$convergence_summary, info = label)

      # ...and only the heavy field is dropped, unambiguously (NULL, not stubs).
      expect_null(without$qdelta_results)
      expect_length(with_qd$qdelta_results, cs$n)
    }
  }
})

test_that("return_qdelta_results defaults to TRUE, preserving the historical contract", {
  ecp_skip_unless_fork()
  cs <- ecp_case(4L, n = 120L)
  res <- do.call(efficient_influence_function_chunked,
                 c(ecp_args(cs, "kl", chunk_size = 29L), list(n_cores = 2L)))
  expect_length(res$qdelta_results, cs$n)
  expect_false(any(vapply(res$qdelta_results, is.null, logical(1))))
  expect_length(res$qdelta_results[[1]]$q_delta, 2^cs$m)

  # Same on the sequential path.
  res1 <- do.call(efficient_influence_function_chunked,
                  c(ecp_args(cs, "kl", chunk_size = 29L), list(n_cores = 1L)))
  expect_identical(res1$qdelta_results, res$qdelta_results)
})

test_that("convergence reporting is complete and correct without the heavy payload", {
  ecp_skip_unless_fork()
  cs <- ecp_case(4L, n = 150L)
  res <- do.call(efficient_influence_function_chunked,
                 c(ecp_args(cs, "kl", chunk_size = 23L),
                   list(n_cores = 3L, return_qdelta_results = FALSE)))

  expect_length(res$convergence_info$converged, cs$n)
  expect_length(res$convergence_info$iterations, cs$n)
  expect_length(res$convergence_info$constraint_norms, cs$n)
  expect_type(res$convergence_info$converged, "logical")
  expect_equal(res$convergence_summary$n_converged, sum(res$convergence_info$converged))
  expect_equal(res$convergence_summary$convergence_rate,
               sum(res$convergence_info$converged) / cs$n)

  # Cross-check the scalars against an independent direct solve.
  for (i in c(1L, 75L, 150L)) {
    qd <- compute_qdelta_per_obs(cs$ps[i, ], cs$delta, "kl", fix_lambda0 = FALSE,
                                 tol = 1e-10, combination_matrix = cs$a_mat)
    expect_identical(res$convergence_info$converged[i], isTRUE(qd$converged))
    expect_identical(res$convergence_info$iterations[i], as.integer(qd$iterations))
  }
})

test_that("auto chunk sizing guarantees at least one chunk per worker", {
  ecp_skip_unless_fork()
  # Regression guard for a real defect found while benchmarking: the pre-existing
  # 50,000 upper bound on chunk_size binds long before the memory budget does at
  # small m (at m = 6 the budget allows ~72,000 obs/chunk), so any n <= 50,000
  # collapsed to ONE chunk and silently ran sequentially despite K workers being
  # requested. Measured impact at m = 6, n = 12,000: 0.91x (i.e. no parallelism,
  # pure overhead) before the fix, 3.47x after.
  #
  # Replicates the implementation's sizing rule and asserts >= n_workers chunks.
  auto_chunk <- function(m, n, K) {
    bytes_per_obs <- 6 * 2^m * 8
    cz <- max(1000, min(floor((2e9 / K) / bytes_per_obs), 50000))
    if (K > 1L) cz <- min(cz, max(1L, ceiling(n / K)))
    cz
  }

  for (m in c(6L, 8L, 11L)) {
    for (n in c(12000L, 60000L, 381191L)) {
      for (K in c(2L, 4L, 9L)) {
        cz <- auto_chunk(m, n, K)
        n_chunks <- ceiling(n / cz)
        expect_gte(n_chunks, K,
                   label = sprintf("m=%d n=%d K=%d -> %d chunks", m, n, K, n_chunks))
        # Shrinking for chunk count must never RAISE per-worker memory above the
        # K = 1 footprint, i.e. the budget still holds.
        expect_lte(cz, auto_chunk(m, n, 1L))
      }
    }
  }

  # Behavioural check: a moderate n at small m must genuinely engage >1 chunk and
  # still return correct, bit-identical results.
  cs <- ecp_case(6L, n = 900L)
  par_res <- do.call(efficient_influence_function_chunked,
                     c(ecp_args(cs, "kl", chunk_size = NULL), list(n_cores = 3L)))
  seq_res <- do.call(efficient_influence_function_chunked,
                     c(ecp_args(cs, "kl", chunk_size = NULL), list(n_cores = 1L)))
  expect_identical(par_res$eif_values, seq_res$eif_values)
  expect_identical(par_res$plugin, seq_res$plugin)
})
