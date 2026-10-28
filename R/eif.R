#' @file eif.R
#' @title Efficient Influence Function Computation
#' @description
#' Implements efficient influence function (EIF) for vector incremental treatment
#' effects. Combines outcome regression and inverse probability weighting
#' components into the one-step estimator analyzed in the accompanying paper
#' (nuisance error enters only through a second-order product-rate remainder;
#' this is not double robustness -- see Theorem 1 of the paper). Includes
#' per-observation Q_delta solver and EIF lambda derivatives.
#' @author Denis Agniel

# Constants used from constants.R:
# - PROPENSITY_SUM_TOL: tolerance for propensity scores sum-to-1 check (1e-6)
# - PROPENSITY_MIN_THRESHOLD: clipping threshold for near-zero propensities (1e-8)
# - MIN_CONVERGENCE_RATE: warn if Q_delta convergence < 0.95
# - NEAR_ZERO_THRESHOLD: general "effectively zero" threshold (1e-10)
# - OMEGA_REGULARIZATION_EPS: regularization for near-singular Omega (1e-4)
#
# h(x)/h'(x) (divergence-specific, boundary-clipped) are computed via the
# canonical divergence_h_pair() in divergence_h_functions.R, not inline here.

#' Default EIF chunk size from a per-worker memory budget (internal)
#'
#' `eif_compute_chunk()` holds about `EIF_LIVE_MATRICES` n_chunk x 2^m double matrices at its peak
#' (propensity and outcome slices, indicator, h, h', Q, Y, Phi, the residual and its temporaries),
#' counted from that function's body. The earlier model assumed 6 and a floor of 1000 observations,
#' which at m >= 12 allocated several times the stated budget. Chunking only partitions the
#' observations, so the chunk size never changes the estimate.
#'
#' @param m Integer, number of treatments.
#' @param n_workers Integer, number of parallel workers (>= 1); K chunks are resident at once.
#' @return Integer chunk size in `[EIF_MIN_CHUNK, EIF_MAX_CHUNK]`.
#' @keywords internal
eif_default_chunk_size <- function(m, n_workers = 1L) {
  bytes_per_obs <- EIF_LIVE_MATRICES * 2^m * 8
  per_worker_bytes <- EIF_CHUNK_BUDGET_BYTES / max(1L, n_workers)
  chunk_size <- floor(per_worker_bytes / bytes_per_obs)
  as.integer(max(EIF_MIN_CHUNK, min(chunk_size, EIF_MAX_CHUNK)))
}

#' Compute EIF for a chunk of observations (internal helper)
#'
#' @description
#' Internal function that computes EIF components for a subset of observations.
#' Used by efficient_influence_function_chunked() for memory-efficient computation.
#'
#' @param chunk_indices Integer vector of observation indices for this chunk
#' @param A Full treatment matrix (n x m)
#' @param Y Full outcome vector (length n)
#' @param R Full missingness indicator (length n, or NULL)
#' @param pR Full missingness probabilities (n x 2^m, or NULL)
#' @param propensity_scores Full propensity score matrix (n x 2^m)
#' @param outcome_means Full outcome means matrix (n x 2^m)
#' @param qdelta_results Full list of Q_delta results (length n)
#' @param a_mat Treatment combination matrix (2^m x m)
#' @param a_extended_matrix Extended treatment matrix with intercept (2^m x (m+1))
#' @param divergence Character, divergence type
#' @param m Integer, number of treatments
#'
#' @return List with eif_values, plugin, and convergence_info for the chunk
#' @keywords internal
eif_compute_chunk <- function(chunk_indices, A, Y, R, pR, propensity_scores, outcome_means,
                              qdelta_results, a_mat, a_extended_matrix, divergence, m) {

  n_chunk <- length(chunk_indices)
  n_combinations <- 2^m

  # Extract chunk data
  A_chunk <- A[chunk_indices, , drop = FALSE]
  Y_chunk <- Y[chunk_indices]
  propensity_scores_chunk <- propensity_scores[chunk_indices, , drop = FALSE]
  outcome_means_chunk <- outcome_means[chunk_indices, , drop = FALSE]
  qdelta_results_chunk <- qdelta_results[chunk_indices]

  if (!is.null(R)) {
    R_chunk <- R[chunk_indices]
    pR_chunk <- pR[chunk_indices, , drop = FALSE]
  } else {
    R_chunk <- NULL
    pR_chunk <- NULL
  }

  # Compute indicator matrix for chunk using optimized bit-pattern matching
  A_patterns <- as.vector(A_chunk %*% (2^(0:(m-1))))
  a_patterns <- as.vector(a_mat %*% (2^(0:(m-1))))
  I_chunk <- outer(A_patterns, a_patterns, "==") * 1.0

  # Extract lambda and compute H matrices for chunk
  LAMBDA_chunk <- as.matrix(do.call(rbind, lapply(qdelta_results_chunk, function(x) x$lambda)))
  H_result <- compute_H_vectorized_general(LAMBDA_chunk, a_extended_matrix, divergence)
  H_mat_chunk <- H_result$H_mat
  H_prime_chunk <- H_result$H_prime

  # Compute EIF_lambda for each observation in chunk
  EIF_LAMBDA_chunk <- list()
  for (i in 1:n_chunk) {
    EIF_LAMBDA_chunk[[i]] <- calculate_eif_lambda(
      A_chunk[i, ], propensity_scores_chunk[i, ],
      qdelta_results_chunk[[i]],
      H_mat_chunk[i, ], H_prime_chunk[i, ],
      a_extended_matrix, divergence
    )
  }
  EIF_LAMBDA_chunk <- as.matrix(do.call(rbind, EIF_LAMBDA_chunk))

  # Extract Q_delta values for chunk
  Q_chunk <- as.matrix(do.call(rbind, lapply(qdelta_results_chunk, function(x) x$q_delta)))

  # h_chunk is the same h(u) already computed above as H_mat_chunk via
  # compute_H_vectorized_general() (single canonical source of h/h', see
  # divergence_h_functions.R). Previously this was recomputed a second time
  # here with its own (differently-clipped) formula -- redundant work, and
  # a source of eif1-vs-eif2/eif3 inconsistency near domain boundaries for
  # reverse_kl/hellinger. See [LEARN 2026-08-10] in divergence_h_functions.R.
  h_chunk <- H_mat_chunk

  # Expand Y_chunk to matrix for broadcasting (n_chunk x n_combinations)
  Y_mat <- matrix(Y_chunk, nrow = n_chunk, ncol = n_combinations)

  # Compute Phi component using h instead of dividing by pi
  # Phi = h * (I * (Y - mu) + pi * mu)
  if (is.null(R_chunk) && is.null(pR_chunk)) {
    # No missingness case
    Phi_chunk <- h_chunk * (I_chunk * (Y_mat - outcome_means_chunk) +
                            propensity_scores_chunk * outcome_means_chunk)
  } else {
    # With missingness
    # Phi = h * (R * I * (Y - mu) / pR + pi * mu)
    safe_pR <- pmax(pR_chunk, 1e-100)
    Phi_chunk <- h_chunk * (R_chunk * I_chunk * (Y_mat - outcome_means_chunk) / safe_pR +
                            propensity_scores_chunk * outcome_means_chunk)
  }

  # Compute residual component
  HResid_chunk <- (I_chunk - propensity_scores_chunk) * H_mat_chunk

  # Compute EIF components
  eif1_chunk <- rowSums(Phi_chunk)
  eif2_chunk <- rowSums(outcome_means_chunk * HResid_chunk)
  eif3_chunk <- rowSums(outcome_means_chunk * propensity_scores_chunk * H_prime_chunk *
                        (EIF_LAMBDA_chunk %*% t(a_extended_matrix)))

  plugin_chunk <- rowSums(outcome_means_chunk * Q_chunk)
  eif_values_chunk <- eif1_chunk + eif2_chunk - eif3_chunk

  return(list(
    eif_values = eif_values_chunk,
    plugin = plugin_chunk
  ))
}


#' Resolve the worker count for chunk-level parallelism
#'
#' @description
#' Central policy for how many parallel workers the chunked EIF pass uses, so the
#' rule lives in one place and is testable.
#'
#' Resolution order: explicit `n_cores` argument, then
#' `getOption("incremental.n_cores")` (a global escape hatch for scripts that must
#' force sequential execution without editing every call site), then the spec's
#' default of `parallel::detectCores() - 1`.
#'
#' **Forking is Unix-only.** `parallel::mclapply()` does not fork on Windows — it
#' silently degrades to `lapply()`. Rather than let that happen implicitly (and
#' then size chunks for workers that will never exist), this function returns 1 on
#' Windows, so the sequential path and the single-chunk memory budget are chosen
#' deliberately. That is a platform limitation, not a defect.
#'
#' @param n_cores Integer or NULL. NULL means "auto" (see above). Values < 1 are
#'   treated as 1.
#' @return Integer >= 1, the number of workers to use.
#' @keywords internal
resolve_eif_workers <- function(n_cores = NULL) {
  if (is.null(n_cores)) n_cores <- getOption("incremental.n_cores", default = NULL)

  if (is.null(n_cores)) {
    detected <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    if (is.na(detected)) detected <- 1L
    n_cores <- max(1L, detected - 1L)
  }

  n_cores <- suppressWarnings(as.integer(n_cores))
  if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

  # mclapply cannot fork on Windows; do not pretend otherwise.
  if (.Platform$OS.type != "unix") n_cores <- 1L

  n_cores
}

#' Pin BLAS to a single thread inside a forked worker
#'
#' @description
#' The stage-1/stage-2 vectorization made BLAS (`crossprod`, `%*%`) the hot path
#' inside every chunk. Combining `fork()` with a multithreaded BLAS
#' (OpenBLAS/MKL/Accelerate) is a classic hang: the BLAS thread pool does not
#' survive the fork, and a child that tries to use it can deadlock. Pinning each
#' child to one BLAS thread avoids that, and also prevents `K` workers x `T` BLAS
#' threads from oversubscribing the machine.
#'
#' `RhpcBLASctl` is the only reliable *runtime* control (env vars such as
#' `OMP_NUM_THREADS` are read when the BLAS loads, i.e. at R startup, so setting
#' them in a child is too late). It is therefore a `Suggests`, used when present.
#' When absent this is a no-op: on a single-threaded reference BLAS -- what this
#' package's own R build uses, `blas_get_num_procs()` == 1 -- there is nothing to
#' pin and no hazard.
#'
#' @return Invisibly, `TRUE` if the thread count was actually set, `FALSE` if
#'   `RhpcBLASctl` was unavailable.
#' @keywords internal
pin_blas_single_thread <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

#' Compute one observation chunk end to end (parallel unit of work)
#'
#' @description
#' The unit dispatched to each parallel worker. It owns BOTH phases for its own
#' observations:
#'   1. the per-observation Q_delta Newton solves (~82% of per-observation cost at
#'      m = 11 after stage 1), when `qdelta_results` was not supplied by the caller;
#'   2. the vectorized EIF matrix assembly (~18%), via `eif_compute_chunk()`.
#'
#' Both phases must live inside the worker for chunk parallelism to be worth
#' anything: before this stage the Q_delta loop ran serially in the parent, which
#' by Amdahl's law capped any chunk-level speedup at ~1.19x regardless of core
#' count (measured). Nothing about either phase's math is changed here -- this
#' function only slices, dispatches and re-assembles.
#'
#' All matrices are sliced to the chunk INSIDE the worker. Under `fork()` the
#' parent's full matrices are shared copy-on-write and only read, so each worker's
#' resident footprint is its own chunk, not the whole cohort.
#'
#' @param chunk_indices Integer vector of observation indices for this chunk
#' @param A,Y,R,pR,propensity_scores,outcome_means Full-cohort inputs (sliced here)
#' @param qdelta_results Full-length list of pre-computed Q_delta results, or NULL
#'   to solve this chunk's observations inside the worker
#' @param delta Numeric vector (length m), required when `qdelta_results` is NULL
#' @param a_mat,a_extended_matrix Combination matrices, built once by the caller
#' @param divergence Character, divergence type
#' @param m Integer, number of treatments
#' @param qdelta_tol Numeric, solver tolerance (the caller's 1e-10 contract)
#' @param pin_blas Logical, pin BLAS to one thread (TRUE inside forked workers)
#' @param return_qdelta Logical, ship the full per-observation Q_delta results
#'   (each containing a `2^m` weight vector) back to the parent. Setting this FALSE
#'   returns only the three scalar convergence fields per observation, which is
#'   dramatically cheaper across a fork pipe -- see the "Return payload dominates
#'   fork speedup" section of `efficient_influence_function_chunked()`.
#' @return List with `eif_values`, `plugin`, `chunk_indices`, `qdelta_conv` (three
#'   scalar vectors, present whenever this worker solved Q_delta), and
#'   `qdelta_results` (the full list; present only when solved here AND
#'   `return_qdelta` is TRUE, so already-known or unwanted results are never
#'   serialized back from the worker)
#' @keywords internal
eif_chunk_worker <- function(chunk_indices, A, Y, R, pR, propensity_scores,
                             outcome_means, qdelta_results, delta,
                             a_mat, a_extended_matrix, divergence, m,
                             qdelta_tol = 1e-10, pin_blas = FALSE,
                             return_qdelta = TRUE) {

  if (pin_blas) pin_blas_single_thread()

  n_chunk <- length(chunk_indices)

  # Slice to this chunk. Passing chunk_indices = seq_len(n_chunk) to
  # eif_compute_chunk() below makes its own internal slicing an identity op, so
  # that (stage-2 verified) function is used completely unmodified.
  ps_chunk <- propensity_scores[chunk_indices, , drop = FALSE]
  A_chunk <- A[chunk_indices, , drop = FALSE]
  Y_chunk <- Y[chunk_indices]
  mu_chunk <- outcome_means[chunk_indices, , drop = FALSE]
  R_chunk <- if (is.null(R)) NULL else R[chunk_indices]
  pR_chunk <- if (is.null(pR)) NULL else pR[chunk_indices, , drop = FALSE]

  # Phase 1: Q_delta. Solve here if the caller did not supply results.
  solved_here <- is.null(qdelta_results)
  if (solved_here) {
    qd_chunk <- vector("list", n_chunk)
    for (j in seq_len(n_chunk)) {
      qd_chunk[[j]] <- compute_qdelta_per_obs(
        ps_chunk[j, ], delta, divergence,
        fix_lambda0 = FALSE, tol = qdelta_tol,
        combination_matrix = a_mat
      )
    }
  } else {
    qd_chunk <- qdelta_results[chunk_indices]
  }

  # Phase 2: vectorized EIF assembly (unmodified stage-2 code path).
  chunk_result <- eif_compute_chunk(
    chunk_indices = seq_len(n_chunk),
    A = A_chunk, Y = Y_chunk, R = R_chunk, pR = pR_chunk,
    propensity_scores = ps_chunk,
    outcome_means = mu_chunk,
    qdelta_results = qd_chunk,
    a_mat = a_mat,
    a_extended_matrix = a_extended_matrix,
    divergence = divergence,
    m = m
  )

  # Scalar convergence summary for this chunk: three numbers per observation, so
  # this is always cheap enough to return even when the full q_delta vectors are
  # deliberately withheld.
  qdelta_conv <- if (solved_here) {
    list(
      converged = vapply(qd_chunk, function(x) isTRUE(x$converged), logical(1)),
      iterations = vapply(qd_chunk, function(x) {
        as.integer(if (is.null(x$iterations)) 0L else x$iterations)
      }, integer(1)),
      constraint_norms = vapply(qd_chunk, function(x) {
        as.numeric(if (is.null(x$constraint_norm)) 0 else x$constraint_norm)
      }, numeric(1))
    )
  } else {
    NULL
  }

  list(
    eif_values = chunk_result$eif_values,
    plugin = chunk_result$plugin,
    chunk_indices = chunk_indices,
    qdelta_conv = qdelta_conv,
    qdelta_results = if (solved_here && return_qdelta) qd_chunk else NULL
  )
}

#' Chunked EIF computation for large datasets (internal)
#'
#' @description
#' Memory-efficient version of EIF computation that processes observations in chunks.
#' Automatically used when n > 10,000 to avoid memory issues with large matrices.
#'
#' Chunks are independent, so they are optionally computed in parallel across
#' observation chunks (see `n_cores`). Parallelism is deliberately over the
#' OBSERVATION axis only -- never over delta scenarios or coalitions, which are
#' left free for future chaining/attribution work.
#'
#' @section Memory budget is per worker, not global:
#' The auto chunk-size rule (`eif_default_chunk_size()`) targets ~2 GB of working matrices
#' (`EIF_LIVE_MATRICES * 2^m * 8` bytes per observation, 14 live matrices counted from
#' `eif_compute_chunk()`). With `K` workers, `K` chunks are resident simultaneously, so the target
#' is divided by `K`:
#' \preformatted{
#'   bytes_per_obs      <- 14 * 2^m * 8
#'   target_chunk_bytes <- 2e9 / n_workers        # <- per-worker budget
#'   chunk_size         <- floor(target_chunk_bytes / bytes_per_obs)
#'   chunk_size         <- max(25, min(chunk_size, 50000))
#' }
#' so total peak stays near the budget. Dividing by `K` also guarantees roughly `K` times more
#' chunks, which keeps `K` workers busy. The previous model assumed 6 matrices and a floor of
#' 1000 observations, which exceeded the budget at m >= 12; the per-observation `q_delta` vectors
#' held by the parent are not counted here.
#'
#' @section Return payload dominates fork speedup:
#' Measured, not assumed (m = 11, K = 9, `benchmark_eif_chunk_parallel.R` and
#' `log_bench_chunk_parallel.txt`): with the full per-observation `qdelta_results`
#' shipped back from workers the speedup is **0.73x** -- i.e. *slower* than
#' sequential -- while with the same work and only the scalar convergence fields
#' returned it is **4.70x**. The cause is `mclapply()`'s return path: each
#' observation's `q_delta` is a `2^m` double vector (16 KB at m = 11), so a chunk
#' of 500 observations serializes ~8 MB through a pipe, and every worker's payload
#' must be deserialized by the single parent process. That parent-side
#' deserialization is serial, so it becomes the bottleneck and gets *worse* with
#' larger chunks (speedup fell 3.05x -> 1.09x -> 0.66x as chunk size rose
#' 200 -> 500 -> 1000).
#'
#' Hence `return_qdelta_results`. It defaults to TRUE so the documented return
#' contract is unchanged, but no code in `R/` or `application/` actually reads the
#' heavy `q_delta`/`lambda` vectors off an EIF result (only the scalar `converged`
#' flag, which `convergence_info` provides either way). Callers that do not need
#' them should pass FALSE and get the full parallel speedup.
#'
#' @inheritParams efficient_influence_function
#' @param chunk_size Integer, observations per chunk (computed automatically if NULL)
#' @param n_cores Integer, number of parallel workers over observation chunks.
#'   NULL (default) resolves via `getOption("incremental.n_cores")` then
#'   `parallel::detectCores() - 1`. 1 forces the sequential path. Always 1 on
#'   Windows, where `mclapply()` cannot fork. See `resolve_eif_workers()`.
#' @param return_qdelta_results Logical, include the full per-observation Q_delta
#'   results in the output. TRUE (default) preserves the historical return
#'   contract. FALSE returns `qdelta_results = NULL` while keeping
#'   `convergence_info` and `convergence_summary` fully populated, and is what
#'   unlocks the parallel speedup -- see "Return payload dominates fork speedup".
#'
#' @return Same as efficient_influence_function
#' @keywords internal
efficient_influence_function_chunked <- function(X = NULL, A, Y = NULL, R = NULL, pR = NULL,
                                                propensity_scores, delta, outcome_means = NULL,
                                                divergence = "kl", qdelta_results = NULL,
                                                feasibility_check = FALSE,
                                                chunk_size = NULL,
                                                n_cores = NULL,
                                                return_qdelta_results = TRUE) {

  # Determine dimensions
  if (!is.null(X)) {
    n <- nrow(X)
  } else if (!is.null(A)) {
    n <- nrow(A)
  } else {
    n <- nrow(propensity_scores)
  }

  m <- ncol(A)
  n_combinations <- 2^m

  n_workers <- resolve_eif_workers(n_cores)

  # Determine optimal chunk size if not provided
  if (is.null(chunk_size)) {
    # Target: ~2 GB of working matrices in total, i.e. ~2 GB / K per worker, since
    # K chunks are resident at once. See eif_default_chunk_size() for the per-observation model.
    chunk_size <- eif_default_chunk_size(m, n_workers)

    # Ensure there is actually work for every worker. At small m the 50,000 upper
    # bound binds long before the memory budget does, which would otherwise yield a
    # single chunk and silently degrade to sequential execution: at m = 6 the budget
    # allows ~72,000 observations per chunk, so any n <= 50,000 became ONE chunk and
    # got no parallelism at all. Shrinking chunk_size to give >= n_workers chunks
    # only ever REDUCES per-worker memory, so the budget above still holds.
    if (n_workers > 1L) {
      chunk_size <- min(chunk_size, max(1L, ceiling(n / n_workers)))
    }
  }

  # Compute number of chunks
  n_chunks <- ceiling(n / chunk_size)

  # No point forking more workers than there are chunks.
  n_workers <- min(n_workers, n_chunks)

  # Prepare treatment combinations (compute once)
  a_mat <- generate_binary_vectors(m)
  a_extended_matrix <- cbind(1, a_mat)

  # Convert inputs to matrices
  propensity_scores <- as.matrix(propensity_scores)
  outcome_means <- as.matrix(outcome_means)
  A <- as.matrix(A)
  Y <- as.numeric(Y)

  # Initialize result vectors
  eif_values <- numeric(n)
  plugin <- numeric(n)

  # Whether the workers must solve Q_delta themselves (they do unless the caller
  # supplied results). Recorded before dispatch because the workers only ship
  # qdelta_results back when they computed them.
  solve_qdelta_in_workers <- is.null(qdelta_results)
  if (solve_qdelta_in_workers) qdelta_results <- vector("list", n)

  chunk_index_list <- lapply(seq_len(n_chunks), function(k) {
    start_idx <- (k - 1L) * chunk_size + 1L
    end_idx <- min(k * chunk_size, n)
    start_idx:end_idx
  })

  # Optional progress reporting for very large datasets
  if (n > 50000 && n_chunks > 1) {
    cat(sprintf("Computing EIF for %d observations in %d chunks (chunk size: %d, workers: %d)...\n",
                n, n_chunks, chunk_size, n_workers))
  }

  worker_args <- list(
    A = A, Y = Y, R = R, pR = pR,
    propensity_scores = propensity_scores,
    outcome_means = outcome_means,
    qdelta_results = if (solve_qdelta_in_workers) NULL else qdelta_results,
    delta = delta,
    a_mat = a_mat,
    a_extended_matrix = a_extended_matrix,
    divergence = divergence,
    m = m,
    return_qdelta = return_qdelta_results
  )

  if (n_workers > 1L) {
    # Fork-based parallelism over observation chunks. mc.preschedule = FALSE gives
    # one fork per chunk, which keeps each child's peak footprint to a single
    # chunk (the basis of the per-worker memory budget above) and lets a child's
    # memory be reclaimed on exit rather than accumulating across chunks.
    chunk_results <- parallel::mclapply(
      chunk_index_list,
      function(ci) do.call(eif_chunk_worker,
                           c(list(chunk_indices = ci, pin_blas = TRUE), worker_args)),
      mc.cores = n_workers,
      mc.preschedule = FALSE
    )

    # Re-raise worker failures LOUDLY. mclapply returns a "try-error" element for
    # a child that errored, and NULL for a child killed outright (e.g. OOM); both
    # would otherwise flow onward and surface as silent NA/garbage in eif_values.
    failed <- vapply(chunk_results, function(x) inherits(x, "try-error"), logical(1))
    missing_res <- vapply(chunk_results, is.null, logical(1))

    if (length(chunk_results) != n_chunks || any(failed) || any(missing_res)) {
      bad <- which(failed | missing_res)
      if (length(chunk_results) != n_chunks) {
        stop("Parallel EIF computation returned ", length(chunk_results),
             " chunk results but ", n_chunks, " chunks were dispatched. ",
             "At least one worker died without reporting. ",
             "Retry with n_cores = 1 to obtain the underlying error.",
             call. = FALSE)
      }
      first_msg <- if (any(failed)) {
        conditionMessage(attr(chunk_results[[which(failed)[1]]], "condition"))
      } else {
        "worker returned NULL (child process killed, e.g. out of memory)"
      }
      stop("Parallel EIF computation failed in ", length(bad), " of ", n_chunks,
           " observation chunk(s) (chunk index ", paste(utils::head(bad, 5), collapse = ", "),
           if (length(bad) > 5) ", ..." else "", "). First error: ", first_msg,
           "\nThis was NOT silently ignored. Retry with n_cores = 1 for a ",
           "sequential traceback.",
           call. = FALSE)
    }
  } else {
    # Sequential path: identical work, no fork. Also the Windows path.
    chunk_results <- vector("list", n_chunks)
    for (k in seq_len(n_chunks)) {
      chunk_results[[k]] <- do.call(
        eif_chunk_worker,
        c(list(chunk_indices = chunk_index_list[[k]], pin_blas = FALSE), worker_args)
      )
      if (n > 50000 && k %% 5 == 0) {
        cat(sprintf("  Processed %d / %d observations (%.1f%%)\n",
                    max(chunk_index_list[[k]]), n,
                    100 * max(chunk_index_list[[k]]) / n))
      }
      gc(verbose = FALSE)
    }
  }

  # Re-assemble in observation order (each result carries its own indices, so the
  # scatter is order-independent and parallel execution cannot permute anything).
  worker_conv <- if (solve_qdelta_in_workers) {
    list(converged = logical(n), iterations = integer(n), constraint_norms = numeric(n))
  } else {
    NULL
  }
  covered <- logical(n)

  for (res in chunk_results) {
    ci <- res$chunk_indices
    eif_values[ci] <- res$eif_values
    plugin[ci] <- res$plugin
    covered[ci] <- TRUE
    if (solve_qdelta_in_workers) {
      worker_conv$converged[ci] <- res$qdelta_conv$converged
      worker_conv$iterations[ci] <- res$qdelta_conv$iterations
      worker_conv$constraint_norms[ci] <- res$qdelta_conv$constraint_norms
      if (return_qdelta_results) qdelta_results[ci] <- res$qdelta_results
    }
  }
  rm(chunk_results)
  gc(verbose = FALSE)

  # Every observation must have been written by exactly one chunk. Cheap guard
  # against a chunking defect silently leaving zeros in eif_values.
  if (!all(covered)) {
    stop("Chunk assembly covered only ", sum(covered), " of ", n,
         " observations. This indicates a chunking defect.", call. = FALSE)
  }
  if (solve_qdelta_in_workers && !return_qdelta_results) qdelta_results <- NULL

  # Compute convergence summary. When the workers solved Q_delta they returned the
  # three scalar fields directly (cheap across a fork pipe), so use those; that
  # keeps convergence reporting complete even when the heavy per-observation
  # results were deliberately not shipped back.
  convergence_info <- if (!is.null(worker_conv)) {
    worker_conv
  } else {
    list(
      converged = sapply(qdelta_results, function(x) x$converged),
      iterations = sapply(qdelta_results, function(x) ifelse(is.null(x$iterations), 0, x$iterations)),
      constraint_norms = sapply(qdelta_results, function(x) ifelse(is.null(x$constraint_norm), 0, x$constraint_norm))
    )
  }

  convergence_summary <- list(
    n_converged = sum(convergence_info$converged),
    convergence_rate = sum(convergence_info$converged) / n,
    mean_iterations = mean(convergence_info$iterations[convergence_info$converged]),
    max_iterations = max(convergence_info$iterations),
    mean_constraint_norm = mean(convergence_info$constraint_norms[convergence_info$converged]),
    max_constraint_norm = max(convergence_info$constraint_norms)
  )

  # Validate EIF computation results
  if (any(!is.finite(eif_values))) {
    n_nonfinite <- sum(!is.finite(eif_values))
    stop("EIF computation produced ", n_nonfinite, " non-finite values (NA, NaN, or Inf). ",
         "This indicates numerical instability. ",
         "Possible causes: (1) propensity scores too close to 0, (2) outcome means unstable, ",
         "(3) Q_delta solver issues. ",
         "Try: increasing regularization or checking for positivity violations.",
         call. = FALSE)
  }

  # Check if centered EIF has reasonable properties
  centered_eif <- eif_values - mean(eif_values)
  eif_sd <- sd(centered_eif)
  if (!is.na(eif_sd) && eif_sd < NEAR_ZERO_THRESHOLD) {
    warning("EIF has extremely small variance (sd = ", sprintf("%.2e", eif_sd), "). ",
            "This may indicate degenerate estimates. Check data and nuisance functions.")
  }

  return(list(
    eif_values = eif_values,
    centered_eif = centered_eif,
    plugin = plugin,
    qdelta_results = qdelta_results,
    convergence_info = convergence_info,
    convergence_summary = convergence_summary
  ))
}


#' Compute Efficient Influence Function for Vector Incremental Effects
#'
#' @description
#' Computes the efficient influence function (EIF) for vector incremental treatment
#' effects, combining outcome-regression and inverse-probability-weighting
#' components into the one-step estimator analyzed in the paper. The EIF provides
#' the basis for asymptotically efficient estimation and valid inference via the
#' central limit theorem: sqrt(n) * (psi_hat - psi) -> N(0, E[EIF^2]).
#'
#' For large datasets (n > 10,000), automatically uses memory-efficient chunked
#' computation to avoid memory issues.
#'
#' @param X Numeric matrix (n x p) of covariates. Optional but recommended for
#'   diagnostics. If NULL, will be extracted from other inputs.
#' @param A Numeric matrix (n x m) of binary treatments (values 0 or 1).
#'   Required for computing IPW components.
#' @param Y Numeric vector (length n) of observed outcomes. Required unless
#'   R and pR are provided for missingness adjustment.
#' @param R Numeric vector (length n) of missingness indicators (1 = observed,
#'   0 = missing). Optional. If provided, pR must also be provided.
#' @param pR Numeric matrix (n x 2^m) of missingness probabilities P(R=1|A,X).
#'   Optional. If provided, R must also be provided. Used for missing outcome adjustment.
#' @param propensity_scores Numeric matrix (n x 2^m) of propensity scores
#'   P(A=a|X) for all 2^m treatment combinations. Rows must sum to 1.
#' @param delta Numeric vector (length m) of target marginal changes. Positive
#'   values required. Represents odds ratios for marginal treatment probabilities.
#' @param outcome_means Numeric matrix (n x 2^m) of conditional outcome means
#'   E[Y|A=a,X] for all treatment combinations. Required for EIF-based estimation.
#' @param divergence Character, divergence measure for Q_delta optimization.
#'   Options: "kl" (default, most stable), "reverse_kl", "chi2", "hellinger".
#' @param qdelta_results List (length n) of pre-computed Q_delta solver results.
#'   Each element from compute_qdelta_per_obs(). If NULL, computed automatically.
#'   Pre-computing can improve performance when calling multiple times.
#' @param feasibility_check Logical, perform delta feasibility check before
#'   optimization. If TRUE and delta is infeasible, suggests alternatives.
#'   Default: FALSE (faster).
#' @param chunk_size Integer, observations per chunk for the large-n (n > 10,000)
#'   chunked path. NULL (default) auto-sizes to a ~2 GB total working-matrix
#'   budget, divided across parallel workers. Ignored when n <= 10,000.
#' @param n_cores Integer, number of parallel workers for the large-n chunked path,
#'   parallelizing over OBSERVATION chunks (never over delta scenarios or
#'   coalitions). NULL (default) resolves via `getOption("incremental.n_cores")`
#'   and then `parallel::detectCores() - 1`; pass 1 to force sequential execution.
#'   Always sequential on Windows, where `parallel::mclapply()` cannot fork.
#'   Ignored when n <= 10,000. Parallelism re-orders independent chunk work only
#'   and does not change any numeric result.
#' @param return_qdelta_results Logical, include the full per-observation Q_delta
#'   results (`qdelta_results`) in the output. TRUE (default) preserves the
#'   historical return contract. On the large-n chunked path, returning these
#'   `2^m`-length vectors from every forked worker is what limits the parallel
#'   speedup (measured 0.73x with them vs 4.70x without, at m = 11 / K = 9);
#'   pass FALSE if the caller does not need them. `convergence_info` and
#'   `convergence_summary` stay fully populated either way. Ignored when
#'   n <= 10,000.
#'
#' @details
#' **EIF Structure:**
#' The EIF has three components (see paper for full derivation):
#' \enumerate{
#'   \item \strong{Outcome regression:} E_Q[Y - mu(A,X)] - corrects for outcome model misspecification
#'   \item \strong{Propensity score:} sum_a [I(A=a) - p(a|X)] * h(a, lambda) * mu(a,X) - corrects for propensity misspecification
#'   \item \strong{Lambda derivative:} -sum_a mu(a,X) * p(a|X) * h'(a, lambda) * EIF_lambda(a) - accounts for Q_delta estimation
#' }
#'
#' **Double robustness:** The estimator is consistent if EITHER the outcome models
#' OR the propensity score models are correctly specified (but not necessarily both).
#'
#' **Computation details:**
#' - If qdelta_results not provided, uses Newton-Raphson to solve m+1 constraints per observation
#' - Convergence diagnostics returned: converged flags, iteration counts, constraint norms
#' - Warnings issued if convergence rate < 0.95 or if EIF variance is degenerate
#' - Input validation checks: propensity sum-to-1, non-finite values, feasibility
#' - For large n (> 10,000), automatically uses chunked computation to reduce memory usage
#'
#' **Performance tips:**
#' - Pre-compute qdelta_results when calling repeatedly with same delta
#' - Use feasibility_check = FALSE unless delta feasibility is uncertain
#' - For very large datasets, chunked computation is automatic and transparent
#'
#' @return List containing:
#' \describe{
#'   \item{eif_values}{Numeric vector (length n) of EIF values}
#'   \item{centered_eif}{Numeric vector (length n) of centered EIF values}
#'   \item{plugin}{Numeric vector (length n) of plugin estimates}
#'   \item{qdelta_results}{List of per-observation Q_delta solver results}
#'   \item{convergence_info}{List with converged, iterations, constraint_norms vectors}
#'   \item{convergence_summary}{List with n_converged, convergence_rate, mean_iterations, etc.}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data
#' sim_data <- simulate_vector_treatment_data(n = 500, m = 2, p = 2, seed = 123)
#'
#' # Estimate nuisance functions
#' ps_est <- estimate_propensity_scores(sim_data$data, method = "glm")
#' outcome_est <- estimate_outcome_models(sim_data$data, method = "glm")
#'
#' # Compute EIF for delta = (0.1, 0.1)
#' eif_result <- efficient_influence_function(
#'   A = as.matrix(sim_data$data[, c("A1", "A2")]),
#'   Y = sim_data$data$Y,
#'   propensity_scores = ps_est,
#'   outcome_means = outcome_est,
#'   delta = c(0.1, 0.1),
#'   divergence = "kl"
#' )
#'
#' # Check convergence
#' print(eif_result$convergence_summary)
#' }
efficient_influence_function <- function(X = NULL, A, Y = NULL, R = NULL, pR = NULL,
                                         propensity_scores,
                                         delta,
                                         outcome_means = NULL,
                                         divergence = "kl",
                                         qdelta_results = NULL,
                                         feasibility_check = FALSE,
                                         chunk_size = NULL,
                                         n_cores = NULL,
                                         return_qdelta_results = TRUE) {

  # Determine dimensions
  if (!is.null(X)) {
  n <- nrow(X)
  } else if (!is.null(A)) {
    n <- nrow(A)
  } else {
    n <- nrow(propensity_scores)
  }

  m <- ncol(A)
  n_combinations <- 2^m

  # Validate inputs
  if (sum(c(is.null(R), is.null(pR))) == 1) stop("Either both R and pR must be null or non-null")
  if (is.null(outcome_means)) stop("Outcome means must be provided")
  if (is.null(Y)) stop("Y must be provided")

  # Validate propensity scores
  if (any(!is.finite(propensity_scores))) {
    stop("Propensity scores contain non-finite values (NA, NaN, or Inf). ",
         "Check data quality and propensity score estimation.")
  }

  # Check propensity scores sum to 1 (within tolerance)
  ps_row_sums <- rowSums(propensity_scores)
  if (!all(abs(ps_row_sums - 1) < PROPENSITY_SUM_TOL)) {
    max_deviation <- max(abs(ps_row_sums - 1))
    stop("Propensity scores do not sum to 1 (max deviation: ", sprintf("%.2e", max_deviation), "). ",
         "Normalize propensity scores before passing to EIF computation.")
  }

  # Validate outcome means
  if (any(!is.finite(outcome_means))) {
    stop("Outcome means contain non-finite values (NA, NaN, or Inf). ",
         "Check data quality and outcome model estimation.")
  }

  # Validate Y
  if (any(!is.finite(Y))) {
    stop("Outcome Y contains non-finite values (NA, NaN, or Inf). ",
         "Handle missing outcomes before EIF computation or provide R and pR for missingness adjustment.")
  }

  # Validate delta
  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a positive numeric vector. ",
         "Received: ", paste(delta, collapse = ", "))
  }

  if (feasibility_check) {
    # Early feasibility check across all observations (with early stopping for fast failure)
    feasibility_analysis <- analyze_delta_feasibility_all_obs(propensity_scores, delta,
                                                              early_stop = TRUE)

    if (feasibility_analysis$n_feasible < n) {
      # Some observations are infeasible

      # Check if early stopped
      if (!is.null(feasibility_analysis$early_stopped) && feasibility_analysis$early_stopped) {
        # Early stop: only checked first few observations
        stop(sprintf(
          "Delta infeasible for observation %d (first infeasible found).\nOnly checked %d / %d observations before stopping early.\nAll %d preceding observations were feasible.\nConsider using a smaller delta (all components < %.2f).",
          feasibility_analysis$first_infeasible,
          feasibility_analysis$first_infeasible,
          n,
          feasibility_analysis$n_feasible,
          min(delta) * 0.9
        ), call. = FALSE)
      }

      # Full check: try to find feasible alternative
      suggested_result <- suggest_feasible_delta(propensity_scores, delta)

      if (suggested_result$n_obs_feasible == n) {
        stop(paste("Delta infeasible for", n - feasibility_analysis$n_feasible, "out of", n, "observations.",
                   "\nSuggested feasible delta:", paste(round(suggested_result$suggested_delta, 4), collapse = ", "),
                   "\nUse this delta or modify your requested delta."))
      } else {
        stop(paste("Delta infeasible for", n - feasibility_analysis$n_feasible, "out of", n, "observations.",
                   "\nNo universally feasible delta found. Consider using different delta values or",
                   "\nchecking the support patterns in your propensity scores."))
      }
    }
  }

  # Combination matrix depends only on m, so it is built AT MOST ONCE per call
  # and shared by the Q_delta solver loop below and the EIF matrix pass further
  # down. NULL here means "not built yet".
  a_mat <- NULL

  # Validate any caller-supplied qdelta_results up front. This must happen before
  # the chunked dispatch below, because both paths rely on it.
  if (!is.null(qdelta_results)) {
    if (length(qdelta_results) != n) {
      stop("qdelta_results must have length equal to number of observations (", n, ")")
    }
    for (i in 1:n) {
      if (is.null(qdelta_results[[i]])) {
        stop("qdelta_results[[", i, "]] is NULL")
      }
    }
  }

  # Dispatch to chunked computation for large datasets (memory efficient, and the
  # path that parallelizes over observation chunks).
  # Threshold: 10,000 observations (conservative to avoid memory issues)
  #
  # NOTE (2026-08-11): the Q_delta solve is deliberately NOT performed before this
  # dispatch. It used to be -- a serial `for (i in 1:n)` loop right here -- which
  # left the chunked pass covering only the EIF assembly. Measured at m = 11 that
  # is 17.6% of per-observation cost, so parallelizing chunks would have been
  # capped by Amdahl's law at ~1.19x no matter how many cores were available.
  # Passing qdelta_results through (possibly NULL) lets each chunk worker own BOTH
  # phases for its own observations, putting ~100% of the per-observation work
  # inside the parallel region. For the small-n path below, the serial solve is
  # unchanged.
  if (n > 10000) {
    return(efficient_influence_function_chunked(
      X = X, A = A, Y = Y, R = R, pR = pR,
      propensity_scores = propensity_scores,
      delta = delta,
      outcome_means = outcome_means,
      divergence = divergence,
      qdelta_results = qdelta_results,
      feasibility_check = FALSE,  # Already checked above if needed
      chunk_size = chunk_size,
      n_cores = n_cores,
      return_qdelta_results = return_qdelta_results
    ))
  }

  # Use provided qdelta_results or compute them (small-n, non-chunked path)
  if (is.null(qdelta_results)) {
    # Build once and thread it into every per-observation solve, instead of each
    # solve regenerating it (and regenerating it again inside
    # handle_zero_propensities()).
    a_mat <- generate_binary_vectors(m)

    # Compute Q_delta for each observation
    qdelta_results <- list()
    convergence_info <- list(
      converged = logical(n),
      iterations = integer(n),
      constraint_norms = numeric(n)
    )

    for (i in 1:n) {
      # Compute Q_delta for this observation using internal normalization with tighter tolerance
      qdelta_results[[i]] <- compute_qdelta_per_obs(propensity_scores[i, ], delta, divergence,
                                              fix_lambda0 = FALSE, tol = 1e-10,
                                              combination_matrix = a_mat)

      # Store convergence info
      convergence_info$converged[i] <- qdelta_results[[i]]$converged
      convergence_info$iterations[i] <- qdelta_results[[i]]$iterations
      convergence_info$constraint_norms[i] <- qdelta_results[[i]]$constraint_norm
    }
  } else {
    # Extract convergence info from provided results (already validated above)
    convergence_info <- list(
      converged = logical(n),
      iterations = integer(n),
      constraint_norms = numeric(n)
    )
    for (i in 1:n) {
      convergence_info$converged[i] <- qdelta_results[[i]]$converged
      convergence_info$iterations[i] <- ifelse(is.null(qdelta_results[[i]]$iterations), 0, qdelta_results[[i]]$iterations)
      convergence_info$constraint_norms[i] <- ifelse(is.null(qdelta_results[[i]]$constraint_norm), 0, qdelta_results[[i]]$constraint_norm)
    }
  }

  # Compute EIF values using first form (vectorized for small n)
  eif_values <- numeric(n)

  propensity_scores <- as.matrix(propensity_scores)  # n × nA
  outcome_means <- as.matrix(outcome_means)          # n × nA (if provided as matrix)
  A <- as.matrix(A)                                  # n × m
  # Reuse the matrix built for the Q_delta loop above; only the
  # qdelta_results-supplied branch reaches here without one.
  if (is.null(a_mat)) {
    a_mat <- generate_binary_vectors(m)   # nA × q (rows are patterns)
  }
  a_extended_matrix <- cbind(1, a_mat)                # nA × (m+1) - compute once!
  Y <- as.numeric(Y)

  # Check for missing Y values - do not silently impute
  if (any(is.na(Y))) {
    if (is.null(R) || is.null(pR)) {
      stop("Y contains missing values but R and pR are not provided. ",
           "Either provide R (missingness indicator) and pR (missingness probabilities), ",
           "or remove observations with missing Y values.",
           call. = FALSE)
    }
  }

  # Dimension validation
  expected_combinations <- 2^m
  if (ncol(propensity_scores) != expected_combinations) {
    stop(sprintf(
      "propensity_scores has %d columns but should have %d (2^%d) for m=%d treatments",
      ncol(propensity_scores), expected_combinations, m, m
    ), call. = FALSE)
  }
  if (ncol(outcome_means) != expected_combinations) {
    stop(sprintf(
      "outcome_means has %d columns but should have %d (2^%d) for m=%d treatments",
      ncol(outcome_means), expected_combinations, m, m
    ), call. = FALSE)
  }
  if (nrow(propensity_scores) != n) {
    stop(sprintf(
      "propensity_scores has %d rows but n=%d",
      nrow(propensity_scores), n
    ), call. = FALSE)
  }
  if (nrow(outcome_means) != n) {
    stop(sprintf(
      "outcome_means has %d rows but n=%d",
      nrow(outcome_means), n
    ), call. = FALSE)
  }

  q   <- ncol(A)
  nA  <- nrow(a_mat)

  # Optimized indicator matrix using bit-pattern matching (faster than string-based outer)
  # Convert rows to integer patterns for fast matching
  A_patterns <- as.vector(A %*% (2^(0:(m-1))))
  a_patterns <- as.vector(a_mat %*% (2^(0:(m-1))))
  I <- outer(A_patterns, a_patterns, "==") * 1.0

  LAMBDA <- as.matrix(do.call(rbind, lapply(qdelta_results, function(x) x$lambda)))  # n × (m+1)
  H_result <- compute_H_vectorized_general(LAMBDA, a_extended_matrix, divergence)
  H_mat <- H_result$H_mat
  H_prime <- H_result$H_prime

  EIF_LAMBDA <- list()
  for (i in 1:n) {
    EIF_LAMBDA[[i]] <- calculate_eif_lambda(A[i,], propensity_scores[i,],
                                            qdelta_results[[i]],
                                            H_mat[i,], H_prime[i,],
                                            a_extended_matrix, divergence)
  }
  EIF_LAMBDA <- as.matrix(do.call(rbind, EIF_LAMBDA))

  # --- Core matrices ---------------------------------------------------------
  # Track whether missingness is present before creating defaults
  has_missingness <- !is.null(R) && !is.null(pR)

  if (is.null(R)) {
    R  <- rep(1, length(Y))
    pR <- matrix(rep(1, nrow(propensity_scores)*ncol(propensity_scores)),
                              nrow(propensity_scores),
                              ncol(propensity_scores))
  }

  # Extract all Q_delta values upfront
  Q <- as.matrix(do.call(rbind, lapply(qdelta_results, function(x) x$q_delta)))

  # Extract lambda values
  LAMBDA <- as.matrix(do.call(rbind, lapply(qdelta_results, function(x) x$lambda)))

  # h is the same h(u) already computed above as H_mat via
  # compute_H_vectorized_general() (single canonical source of h/h', see
  # divergence_h_functions.R). Previously this was recomputed a second time
  # here with its own (differently-clipped) formula -- redundant work, and
  # a source of eif1-vs-eif2/eif3 inconsistency near domain boundaries for
  # reverse_kl/hellinger. See [LEARN 2026-08-10] in divergence_h_functions.R.
  h <- H_mat

  # Expand Y to matrix for broadcasting (n x n_combinations)
  n_combinations <- ncol(outcome_means)
  Y_mat <- matrix(Y, nrow = n, ncol = n_combinations)

  # Compute Phi component using h instead of dividing by pi
  # Phi = h * (I * (Y - mu) + pi * mu)
  if (!has_missingness) {
    # No missingness
    Phi <- h * (I * (Y_mat - outcome_means) + propensity_scores * outcome_means)
  } else {
    # With missingness
    safe_pR <- pmax(pR, 1e-100)
    Phi <- h * (R * I * (Y_mat - outcome_means) / safe_pR + propensity_scores * outcome_means)
  }

  HResid <- (I - propensity_scores) * H_mat

  eif1 <- rowSums(Phi)
  eif2 <- rowSums(outcome_means * HResid)
  eif3 <- rowSums(outcome_means * propensity_scores * H_prime * (EIF_LAMBDA %*% t(a_extended_matrix)))

  plugin <- rowSums(outcome_means * Q)
  eif_values <- eif1 + eif2 - eif3

  # Compute summary statistics
  convergence_summary <- list(
    n_converged = sum(convergence_info$converged),
    convergence_rate = sum(convergence_info$converged) / n,
    mean_iterations = mean(convergence_info$iterations[convergence_info$converged]),
    max_iterations = max(convergence_info$iterations),
    mean_constraint_norm = mean(convergence_info$constraint_norms[convergence_info$converged]),
    max_constraint_norm = max(convergence_info$constraint_norms)
  )

  # Warn about low convergence rate
  if (convergence_summary$convergence_rate < MIN_CONVERGENCE_RATE) {
    n_failed <- n - convergence_summary$n_converged
    warning(sprintf(
      "Q_delta solver failed to converge for %d out of %d observations (%.1f%% convergence). ",
      n_failed, n, convergence_summary$convergence_rate * 100,
      "This may indicate: (1) infeasible delta values, (2) numerical instability, ",
      "(3) insufficient max_iter. Consider checking delta feasibility or increasing tolerance."
    ), call. = FALSE)
  }

  # Validate EIF computation results
  if (any(!is.finite(eif_values))) {
    n_nonfinite <- sum(!is.finite(eif_values))
    stop("EIF computation produced ", n_nonfinite, " non-finite values (NA, NaN, or Inf). ",
         "This indicates numerical instability. ",
         "Possible causes: (1) propensity scores too close to 0, (2) outcome means unstable, ",
         "(3) Q_delta solver issues. ",
         "Try: increasing regularization or checking for positivity violations.")
  }

  # Check if centered EIF has reasonable properties
  centered_eif <- eif_values - mean(eif_values)
  eif_sd <- sd(centered_eif)
  if (!is.na(eif_sd) && eif_sd < NEAR_ZERO_THRESHOLD) {
    warning("EIF has extremely small variance (sd = ", sprintf("%.2e", eif_sd), "). ",
            "This may indicate degenerate estimates. Check data and nuisance functions.")
  }

  return(list(
    eif_values = eif_values,
    centered_eif = centered_eif,
    plugin = plugin,
    qdelta_results = qdelta_results,
    convergence_info = convergence_info,
    convergence_summary = convergence_summary
  ))
}

compute_H_vectorized_general <- function(LAMBDA, a_extended_matrix, divergence) {
  # LAMBDA: n × (m+1) matrix of lambda values
  # a_extended_matrix: nA × (m+1) matrix of treatment patterns with intercept
  # divergence: character, one of "kl", "reverse_kl", "chi2", "hellinger"
  # Returns: list with H_mat and H_prime, both n × nA

  linear_terms <- LAMBDA %*% t(a_extended_matrix)  # n × nA

  # Canonical h/h' (see divergence_h_functions.R) -- this used to be
  # reimplemented here with its own boundary-clipping constants, which
  # disagreed with fast_qdelta_solver.R's convention despite the comments
  # above claiming to match it. See [LEARN 2026-08-10] in
  # divergence_h_functions.R.
  hp <- divergence_h_pair(linear_terms, divergence)

  return(list(H_mat = hp$h, H_prime = hp$h_prime))
}

#' Per-observation EIF contribution for the Lagrange multipliers lambda
#'
#' @description
#' Computes \eqn{\Omega^{-1} g(\bar A; \lambda)} for one observation, the
#' lambda-derivative term of the efficient influence function.
#'
#' @section Omega assembly (matrix form):
#' \eqn{\Omega = \sum_a \pi_a h'(u_a)\, \bar a\, \bar a^\top}, i.e. in matrix form
#' \eqn{\Omega = \bar A^\top \mathrm{diag}(w) \bar A} with
#' \eqn{w = \pi \odot h'(u)} and \eqn{\bar A} the extended
#' \eqn{2^m \times (m+1)} design. Before 2026-08-11 this was coded literally:
#' `diag(weights)` materialized a dense \eqn{2^m \times 2^m} matrix (33.6 MB at
#' m = 11) **per observation** and the triple product multiplied through it,
#' costing \eqn{O(2^{2m} m)} ~ 50M flops instead of the \eqn{O(2^m m^2)} ~ 0.3M
#' actually required. Profiling put those two lines at 95.8% of the whole EIF pass
#' at m = 11 (see `explorations/investigations/profile_eif_pass.md`).
#'
#' It is now assembled as `crossprod(a_mat * weights, a_mat)`: scale each ROW of
#' \eqn{\bar A} by its weight (`weights` has length \eqn{2^m} = `nrow(a_mat)`, so
#' `*` recycles down columns = row scaling), then form the
#' \eqn{(m+1) \times 2^m} by \eqn{2^m \times (m+1)} product. This is an exact
#' algebraic identity, not an approximation --
#' \eqn{\sum_k w_k \bar a_{ki} \bar a_{kj}} either way -- differing only by
#' floating-point summation order. The same identity is already used by
#' `compute_omega_reduced()` (`R/treatment_subset_utils.R`), which is the
#' active-set path's Omega builder.
#'
#' @param A Numeric vector of length m, observed treatment pattern for this
#'   observation (without intercept).
#' @param propensity_scores Numeric vector of length 2^m, \eqn{p_0(a \mid X_i)}.
#' @param qdelta_result List for this observation from `compute_qdelta_per_obs()`
#'   (needs `lambda`, `q_delta`, `delta`).
#' @param h_values Numeric vector of length 2^m, \eqn{h(u_a)}.
#' @param h_dot_values Numeric vector of length 2^m, \eqn{h'(u_a)}.
#' @param a_mat The **extended** \eqn{2^m \times (m+1)} design `cbind(1, combos)`.
#'   Note the parameter name: both call sites (`eif_compute_chunk()` and the
#'   non-chunked path in `efficient_influence_function()`) pass
#'   `a_extended_matrix` here, and the body relies on that (`a_mat[, -1]` for the
#'   marginals, `ncol(a_mat)` as m+1, rows compared against `c(1, A)`).
#' @param divergence Character, divergence type (unused directly; h/h' arrive
#'   pre-computed via `divergence_h_pair()`).
#' @return Numeric vector of length m+1, zero-padded over inactive constraints.
#' @keywords internal
calculate_eif_lambda <- function(A, propensity_scores, qdelta_result, h_values, h_dot_values, a_mat, divergence) {
  # Get treatment assignment for this observation
  #treatment_combinations <- generate_treatment_combinations(m = length(A))

  # Get Q_delta result for this observation
  lambda <- qdelta_result$lambda
  q_delta <- qdelta_result$q_delta
  delta <- qdelta_result$delta

  # Determine active constraints for this observation.
  #
  # Marginals come from the FULL extended design, dropping the intercept from the
  # RESULT rather than slicing it out of the INPUT: compute_marginals() returns
  # `joint_dist %*% vectors`, so element 1 is sum(propensity_scores) and elements
  # 2..(m+1) are exactly the marginals `a_mat[, -1]` would have produced (verified
  # bitwise identical for m in 1..11 over random propensity vectors). The old
  # `a_mat[, -1, drop = FALSE]` copied a 2^m x m matrix on every observation
  # (184 KB at m = 11) and was the dominant residual cost -- 49% of
  # calculate_eif_lambda() -- once the Omega and pattern-lookup fixes landed.
  # Slicing nothing also removes the reason the old `drop = FALSE` guard existed
  # (m = 1, where `a_mat[, -1]` collapsed to a vector and broke compute_marginals).
  E_A <- compute_marginals(propensity_scores, a_mat, m = ncol(a_mat))[-1]
  active_constraints <- find_active_constraints(E_A, tol = 1e-8)
  active_constraints <- c(1, (active_constraints + 1))

  # Calculate the Omega matrix (e.g. E[h(A)h(A)'] for KL)
  #a_extended_matrix <- cbind(1, treatment_combinations)
  n_combinations <- length(propensity_scores)
  weights <- as.numeric(propensity_scores * h_dot_values)

  # Omega = t(a_mat) %*% diag(weights) %*% a_mat, WITHOUT materializing the dense
  # 2^m x 2^m diagonal matrix or multiplying through it. `a_mat * weights` scales
  # row k by weights[k] (recycling down columns), so this is the same weighted
  # second-moment sum, exact up to floating-point summation order. See the "Omega
  # assembly" section above.
  omega_matrix <- crossprod(a_mat * weights, a_mat)
  omega_active <- omega_matrix[active_constraints, active_constraints]

  # Number of active constraints (including intercept)
  n_active <- length(active_constraints)

  # Compute inverse with regularization if needed; is there a faster function than solve?
  regularization_applied <- FALSE
  tryCatch({
    omega_inv <- solve(omega_active)
    # Check for Inf/NaN in result (can happen with near-singular matrices)
    if (any(is.infinite(omega_inv)) || any(is.na(omega_inv))) {
      stop("Matrix inverse contains Inf/NaN")
    }
  }, error = function(e) {
    # Add stronger regularization for near-singular matrices
    regularization_applied <<- TRUE
    omega_inv <<- solve(omega_active + OMEGA_REGULARIZATION_EPS * diag(n_active))
  })

  if (regularization_applied) {
    warning(
      "Omega matrix regularization applied in EIF lambda calculation. ",
      "This indicates the Omega matrix is near-singular, which may occur with: ",
      "(1) extreme propensity scores, (2) redundant constraints, (3) numerical instability. ",
      "Results may be less reliable.",
      call. = FALSE
    )
  }

  # Compute g(bar_A; lambda) from paper equation (199)
  # g = h(bar_A) * bar_A - E_Q(bar_A|X) - Delta^2 * (bar_A - E(bar_A|X))

  # bar_A = (1, A) - observed treatment extended
  a_prime_extended <- c(1, A)

  # Find h value for observed pattern bar_A.
  #
  # Two-step exact lookup, replacing `apply(a_mat, 1, function(row) all(row ==
  # a_prime_extended))` (2^m closure calls per observation; 2.6% of the pre-fix
  # EIF pass at m = 11 but ~70% of what remained after the Omega fix, and the
  # dominant residual line):
  #   1. PREFILTER on an integer bit key. Rows of the extended design are binary,
  #      so `row %*% 2^(0:m)` is an exact integer (sums of distinct powers of two,
  #      exactly representable in double for m+1 <= 53) and is injective on binary
  #      rows -- one BLAS dgemv instead of 2^m interpreted closure calls.
  #   2. VERIFY the surviving candidate(s) with the original elementwise test, so
  #      semantics are preserved exactly even if a caller passes a non-binary
  #      `A` (where the key could collide) or a design with duplicate rows: no
  #      match and multiple matches both still fall through to the stop() below,
  #      as before.
  # Measured at m = 11: 1653 us (apply) -> 24 us, identical return value.
  pattern_pow2 <- 2^(seq_len(ncol(a_mat)) - 1L)
  candidate_idx <- which(as.vector(a_mat %*% pattern_pow2) ==
                           sum(a_prime_extended * pattern_pow2))
  pattern_idx <- candidate_idx[vapply(
    candidate_idx,
    function(k) all(a_mat[k, ] == a_prime_extended),
    logical(1)
  )]
  if (length(pattern_idx) != 1) {
    stop("Could not uniquely identify observed treatment pattern in calculate_eif_lambda")
  }
  h_A <- h_values[pattern_idx]  # h(bar_A; lambda) - scalar

  # E_Q(bar_A | X) - extended marginals under Q_delta
  E_QA <- as.numeric(q_delta %*% a_mat)

  # E(bar_A | X) - extended marginals under P_0
  E_A <- c(1, E_A)

  # Diagonal of \bar{Delta} (definition: Delta_{kk} = delta_k / (1 + P_0(A_k=1|X)(delta_k - 1)))
  delta_extended <- c(1, delta)  # (1, δ_1, ..., δ_m)
  Delta_diag <- delta_extended / (1 + E_A * (delta_extended - 1))
  Delta2_diag <- Delta_diag^2

  # From implicit function theorem derivation:
  # ∂ζ_k/∂π_a = (Δ_{kk}^2 / δ_k) · a_k
  # So g should have term: Δ^2/δ * (bar_A - E(bar_A))
  Delta2_over_delta <- Delta2_diag / delta_extended  # Δ^2 / δ

  # g(bar_A; lambda) from my derivation
  g <- h_A * a_prime_extended - E_QA - Delta2_over_delta * (a_prime_extended - E_A)

  # EIF_lambda = Omega^{-1} * g (for active constraints)
  # Note: minus sign applied in EIF computation (line 338: eif_values = eif1 + eif2 - eif3)
  eif_active <- as.numeric(omega_inv %*% g[active_constraints])

  # Pad to full dimension (m+1) with zeros for inactive constraints
  # This ensures consistent dimensions when rbind-ing across observations
  m_plus_1 <- ncol(a_mat)
  if (length(active_constraints) < m_plus_1) {
    eif_full <- rep(0, m_plus_1)
    eif_full[active_constraints] <- eif_active
    return(eif_full)
  }
  eif_active
}

#' Compute Q_delta distribution for a single observation
#'
#' @description
#' Computes the observation-specific Q_delta distribution q_delta(a|X_i) by solving
#' a constrained optimization problem. Finds lambda parameters that tilt the
#' propensity distribution p_0(a|X_i) to achieve target marginal probabilities
#' E_Q[A_k] = E_P[A_k] + delta_k while minimizing a specified divergence.
#'
#' @param propensity_scores_i Numeric vector (length 2^m) of propensity scores
#'   P(A=a|X_i) for all treatment combinations. Must sum to 1. Combinations with
#'   zero probability are handled automatically.
#' @param delta Numeric vector (length m) of target marginal changes. Positive
#'   odds ratios. E.g., delta = c(1.5, 2.0) increases P(A1=1) by 1.5x and P(A2=1) by 2x.
#' @param divergence Character, divergence measure to minimize. Options:
#'   - "kl" (default): Kullback-Leibler D(Q||P), most stable, always convex
#'   - "reverse_kl": Reverse KL D(P||Q), can be unstable near boundary
#'   - "chi2": Chi-squared divergence, quadratic penalty
#'   - "hellinger": Hellinger distance, bounded divergence
#' @param fix_lambda0 Logical, normalization approach:
#'   - FALSE (default): Internal normalization, solves (m+1) x (m+1) system, more stable
#'   - TRUE: External normalization, fixes lambda_0, solves m x m system
#' @param tol Numeric, convergence tolerance for constraint norm. Default: 1e-8.
#'   Smaller values = tighter convergence but may require more iterations.
#' @param check_feasibility Logical, check if delta is feasible before optimization.
#'   If TRUE and infeasible, provides diagnostic error. Default: FALSE (faster).
#' @param adaptive_constraints Logical, use adaptive constraint dropping when
#'   Jacobian is singular (e.g., at polytope boundary). Default: FALSE. Set TRUE
#'   for robust handling of near-boundary cases.
#' @param max_constraint_drops Integer, maximum constraints to drop if
#'   adaptive_constraints = TRUE. Default: 2.
#' @param verbose_adaptive Logical, print diagnostics when dropping constraints.
#'   Default: FALSE.
#' @param combination_matrix Optional `2^m x m` matrix of all binary treatment
#'   combinations, as returned by [generate_binary_vectors()]. Depends only on
#'   `m`, so callers looping over many observations at a fixed `m` should compute
#'   it once and pass it in; `NULL` (default) computes it internally.
#'
#' @details
#' **Optimization problem:**
#' Solves: min_Q D(Q||P) subject to E_Q[A_k] = E_P[A_k] + delta_k for k=1,...,m
#'
#' The solution has form: q(a) = p_0(a) * h(a; lambda) where h depends on divergence:
#' - KL: h(a) = exp(lambda_0 + sum_k lambda_k * a_k)
#' - Reverse KL: h(a) = -1 / (lambda_0 + sum_k lambda_k * a_k)
#' - Chi2: h(a) = (lambda_0 + sum_k lambda_k * a_k) / 2 + 1
#' - Hellinger: h(a) = 1 / (1 - lambda_0 - sum_k lambda_k * a_k)^2
#'
#' **Algorithm:**
#' 1. Initialize lambda (divergence-specific starting values)
#' 2. Newton-Raphson: lambda <- lambda - J^{-1} * g where g is constraint vector
#' 3. Uses analytical Jacobians for speed (falls back to numerical if needed)
#' 4. Convergence when ||g|| < tol or max iterations reached
#'
#' **Convergence diagnostics:**
#' - converged: TRUE if ||constraints|| < tol
#' - iterations: Number of Newton steps taken
#' - constraint_norm: Final ||E_Q[A] - target|| norm
#'
#' **Special cases handled:**
#' - Zero propensities: Removes zero-probability combinations, works on active support
#' - Singular Jacobians: Uses regularization or adaptive constraint dropping
#' - Domain violations: Safeguards prevent h(a) from leaving valid domain
#'
#' @return List containing:
#' \describe{
#'   \item{q_delta}{Numeric vector (length 2^m), Q_delta weights summing to 1}
#'   \item{lambda}{Numeric vector (length m+1 or m), fitted Lagrange multipliers}
#'   \item{converged}{Logical, TRUE if converged within tolerance}
#'   \item{constraint_norm}{Numeric, final ||constraint|| value}
#'   \item{iterations}{Integer, number of Newton-Raphson iterations}
#'   \item{delta}{Numeric vector, input delta (for reference)}
#' }
#'
#' @seealso
#' \code{\link{efficient_influence_function}} for EIF computation using Q_delta,
#' \code{\link{estimate_qdelta_mean}} for end-to-end estimation
#'
#' @return List containing:
#' \describe{
#'   \item{q_delta}{Numeric vector (length 2^m) of Q_delta probabilities for each treatment combination}
#'   \item{lambda}{Numeric vector (length m+1 or m) of observation-specific lambda parameters}
#'   \item{converged}{Logical, whether optimization converged within tolerance}
#'   \item{constraint_norm}{Numeric, final constraint violation norm (should be < tol)}
#'   \item{iterations}{Integer, number of Newton-Raphson iterations}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # For a single observation with 2 treatments
#' propensity_scores_i <- c(0.25, 0.25, 0.25, 0.25)  # 4 combinations
#' result <- compute_qdelta_per_obs(
#'   propensity_scores_i = propensity_scores_i,
#'   delta = c(0.1, 0.1),
#'   divergence = "kl"
#' )
#' print(result$q_delta)  # Q_delta distribution
#' print(result$converged)  # Check convergence
#' }
compute_qdelta_per_obs <- function(propensity_scores_i, delta, divergence = "kl",
                                  fix_lambda0 = FALSE, tol = 1e-8, check_feasibility = FALSE,
                                  adaptive_constraints = FALSE, max_constraint_drops = 2,
                                  verbose_adaptive = FALSE, combination_matrix = NULL) {

  # Check feasibility before attempting optimization (optional)
  if (check_feasibility) {
    m <- length(delta)
    vectors <- if (is.null(combination_matrix)) generate_binary_vectors(m) else combination_matrix
    p0_marginals <- compute_marginals(propensity_scores_i, vectors)
    support_matrix <- propensity_scores_i > tol

    feasibility_result <- check_delta_feasibility(delta, p0_marginals, support_matrix, tol)

    if (!feasibility_result$feasible) {
      stop(paste("Delta infeasible for this observation. Reason:", feasibility_result$reason,
                 "\nTarget marginals:", paste(round(feasibility_result$target_marginals, 4), collapse = ", "),
                 "\nUse suggest_feasible_delta() to find feasible alternatives."))
    }
  }

  # Use fast Newton-Raphson solver (newton_qdelta)
    result <- newton_qdelta(propensity_scores_i, delta, divergence, fix_lambda0 = fix_lambda0,
                           tol = tol, adaptive_constraints = adaptive_constraints,
                           max_constraint_drops = max_constraint_drops,
                           verbose_adaptive = verbose_adaptive,
                           combination_matrix = combination_matrix)

    return(list(
      q_delta = result$q_delta_weights,
      lambda = result$lambda,
      converged = result$converged,
      constraint_norm = result$constraint_norm,
      iterations = result$iterations,
      delta = delta
    ))
}



#' Compute scalar derivative ḣ(a,λ) = dh/dx where x is the linear argument
#'
#' Computes the scalar derivative of h(a,λ) with respect to its linear argument
#' x = λ₀ + Σ λₖaₖ. This is needed for the EIF correction term.
#'
#' @param a Numeric vector, treatment combination
#' @param lambda Numeric vector, Lagrange multipliers
#' @param divergence Character, divergence type
#' @return Numeric scalar, ḣ(a,λ) = dh/dx
