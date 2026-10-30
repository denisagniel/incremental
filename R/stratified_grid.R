#' Pattern-major stratified estimation for several outcomes and delta scenarios at once
#'
#' @description
#' Estimates the eligibility-stratified counterfactual mean for every requested (outcome, delta
#' scenario) run, going through the eligibility patterns ONE AT A TIME. For each pattern it fits the
#' propensity model once and one outcome model per outcome, and for each scenario solves the
#' per-observation Q_delta problems once and reuses them for every outcome (they do not depend on the
#' outcome). Everything fitted for a pattern is dropped before the next pattern starts, so peak memory
#' is set by the largest single pattern, not by the sum over patterns. The per-pattern results are
#' written to disk as they finish (one file per pattern), so a crash loses at most one pattern, and
#' adding a scenario or an outcome later computes only what is missing.
#'
#' @details
#' The numbers equal those of \code{estimate_qdelta_mean_stratified()} run separately for each
#' (outcome, scenario) with a shared \code{nuisance_cache}: the per-fit seeds depend only on
#' \code{seed_base}, the outcome's position in \code{outcome_cols} and the pattern's rank
#' (\code{pattern_pos}), none of which change here. A solver failure or a failed pattern fails only the
#' runs it touches; the other runs are still returned. A cache mismatch (a programming error) aborts.
#'
#' @param data Data frame with the A1..Am, eligibility, covariate and outcome columns.
#' @param outcome_cols Named character vector: names are the outcome labels (their ORDER fixes the
#'   outcome-fit seeds, keep it stable), values the column names in \code{data}.
#' @param deltas Named list of delta vectors (length m), the scenarios.
#' @param runs Optional data frame with columns \code{outcome} and \code{scenario}: which runs to
#'   compute. Default: all combinations.
#' @param seed_base Integer base of the per-fit seeds (see \code{new_nuisance_cache()}).
#' @param checkpoint_dir Directory for the per-pattern result files, or NULL for none.
#' @param checkpoint_context List hashed into every unit signature (for example the data and package
#'   fingerprints of the calling script).
#' @inheritParams estimate_qdelta_mean_stratified
#' @return A list: \code{runs} (named \code{"<outcome>__<scenario>"}; each element is a stratified
#'   result list like \code{estimate_qdelta_mean_stratified()}'s, minus \code{pattern_results}, or
#'   \code{list(success = FALSE, error_message = ...)}), \code{eligibility_patterns},
#'   \code{patterns_to_use}, \code{feasibility_check}, \code{timing}.
#' @export
estimate_qdelta_stratified_grid <- function(data, outcome_cols, deltas, runs = NULL,
                                            covariate_cols, eligibility_cols,
                                            min_stratum_size = 100, min_realized_coverage = 0.8,
                                            max_patterns = NULL,
                                            ps_method = "glm", binary_method = "glm",
                                            outcome_method = "glm", estimator = "one_step",
                                            cross_fit = TRUE, n_folds = 5, divergence = "kl",
                                            xgb_params = NULL, max_qdelta_failure_rate = 0,
                                            qdelta_tol = NULL, skip_null_shift = FALSE,
                                            seed_base, checkpoint_dir = NULL,
                                            checkpoint_context = list(),
                                            confidence_level = 0.95, verbose = TRUE) {
  if (is.null(names(outcome_cols)) || anyDuplicated(names(outcome_cols))) {
    stop("`outcome_cols` must be a uniquely named character vector.", call. = FALSE)
  }
  if (is.null(names(deltas)) || anyDuplicated(names(deltas))) {
    stop("`deltas` must be a uniquely named list.", call. = FALSE)
  }
  need <- c(covariate_cols, eligibility_cols, unname(outcome_cols))
  if (length(miss <- setdiff(need, names(data)))) {
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (is.null(runs)) {
    runs <- expand.grid(outcome = names(outcome_cols), scenario = names(deltas), stringsAsFactors = FALSE)
  }
  stopifnot(all(runs$outcome %in% names(outcome_cols)), all(runs$scenario %in% names(deltas)))
  run_key <- paste(runs$outcome, runs$scenario, sep = "__")

  pattern_info <- identify_eligibility_patterns(data, eligibility_cols)
  m <- pattern_info$m
  for (d in deltas) if (length(d) != m) stop("every delta must have length ", m, call. = FALSE)
  feasibility_check <- validate_stratification_feasibility(pattern_info$patterns,
                                                           min_stratum_size = min_stratum_size,
                                                           min_coverage = 0.5)
  if (!feasibility_check$feasible) stop("Stratification not feasible: ", feasibility_check$reason, call. = FALSE)
  patterns_to_use <- pattern_info$patterns[pattern_info$patterns$pattern_id %in%
                                             feasibility_check$feasible_pattern_ids, ]
  if (!is.null(max_patterns) && nrow(patterns_to_use) > max_patterns) {
    patterns_to_use <- patterns_to_use[1:max_patterns, ]
  }
  n_pat <- nrow(patterns_to_use)
  if (verbose) message(sprintf("Grid: %d patterns x %d runs (%d outcomes, %d scenarios)",
                               n_pat, nrow(runs), length(outcome_cols), length(deltas)))

  a_cols <- grep("^A[0-9]+$", names(data), value = TRUE)
  if (!is.null(checkpoint_dir)) fs::dir_create(checkpoint_dir)
  pkg_version <- as.character(utils::packageVersion("incremental"))
  unit_results <- stats::setNames(vector("list", n_pat), seq_len(n_pat))   # per pattern: named list of results
  timing <- data.frame(pattern = seq_len(n_pat), n_obs = NA_integer_, n_eligible = NA_integer_,
                       seconds = NA_real_, computed = NA)

  for (i in seq_len(n_pat)) {
    t0 <- Sys.time()
    row <- patterns_to_use[i, ]
    pat <- list(pattern_id = row$pattern_id, pattern_string = row$pattern_string,
                eligibility = row$eligibility[[1]], n_obs = row$n_obs, prevalence = row$prevalence,
                n_eligible = row$n_eligible, obs_indices = row$obs_indices[[1]])
    shared_md5 <- rlang::hash(data[pat$obs_indices, c(covariate_cols, a_cols, eligibility_cols), drop = FALSE])
    y_md5 <- vapply(outcome_cols, function(cn) rlang::hash(data[[cn]][pat$obs_indices]), "")
    signature <- function(oc, sc) {
      list(version = 3L, pattern = pat$pattern_string, n_obs = pat$n_obs, n_full = nrow(data),
           pattern_pos = i, outcome = oc, outcome_index = match(oc, names(outcome_cols)),
           outcome_col = unname(outcome_cols[oc]), delta = as.numeric(deltas[[sc]]),
           divergence = divergence, estimator = estimator, ps_method = ps_method,
           binary_method = binary_method, outcome_method = outcome_method, cross_fit = cross_fit,
           n_folds = n_folds, xgb_params = xgb_params, covariate_cols = covariate_cols,
           qdelta_tol = resolve_qdelta_tol(qdelta_tol),
           max_qdelta_failure_rate = max_qdelta_failure_rate, skip_null_shift = skip_null_shift,
           seed_base = as.integer(seed_base), package_version = pkg_version,
           shared_md5 = shared_md5, y_md5 = unname(y_md5[oc]), context = checkpoint_context)
    }
    file <- if (is.null(checkpoint_dir)) NULL else
      file.path(checkpoint_dir, sprintf("pattern_%04d.rds", i))
    saved <- list()
    if (!is.null(file) && file.exists(file)) {
      ck <- tryCatch(readRDS(file), error = function(e) NULL)
      if (is.list(ck) && is.list(ck$units)) saved <- ck$units
    }
    have <- stats::setNames(vector("list", length(run_key)), run_key)
    todo <- logical(length(run_key))
    for (r in seq_along(run_key)) {
      u <- saved[[run_key[r]]]
      if (is.list(u) && identical(u$signature, signature(runs$outcome[r], runs$scenario[r])) &&
          isTRUE(u$result$success)) {
        have[[r]] <- u$result
      } else {
        todo[r] <- TRUE
      }
    }

    if (any(todo)) {
      cache <- new_nuisance_cache(seed_base, tags = names(outcome_cols))   # fresh: freed with this pattern
      for (sc in unique(runs$scenario[todo])) {
        memo <- new.env()                                                   # Q_delta solves shared by outcomes
        for (r in which(todo & runs$scenario == sc)) {
          oc <- runs$outcome[r]
          # downstream nuisance code refers to the outcome as "Y": present this outcome under that name
          # (a shallow copy of the data frame; only the Y vector is new)
          dd <- data
          dd[["Y"]] <- data[[unname(outcome_cols[oc])]]
          res <- tryCatch(
            estimate_pattern_effect(
              data = dd, pattern_info = pat, delta = deltas[[sc]], covariate_cols = covariate_cols,
              outcome_col = "Y", ps_method = ps_method, binary_method = binary_method,
              outcome_method = outcome_method, estimator = estimator, cross_fit = cross_fit,
              n_folds = n_folds, divergence = divergence, verbose = FALSE, xgb_params = xgb_params,
              nuisance_cache = cache, nuisance_tag = oc, pattern_pos = i,
              max_qdelta_failure_rate = max_qdelta_failure_rate, qdelta_tol = qdelta_tol,
              skip_null_shift = skip_null_shift, qdelta_memo = memo),
            error = function(e) {
              if (inherits(e, "nuisance_cache_mismatch")) stop(e)
              list(success = FALSE, error_message = conditionMessage(e),
                   solver_failure = inherits(e, "qdelta_solver_failure"))
            })
          if (isTRUE(res$success)) {
            # keep only the pattern's own influence values (the full-length vector is rebuilt on assembly)
            res$eif_k <- res$eif_values_extended[pat$obs_indices]
            res$eif_values_extended <- NULL
          }
          have[[r]] <- res
          saved[[run_key[r]]] <- list(signature = signature(oc, sc), result = res)
        }
        rm(memo)
      }
      rm(cache); invisible(gc(verbose = FALSE))
      if (!is.null(file)) {
        tmp <- paste0(file, ".", Sys.getpid(), ".tmp")
        saveRDS(list(version = 3L, units = saved), tmp)
        ok <- suppressWarnings(file.rename(tmp, file)) || isTRUE(file.copy(tmp, file, overwrite = TRUE))
        unlink(tmp)
        if (!ok) stop("could not write the pattern file ", file, call. = FALSE)
      }
    }
    unit_results[[i]] <- have
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    timing[i, ] <- list(i, pat$n_obs, pat$n_eligible, secs, any(todo))
    if (verbose) message(sprintf("  pattern %d/%d (n=%d, %d eligible): %s, %.0f s", i, n_pat, pat$n_obs,
                                 pat$n_eligible, if (any(todo)) sprintf("%d run(s) computed", sum(todo)) else "reloaded",
                                 secs))
  }

  # ---- assemble every run from the per-pattern results ------------------------------------------
  n_full <- nrow(data)
  out <- stats::setNames(vector("list", length(run_key)), run_key)
  for (r in seq_along(run_key)) {
    pr <- lapply(seq_len(n_pat), function(i) {
      u <- unit_results[[i]][[r]]
      if (isTRUE(u$success)) {
        u$eif_values_extended <- extend_eif_to_full(u$eif_k, u$obs_indices, n_full)
      }
      u
    })
    succeeded <- vapply(pr, function(x) isTRUE(x$success), logical(1))
    out[[r]] <- tryCatch({
      if (!any(succeeded)) stop("all patterns failed: ", paste(unique(vapply(pr, function(x) as.character(x$error_message), "")), collapse = "; "), call. = FALSE)
      agg <- suppressWarnings(aggregate_stratified_results(pr, patterns_to_use, confidence_level = confidence_level))
      n_estimated <- sum(vapply(pr[succeeded], function(x) x$n_obs, numeric(1)))
      n_stratified <- sum(patterns_to_use$n_obs)
      realized <- if (n_stratified > 0) n_estimated / n_stratified else 0
      if (realized < min_realized_coverage) {
        dropped <- patterns_to_use[!succeeded, , drop = FALSE]
        stop(sprintf("Realized pattern coverage %.1f%% is below min_realized_coverage %.1f%%: %d of %d patterns failed (%s)",
                     100 * realized, 100 * min_realized_coverage, sum(!succeeded), length(succeeded),
                     paste(sprintf("'%s' n=%s: %s", dropped$pattern_string, dropped$n_obs,
                                   vapply(pr[!succeeded], function(x) as.character(x$error_message), "")),
                           collapse = "; ")), call. = FALSE)
      }
      list(success = TRUE, estimate = agg$estimate, se = agg$se, ci_lower = agg$ci_lower,
           ci_upper = agg$ci_upper, p_value = agg$p_value, eif_values = agg$eif_values,
           pattern_summaries = agg$pattern_summaries, n_patterns_used = agg$n_patterns_used,
           coverage = agg$coverage, realized_coverage = realized, n_patterns_failed = sum(!succeeded),
           n_qdelta_fallback = sum(vapply(pr[succeeded], function(x) as.integer(x$n_qdelta_fallback %||% 0L), integer(1))),
           pattern_contributions = agg$pattern_contributions, outcome = runs$outcome[r],
           scenario = runs$scenario[r], delta = deltas[[runs$scenario[r]]], n = n_full, m = m)
    }, error = function(e) list(success = FALSE, error_message = conditionMessage(e),
                                outcome = runs$outcome[r], scenario = runs$scenario[r]))
    rm(pr)
  }
  list(runs = out, eligibility_patterns = pattern_info, patterns_to_use = patterns_to_use,
       feasibility_check = feasibility_check, timing = timing)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
