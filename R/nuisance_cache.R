#' ============================================================
#' Nuisance-model-fit caching for stratified estimation
#' ============================================================
#'
#' The eligibility-stratified pipeline
#' (\code{\link{estimate_qdelta_mean_stratified}}) is called once per
#' (cost outcome x delta scenario) cell. Inside each call it loops over
#' eligibility patterns and fits, per pattern, a cross-fitted propensity model
#' and a cross-fitted outcome model. Neither nuisance depends on \code{delta},
#' and the propensity nuisance does not depend on the outcome (cost type) at
#' all. In the real application (4 cost types x 14 delta scenarios) that means
#' every pattern's propensity fit is repeated 56 times and every pattern's
#' outcome fit 14 times, for identical inputs.
#'
#' This file implements an in-memory cache that removes that redundancy:
#' propensity fits are shared across ALL cells, outcome fits across the delta
#' scenarios of one cost type. The cache stores only the small fitted
#' \emph{matrices} (n_k x 2^m_k propensity and outcome-mean matrices), never the
#' bulk dataset, so the caller's "one cost-type dataset resident at a time"
#' memory discipline is preserved.
#'
#' @section Reference semantics:
#' The cache is a classed \strong{environment}, deliberately. Two properties
#' follow and both are load-bearing:
#' \itemize{
#'   \item Partial fills survive an error. Per-pattern estimation is wrapped in
#'     \code{tryCatch} by \code{estimate_pattern_effect()}; with value semantics
#'     a pattern that failed after its propensity fit would discard that fit.
#'   \item Large matrices are not copied when an entry is added.
#' }
#'
#' @section Reproducibility (per-fit namespaced seeds):
#' A cached fit must be bit-identical to a freshly computed one regardless of
#' how much RNG the surrounding loop happened to consume. Fold assignment for
#' cross-fitting is drawn with \code{sample()}, so each fit is executed inside
#' \code{withr::with_seed()} with a seed derived from the cache's
#' \code{seed_base}, the fit kind, the cost-type index and the pattern's
#' position in the pattern loop:
#' \preformatted{
#'   propensity seed = seed_base + 10000 + pattern_pos
#'   outcome    seed = seed_base + 20000 + 1000 * cost_type_index + pattern_pos
#' }
#' Two consequences: (1) seeds do not depend on iteration-order-consumed RNG,
#' so a run in which the propensity fit is a cache hit still seeds its outcome
#' fit exactly as a fresh run would; (2) because \code{withr::with_seed()}
#' saves and restores \code{.Random.seed}, the cache is RNG-neutral to the
#' caller's stream (e.g. the per-run \code{set.seed(SEED_BASE + run_id)} in
#' \code{application/analysis/07_stratified_estimation.R}).
#'
#' @section Fingerprints:
#' Every entry stores a fingerprint of the rows it was fitted on (number of
#' observations plus the sum of the row indices) together with the matrix shape
#' and the nuisance configuration it was fitted under. A mismatch on lookup is
#' a \strong{hard error}: silently returning a matrix aligned to different rows
#' would corrupt estimates in a way no downstream check would catch.
#'
#' @name nuisance_cache
#' @keywords internal
NULL

# Seed namespacing constants (see "Reproducibility" above). Fixed, not tunable:
# changing them changes every estimate, so they belong in one visible place.
NUISANCE_SEED_OFFSET_PS      <- 10000L  # propensity block
NUISANCE_SEED_OFFSET_OUTCOME <- 20000L  # outcome block base
NUISANCE_SEED_STRIDE_TAG     <- 1000L   # per-cost-type stride within outcome block

#' Create a nuisance-fit cache
#'
#' @description
#' Creates the classed environment used to share cross-fitted nuisance
#' estimates across calls to \code{\link{estimate_qdelta_mean_stratified}}.
#' Propensity matrices are keyed by eligibility pattern only (they do not
#' depend on the outcome); outcome-mean matrices are keyed by
#' (cost-type tag, eligibility pattern).
#'
#' @details
#' Intended usage in a (cost type x delta scenario) sweep:
#' \preformatted{
#' cache <- new_nuisance_cache(SEED_BASE, tags = cost_types)
#' for (ct in cost_types) {
#'   .clear_outcomes(cache)          # outcome fits are cost-type specific
#'   for (scen in scenarios) {
#'     estimate_qdelta_mean_stratified(..., nuisance_cache = cache,
#'                                     nuisance_tag = ct)
#'   }
#' }
#' }
#'
#' Internally the cached matrices are handed to
#' \code{\link{estimate_qdelta_mean}} through its existing
#' \code{oracle_propensities} / \code{oracle_outcomes} pre-supplied-nuisance
#' bypass. Those argument names predate this cache and are reused as-is: they
#' mean "nuisances supplied by the caller, skip fitting", which is exactly what
#' a cache hit needs.
#'
#' @param seed_base Integer scalar. Base for the per-fit seeds
#'   (see \code{\link{nuisance_cache}}). Use the same value that governs the
#'   surrounding sweep (e.g. \code{SEED_BASE}) so fits are reproducible.
#' @param tags Optional character vector of cost-type tags, in a FIXED order
#'   (e.g. the script's \code{cost_types}). Registering them up front makes the
#'   outcome-fit seeds depend only on the tag's position in this vector, so a
#'   crash-resumed sweep that skips earlier cost types still reproduces the same
#'   fits. An unregistered tag is then a hard error. If \code{NULL}, tags are
#'   registered in order of first use, which is only reproducible if the order
#'   of first use is.
#'
#' @return An environment of class \code{nuisance_cache}.
#'
#' @seealso \code{\link{.clear_outcomes}}, \code{\link{nuisance_cache_stats}},
#'   \code{\link{estimate_qdelta_mean_stratified}}
#'
#' @export
#' @examples
#' cache <- new_nuisance_cache(20260706L, tags = c("total", "inpatient"))
#' nuisance_cache_stats(cache)$n_ps_entries  # 0
new_nuisance_cache <- function(seed_base, tags = NULL) {
  if (!is.numeric(seed_base) || length(seed_base) != 1 || !is.finite(seed_base)) {
    stop("`seed_base` must be a single finite number.", call. = FALSE)
  }
  if (!is.null(tags)) {
    if (!is.character(tags) || length(tags) == 0 || anyNA(tags)) {
      stop("`tags` must be a non-empty character vector without NAs (or NULL).",
           call. = FALSE)
    }
    if (anyDuplicated(tags) > 0) {
      stop("`tags` must be unique; duplicated tags would collide in the cache.",
           call. = FALSE)
    }
  }

  cache <- new.env(parent = emptyenv())
  cache$seed_base  <- as.integer(seed_base)
  cache$tags       <- if (is.null(tags)) character(0) else as.character(tags)
  cache$tags_fixed <- !is.null(tags)

  # Nested (not flat) stores: outcome[[tag]][[pattern_string]], ps[[pattern_string]].
  # Nesting is what makes .clear_outcomes() a single assignment and keeps the
  # cost-type-independence of the propensity fits structurally obvious.
  cache$ps      <- list()
  cache$outcome <- list()

  # Diagnostics (cumulative across the whole sweep; not reset by .clear_outcomes).
  cache$n_ps_hit            <- 0L
  cache$n_ps_miss           <- 0L
  cache$n_ps_fit            <- 0L
  cache$n_outcome_hit       <- 0L
  cache$n_outcome_miss      <- 0L
  cache$n_outcome_fit       <- 0L
  cache$ps_fit_seconds      <- 0
  cache$outcome_fit_seconds <- 0
  cache$n_outcome_clears    <- 0L

  class(cache) <- c("nuisance_cache", "environment")
  cache
}

#' Check that an object is a usable nuisance cache
#'
#' @param cache Object to check.
#' @param arg Character, argument name for the error message.
#' @return Invisibly \code{TRUE}; errors otherwise.
#' @keywords internal
.check_nuisance_cache <- function(cache, arg = "nuisance_cache") {
  if (!inherits(cache, "nuisance_cache") || !is.environment(cache)) {
    stop("`", arg, "` must be created by new_nuisance_cache().", call. = FALSE)
  }
  invisible(TRUE)
}

#' Resolve a cost-type tag to its (deterministic) index
#'
#' @description
#' Maps a cost-type tag to the integer index used in the outcome-fit seed. When
#' the cache was created with a fixed \code{tags} vector the index is that
#' vector's position, so it is independent of which cost types a resumed run
#' actually recomputes. Otherwise tags are registered on first use.
#'
#' @param cache A \code{nuisance_cache}.
#' @param tag Character scalar tag, or \code{NULL} for the single-outcome case
#'   (treated as \code{"default"}).
#' @return Integer index (>= 1).
#' @keywords internal
.cache_tag_index <- function(cache, tag) {
  .check_nuisance_cache(cache)
  tag <- .cache_tag_name(tag)
  idx <- match(tag, cache$tags)
  if (!is.na(idx)) return(as.integer(idx))

  if (cache$tags_fixed) {
    stop("nuisance cache: tag '", tag, "' is not one of the registered tags (",
         paste(cache$tags, collapse = ", "), "). Registered tags fix the ",
         "outcome-fit seeds; an unregistered tag would silently change them.",
         call. = FALSE)
  }
  cache$tags <- c(cache$tags, tag)
  length(cache$tags)
}

#' Normalize a cache tag
#' @param tag Character scalar or NULL.
#' @return Character scalar.
#' @keywords internal
.cache_tag_name <- function(tag) {
  if (is.null(tag)) return("default")
  if (!is.character(tag) || length(tag) != 1 || is.na(tag)) {
    stop("`nuisance_tag` must be a single non-NA character string (or NULL).",
         call. = FALSE)
  }
  tag
}

#' Per-fit seeds for one pattern
#'
#' @description
#' Returns the namespaced seeds for a pattern's propensity and outcome fits.
#' The two fits are seeded SEPARATELY and each fit is executed inside
#' \code{withr::with_seed()}: a run whose propensity fit is a cache hit must
#' still seed its outcome fit exactly where a fresh run would, otherwise
#' resume-identity breaks.
#'
#' @param cache A \code{nuisance_cache}.
#' @param tag Cost-type tag (or NULL).
#' @param pattern_pos Integer, 1-based position of the pattern in the pattern
#'   loop. Stable because patterns are ordered by prevalence.
#' @return List with integer elements \code{ps} and \code{outcome}.
#' @keywords internal
.cache_seeds <- function(cache, tag, pattern_pos) {
  .check_nuisance_cache(cache)
  if (!is.numeric(pattern_pos) || length(pattern_pos) != 1 ||
      !is.finite(pattern_pos) || pattern_pos < 1) {
    stop("`pattern_pos` must be a single positive integer when a nuisance cache ",
         "is supplied (it namespaces the per-fit seeds).", call. = FALSE)
  }
  pattern_pos <- as.integer(pattern_pos)
  if (pattern_pos >= NUISANCE_SEED_STRIDE_TAG) {
    stop("`pattern_pos` (", pattern_pos, ") must be < ", NUISANCE_SEED_STRIDE_TAG,
         ": larger values would collide with the next cost type's seed block.",
         call. = FALSE)
  }
  tag_index <- .cache_tag_index(cache, tag)
  list(
    ps = cache$seed_base + NUISANCE_SEED_OFFSET_PS + pattern_pos,
    outcome = cache$seed_base + NUISANCE_SEED_OFFSET_OUTCOME +
      NUISANCE_SEED_STRIDE_TAG * tag_index + pattern_pos
  )
}

#' Build a cache-entry fingerprint
#'
#' @description
#' Cheap identity of the rows and configuration a nuisance matrix was fitted
#' on. Row identity uses \code{length(obs_indices)} and
#' \code{sum(obs_indices)}: patient identifiers are not available as a vector
#' here, and this pair changes for essentially any realistic re-slicing of the
#' cohort while costing one pass over an integer vector.
#'
#' @param obs_indices Integer vector of row indices in the full dataset.
#' @param n_combinations Integer, expected number of matrix columns (2^m_k).
#' @param config Named list of configuration that changes the fit (methods,
#'   cross-fitting, folds, covariates, outcome column).
#' @return List used with \code{identical()} for comparison.
#' @keywords internal
.nuisance_fingerprint <- function(obs_indices, n_combinations, config) {
  list(
    n_obs = length(obs_indices),
    idx_sum = sum(as.numeric(obs_indices)),
    n_combinations = as.integer(n_combinations),
    config = config
  )
}

#' Retrieve a cached nuisance matrix
#'
#' @description
#' Returns the cached matrix for \code{(kind, tag, pattern_string)}, or
#' \code{NULL} on a miss. If an entry exists but its fingerprint differs from
#' \code{fingerprint}, this is a hard error: the cache is being asked for a
#' matrix that does not describe the rows/configuration the caller is working
#' with, and returning it would misalign the estimate silently.
#'
#' @param cache A \code{nuisance_cache}.
#' @param kind Either \code{"ps"} or \code{"outcome"}.
#' @param tag Cost-type tag; ignored (and may be NULL) for \code{kind = "ps"},
#'   since propensity fits do not depend on the outcome.
#' @param pattern_string Character, eligibility-pattern key (e.g. \code{"101101"}).
#' @param fingerprint Fingerprint from \code{\link{.nuisance_fingerprint}}, or
#'   \code{NULL} to skip the check (not recommended).
#' @return Matrix on a hit, \code{NULL} on a miss.
#' @keywords internal
.cache_get <- function(cache, kind, tag, pattern_string, fingerprint = NULL) {
  .check_nuisance_cache(cache)
  kind <- .cache_check_kind(kind)
  pattern_string <- .cache_check_pattern(pattern_string)

  entry <- if (kind == "ps") {
    cache$ps[[pattern_string]]
  } else {
    tag_name <- .cache_tag_name(tag)
    store <- cache$outcome[[tag_name]]
    if (is.null(store)) NULL else store[[pattern_string]]
  }

  if (is.null(entry)) {
    if (kind == "ps") {
      cache$n_ps_miss <- cache$n_ps_miss + 1L
    } else {
      cache$n_outcome_miss <- cache$n_outcome_miss + 1L
    }
    return(NULL)
  }

  if (!is.null(fingerprint) && !identical(entry$fingerprint, fingerprint)) {
    msg <- paste0(
      "nuisance cache fingerprint mismatch for ", kind, " / pattern '",
      pattern_string, "'",
      if (kind == "outcome") paste0(" / tag '", .cache_tag_name(tag), "'") else "",
      ".\n  Cached: ", .fingerprint_str(entry$fingerprint),
      "\n  Requested: ", .fingerprint_str(fingerprint),
      "\n  The cached matrix does not describe these rows/configuration. ",
      "Reusing it would silently misalign the estimate. Create a fresh cache ",
      "(or call .clear_outcomes()) when the data or nuisance configuration changes.")
    # Signalled with a dedicated condition class so per-pattern error handling
    # cannot demote a cache misalignment to a "this pattern failed" result: it is
    # a programming/data error and must abort the run.
    stop(structure(
      class = c("nuisance_cache_mismatch", "error", "condition"),
      list(message = msg, call = NULL)
    ))
  }

  if (kind == "ps") {
    cache$n_ps_hit <- cache$n_ps_hit + 1L
  } else {
    cache$n_outcome_hit <- cache$n_outcome_hit + 1L
  }
  entry$value
}

#' Store a nuisance matrix in the cache
#'
#' @inheritParams .cache_get
#' @param value Matrix (n_k x 2^m_k) to store.
#' @param fingerprint Fingerprint from \code{\link{.nuisance_fingerprint}}.
#' @param fit_seconds Numeric, wall-clock seconds the fit took (diagnostics).
#' @return Invisibly \code{value}.
#' @keywords internal
.cache_put <- function(cache, kind, tag, pattern_string, value, fingerprint,
                       fit_seconds = 0) {
  .check_nuisance_cache(cache)
  kind <- .cache_check_kind(kind)
  pattern_string <- .cache_check_pattern(pattern_string)
  if (!is.matrix(value)) {
    stop("`value` must be a matrix (n_k x 2^m_k).", call. = FALSE)
  }
  if (nrow(value) != fingerprint$n_obs || ncol(value) != fingerprint$n_combinations) {
    stop("nuisance cache: matrix shape (", nrow(value), " x ", ncol(value),
         ") disagrees with its fingerprint (", fingerprint$n_obs, " x ",
         fingerprint$n_combinations, ").", call. = FALSE)
  }

  entry <- list(value = value, fingerprint = fingerprint)

  if (kind == "ps") {
    cache$ps[[pattern_string]] <- entry
    cache$n_ps_fit <- cache$n_ps_fit + 1L
    cache$ps_fit_seconds <- cache$ps_fit_seconds + as.numeric(fit_seconds)
  } else {
    tag_name <- .cache_tag_name(tag)
    if (is.null(cache$outcome[[tag_name]])) cache$outcome[[tag_name]] <- list()
    cache$outcome[[tag_name]][[pattern_string]] <- entry
    cache$n_outcome_fit <- cache$n_outcome_fit + 1L
    cache$outcome_fit_seconds <- cache$outcome_fit_seconds + as.numeric(fit_seconds)
  }
  invisible(value)
}

#' Drop all cached outcome fits (keep propensity fits)
#'
#' @description
#' Called at each cost-type transition of a sweep. Outcome models depend on the
#' outcome column, so their cached matrices are invalid for the next cost type;
#' propensity models depend only on treatments and covariates, which are
#' identical across the cost-type datasets (same cohort, different Y), so
#' \code{ps} is left untouched and its fits are shared across all cost types.
#'
#' @param cache A \code{nuisance_cache}.
#' @return Invisibly the number of outcome entries dropped.
#'
#' @export
#' @examples
#' cache <- new_nuisance_cache(1L)
#' .clear_outcomes(cache)  # no-op on an empty cache
.clear_outcomes <- function(cache) {
  .check_nuisance_cache(cache)
  n_dropped <- sum(vapply(cache$outcome, length, integer(1)))
  cache$outcome <- list()
  cache$n_outcome_clears <- cache$n_outcome_clears + 1L
  invisible(as.integer(n_dropped))
}

#' Cache size, hit-count and fit-time diagnostics
#'
#' @description
#' Reports how much the cache is actually being used. \code{n_ps_fit} /
#' \code{n_outcome_fit} are the number of model fits that were really performed:
#' comparing them against \code{n_ps_hit} / \code{n_outcome_hit} shows the
#' avoided redundancy directly.
#'
#' @param cache A \code{nuisance_cache}.
#' @return A one-row \code{data.frame} with entry counts, hit/miss/fit counts,
#'   cumulative fit seconds, and the approximate size of the stored matrices in MB.
#'
#' @export
#' @examples
#' cache <- new_nuisance_cache(1L)
#' nuisance_cache_stats(cache)
nuisance_cache_stats <- function(cache) {
  .check_nuisance_cache(cache)
  n_ps_entries <- length(cache$ps)
  n_outcome_entries <- sum(vapply(cache$outcome, length, integer(1)))
  # Stored values are double matrices, so 8 bytes/element is the size to within
  # a negligible header (avoids depending on utils::object.size).
  n_elements <- sum(vapply(cache$ps, function(e) as.numeric(length(e$value)), numeric(1))) +
    sum(unlist(lapply(cache$outcome, function(store) {
      vapply(store, function(e) as.numeric(length(e$value)), numeric(1))
    }), use.names = FALSE))
  size_mb <- if (length(n_elements) == 0) 0 else 8 * n_elements / 1024^2
  data.frame(
    n_ps_entries        = as.integer(n_ps_entries),
    n_outcome_entries   = as.integer(n_outcome_entries),
    n_ps_hit            = cache$n_ps_hit,
    n_ps_miss           = cache$n_ps_miss,
    n_ps_fit            = cache$n_ps_fit,
    n_outcome_hit       = cache$n_outcome_hit,
    n_outcome_miss      = cache$n_outcome_miss,
    n_outcome_fit       = cache$n_outcome_fit,
    ps_fit_seconds      = cache$ps_fit_seconds,
    outcome_fit_seconds = cache$outcome_fit_seconds,
    n_outcome_clears    = cache$n_outcome_clears,
    size_mb             = if (n_ps_entries + n_outcome_entries == 0) 0 else size_mb,
    stringsAsFactors    = FALSE
  )
}

#' @param x A \code{nuisance_cache}.
#' @param ... Ignored.
#' @rdname nuisance_cache_stats
#' @export
print.nuisance_cache <- function(x, ...) {
  s <- nuisance_cache_stats(x)
  cat("<nuisance_cache>\n")
  cat(sprintf("  seed_base: %d | tags: %s%s\n", x$seed_base,
              if (length(x$tags) == 0) "(none yet)" else paste(x$tags, collapse = ", "),
              if (x$tags_fixed) " [fixed]" else ""))
  cat(sprintf("  propensity: %d entries | %d fits, %d hits (%.1fs fitting)\n",
              s$n_ps_entries, s$n_ps_fit, s$n_ps_hit, s$ps_fit_seconds))
  cat(sprintf("  outcome:    %d entries | %d fits, %d hits (%.1fs fitting), %d clears\n",
              s$n_outcome_entries, s$n_outcome_fit, s$n_outcome_hit,
              s$outcome_fit_seconds, s$n_outcome_clears))
  cat(sprintf("  stored matrices: %.2f MB\n", s$size_mb))
  invisible(x)
}

#' Validate a cache kind argument
#' @param kind Character.
#' @return The validated kind.
#' @keywords internal
.cache_check_kind <- function(kind) {
  if (!is.character(kind) || length(kind) != 1 || !(kind %in% c("ps", "outcome"))) {
    stop("`kind` must be either 'ps' or 'outcome'.", call. = FALSE)
  }
  kind
}

#' Validate a pattern-string key
#' @param pattern_string Character.
#' @return The validated key.
#' @keywords internal
.cache_check_pattern <- function(pattern_string) {
  if (!is.character(pattern_string) || length(pattern_string) != 1 ||
      is.na(pattern_string) || !nzchar(pattern_string)) {
    stop("`pattern_string` must be a single non-empty character key.", call. = FALSE)
  }
  pattern_string
}

#' Fetch (or fit and store) one pattern's nuisance matrices
#'
#' @description
#' Cache-aware nuisance estimation for a single eligibility pattern. On a hit
#' the stored matrix is returned and no model is fitted; on a miss the model is
#' fitted with exactly the calls \code{\link{estimate_qdelta_mean}} would have
#' made internally, wrapped in \code{withr::with_seed()} with the namespaced
#' per-fit seed for this pattern.
#'
#' The two fits are seeded independently on purpose: if the propensity fit is a
#' cache hit and the outcome fit a miss, the outcome fit must still use the seed
#' it would have used in a fresh run.
#'
#' @param cache A \code{nuisance_cache}.
#' @param tag Cost-type tag (character scalar or NULL).
#' @param pattern_string Character, eligibility-pattern key.
#' @param pattern_pos Integer, 1-based pattern position in the pattern loop.
#' @param obs_indices Integer vector, this pattern's row indices in the full data.
#' @param data_k Data frame, this pattern's rows with treatments renamed to
#'   consecutive \code{A1..A_mk} (i.e. exactly what is passed to
#'   \code{estimate_qdelta_mean}).
#' @param treatment_cols Character vector of (renamed) treatment columns.
#' @param covariate_cols,outcome_col Column specifications.
#' @param ps_method,binary_method,outcome_method Nuisance methods.
#' @param cross_fit,n_folds Cross-fitting configuration ACTUALLY used for this
#'   pattern (i.e. after any per-pattern downgrade of \code{cross_fit}).
#' @param xgb_params Optional xgboost parameter list.
#'
#' @return List with \code{propensity_scores} (n_k x 2^m_k) and
#'   \code{outcome_means} (n_k x 2^m_k) matrices.
#' @keywords internal
.pattern_nuisances <- function(cache, tag, pattern_string, pattern_pos, obs_indices,
                               data_k, treatment_cols, covariate_cols,
                               outcome_col = "Y",
                               ps_method = "glm", binary_method = "glm",
                               outcome_method = "glm",
                               cross_fit = TRUE, n_folds = 5, xgb_params = NULL) {
  .check_nuisance_cache(cache)
  seeds <- .cache_seeds(cache, tag, pattern_pos)
  n_comb <- 2L^length(treatment_cols)

  if (nrow(data_k) != length(obs_indices)) {
    stop("nuisance cache: pattern data has ", nrow(data_k), " rows but ",
         length(obs_indices), " observation indices were supplied.", call. = FALSE)
  }

  # --- Propensity: keyed by pattern only (independent of the outcome) --------
  ps_fp <- .nuisance_fingerprint(obs_indices, n_comb, list(
    ps_method = ps_method, binary_method = binary_method,
    cross_fit = cross_fit, n_folds = n_folds,
    treatment_cols = treatment_cols, covariate_cols = covariate_cols
  ))
  ps_mat <- .cache_get(cache, "ps", NULL, pattern_string, ps_fp)
  if (is.null(ps_mat)) {
    t0 <- Sys.time()
    ps_est <- withr::with_seed(seeds$ps, estimate_propensity_scores(
      data_k,
      treatment_cols = treatment_cols,
      covariate_cols = covariate_cols,
      method = ps_method,
      binary_method = binary_method,
      cross_fit = cross_fit,
      n_folds = n_folds,
      xgb_params = xgb_params
    ))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ps_mat <- ps_est$propensity_scores
    .cache_put(cache, "ps", NULL, pattern_string, ps_mat, ps_fp, secs)
    rm(ps_est)
    gc(verbose = FALSE)  # mirrors estimate_qdelta_mean's post-fit memory discipline
  }

  # --- Outcome: keyed by (cost type, pattern) -------------------------------
  out_fp <- .nuisance_fingerprint(obs_indices, n_comb, list(
    outcome_method = outcome_method, cross_fit = cross_fit, n_folds = n_folds,
    treatment_cols = treatment_cols, covariate_cols = covariate_cols,
    outcome_col = outcome_col
  ))
  out_mat <- .cache_get(cache, "outcome", tag, pattern_string, out_fp)
  if (is.null(out_mat)) {
    t0 <- Sys.time()
    outcome_est <- withr::with_seed(seeds$outcome, estimate_outcome_models(
      data_k,
      treatment_cols = treatment_cols,
      covariate_cols = covariate_cols,
      outcome_col = outcome_col,
      method = outcome_method,
      cross_fit = cross_fit,
      n_folds = n_folds,
      xgb_params = xgb_params
    ))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out_mat <- outcome_est$outcome_means
    .cache_put(cache, "outcome", tag, pattern_string, out_mat, out_fp, secs)
    rm(outcome_est)
    gc(verbose = FALSE)
  }

  list(propensity_scores = ps_mat, outcome_means = out_mat)
}

#' Render a fingerprint compactly for error messages
#' @param fp Fingerprint list.
#' @return Character scalar.
#' @keywords internal
.fingerprint_str <- function(fp) {
  if (is.null(fp)) return("(none)")
  sprintf("n_obs=%s, idx_sum=%s, n_combinations=%s, config=[%s]",
          format(fp$n_obs), format(fp$idx_sum), format(fp$n_combinations),
          paste(names(fp$config),
                vapply(fp$config, function(v) paste(as.character(v), collapse = "+"),
                       character(1)),
                sep = "=", collapse = ", "))
}
