# incremental 0.1.3

## New features

* `estimate_qdelta_stratified_grid()` estimates eligibility-stratified effects
  pattern by pattern across a grid of deltas and outcomes, reusing nuisance fits
  and checkpointing completed patterns.
* `stratified_contrast_vs_observed()` contrasts a stratified fit with the
  observed mean without re-estimation.
* `estimate_qdelta_mean()` gains `qdelta_memo`, an environment that reuses the
  feasibility check and Q_delta solves across calls with identical propensity
  matrix, delta, divergence and tolerance (e.g. one call per outcome).
* `estimate_qdelta_mean_stratified()` gains `skip_null_shift` (opt-in): when every
  eligible delta in a pattern equals 1, return the sample mean and its influence
  function exactly, skipping nuisance fits.
* Stratified fits report `n_qdelta_fallback`.

## Fixes

* A duplicated `treatment_order` is now rejected rather than silently
  mis-evaluating the joint treatment distribution.
* `generate_binary_vectors()` validates `m` (positive integer scalar) and is
  limited to `m <= 20` (2^24 x 24 doubles would be about 3 GB).
* The feasibility shortcut no longer accepts targets within tolerance of the
  boundary; they fall through to the LP as before.

# incremental 0.1.2

## Breaking changes

* Column auto-detection has been removed throughout the package. Estimation and
  nuisance functions no longer infer treatment, covariate, outcome, or
  eligibility columns from name patterns (`grep("^A[0-9]+$", ...)`, etc.).
  Callers must now pass columns explicitly:
  - `estimate_qdelta_mean()`, `estimate_qdelta_vanishing_threshold()`,
    `estimate_qdelta_adaptive()`, `compute_psi_delta()`,
    `estimate_outcome_models()` gain required `treatment_cols` and
    `covariate_cols` arguments and an `outcome_col` argument (default `"Y"`).
  - `estimate_propensity_scores()` gains required `treatment_cols` and
    `covariate_cols` arguments.
  - `estimate_qdelta_mean_stratified()` gains required `covariate_cols`; its
    `eligibility_cols` argument is now required (no NULL auto-detection). Per
    the eligibility-stratified design, treatment columns are still derived per
    pattern internally.
  - `identify_eligibility_patterns()` now requires `eligibility_cols`.
  A new internal helper `validate_columns()` checks that all named columns exist
  and fails loudly otherwise. Because the estimators now use exactly the named
  columns, stray columns in the data (e.g. treatments outside the analysis set)
  are ignored rather than silently changing the number of treatments — this
  fixes a class of "non-conformable arguments" errors that arose when extra
  `A`-numbered columns were present.

* The `estimator = "tmle"` option has been **removed**. The branch computed the
  same one-step EIF estimate as `estimator = "one_step"` under a different
  label -- it never performed a substitution step, so it was not TMLE. Callers
  passing `estimator = "tmle"` now get an error. Use `"one_step"`.
* `sensitivity_analysis()` is renamed to `delta_sweep()`. The old name still
  works but warns, and is reserved for a bounded-risk-ratio sensitivity
  analysis under unmeasured confounding (not yet implemented).
* `estimate_qdelta_sensitivity()` removed: it called
  `estimate_propensity_scores_sequential()`, which never existed, and its own
  TODO conceded that nuisance reuse across the delta grid was unimplemented.
* `estimate_qdelta_mean_with_noisy_true_nuisance()` and
  `compute_noisy_true_nuisance()` removed, superseded by the
  `oracle_propensities` / `oracle_outcomes` arguments of
  `estimate_qdelta_mean()`.
* Roxygen docs no longer describe the estimator as "doubly robust". It is
  rate-robust via a second-order product-rate remainder; both nuisances must
  be consistent (Theorem 1).

## Bug fixes

* EIF standard errors now use the `1/n` plug-in convention
  (`sqrt(mean(phi^2)/n)`) at every site, matching `inference_qdelta_mean()`.
  Reported SEs change by a factor of `sqrt((n-1)/n)`.
* `estimate_qdelta_mean()` divided the EIF variance by `nrow(data)` rather
  than the number of retained EIF values, understating the SE whenever
  non-finite EIF values were dropped.
* `bootstrap_qdelta_mean()` and the delta sweep passed `divergence`
  positionally into the `treatment_cols` slot of `estimate_qdelta_mean()`.
* `to_numeric_matrix()` now sets `colnames()` on its result, so
  `m[, "colname"]` works.
