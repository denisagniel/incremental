# incremental (development version)

## Breaking changes

* Column auto-detection has been removed throughout the package. Estimation and
  nuisance functions no longer infer treatment, covariate, outcome, or
  eligibility columns from name patterns (`grep("^A[0-9]+$", ...)`, etc.).
  Callers must now pass columns explicitly:
  - `estimate_qdelta_mean()`, `estimate_qdelta_vanishing_threshold()`,
    `estimate_qdelta_adaptive()`, `estimate_qdelta_sensitivity()`,
    `compute_psi_delta()`, `estimate_outcome_models()`,
    `estimate_qdelta_mean_with_noisy_true_nuisance()` gain required
    `treatment_cols` and `covariate_cols` arguments and an `outcome_col`
    argument (default `"Y"`).
  - `estimate_propensity_scores()` and `compute_noisy_true_nuisance()` gain
    required `treatment_cols` and `covariate_cols` arguments.
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
