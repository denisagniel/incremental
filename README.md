# incremental: Vector Incremental Treatment Effects

[![R-CMD-check](https://github.com/denisagniel/incremental/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/denisagniel/incremental/actions/workflows/R-CMD-check.yaml)

An R package implementing vector incremental treatment effects for causal inference with multiple binary treatments. The package provides cross-fitted one-step estimation based on the efficient influence function (EIF), support for multiple divergence measures (KL, reverse KL, chi-squared, Hellinger), and tools for inference.

This package accompanies the paper "Vector incremental treatment effects" (Agniel, Rubinstein, Han, Yu, Finnerty, Horvitz-Lennon, and Normand). See [Citation](#citation) below. **The paper is the source of truth for what the method does and what properties it has; this README defers to it wherever the two could be read as disagreeing.**

## Overview

Vector incremental treatment effects extend the concept of incremental propensity score interventions to settings with multiple binary treatments. The method allows researchers to estimate causal effects when realistic interventions may act on many treatments simultaneously, when synergies and interactions among treatments are known to occur, or when the true underlying treatment is difficult to describe comprehensively.

### Key Features

- **Multiple Divergence Measures**: KL-divergence, reverse KL, chi-squared, and Hellinger. (The divergence set is a fixed, hard-coded four-way choice -- there is currently no way to supply an arbitrary custom divergence function.)
- **Efficient-Influence-Function Estimation**: A cross-fitted one-step estimator built from the EIF, with Wald confidence intervals. Nuisance error enters only through a second-order product-rate remainder, so a slowly-converging outcome model is tolerated when the propensity model converges quickly. **This is not double robustness**: the estimator is *not* consistent when the propensity model is inconsistent but the outcome model is correct, because the estimand depends on the propensity score through the divergence tilt and not only through inverse-probability weighting. Both nuisances must be consistent; see Theorem 1 of the accompanying paper.
- **Inference**: Asymptotic (EIF-based) confidence intervals, bootstrap confidence intervals, and a simple Wald hypothesis test.
- **Nuisance Estimation**: GLM (multinomial logistic regression via `nnet::multinom()`), random forest, xgboost, and SuperLearner, plus cross-fitting. (The `"nnet"` method name is an alias for the GLM path -- it is multinomial logistic regression, not a neural network.)
- **Delta Sensitivity Sweeps**: `sensitivity_analysis()` and `estimate_qdelta_sensitivity()` re-run estimation across a grid of delta values to show how the estimate changes with the intervention strength.

### Not yet implemented

The paper also derives bounds on $\psi(\bdelta)$ under a bounded counterfactual risk-ratio model for violations of ignorability (Proposition on sensitivity bounds, with EIF-based estimators for the bounds). **This sensitivity-bounds machinery is not yet implemented in the package** -- only the delta-sweep functions described above are. If you need the risk-ratio bounds, they must currently be implemented by hand from the paper's formulas.

### Experimental

An `estimator = "tmle"` option is implemented in `estimate_qdelta_mean()` but is **not covered by the accompanying paper's theory** -- no asymptotic normality or efficiency result is established for it, and the function emits a warning when used. Its targeting step fluctuates only the outcome model, while the estimand also depends on the propensity score through the divergence tilt, so it is not solving the full efficient score. Use `estimator = "one_step"` for any inference you intend to report; the `"tmle"` path is provided for exploration only.

## Installation

```r
# install.packages("devtools")  # if needed
devtools::install_github("denisagniel/incremental")
```

### Requirements

- R (>= 4.0.0)
- Dependencies (installed automatically): `dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `R6`, `MASS`, `lpSolve`, `withr`, `parallel`, `nnet`

### Optional dependencies

For additional nuisance-estimation methods, install as needed:

```r
install.packages("SuperLearner")   # SuperLearner integration
install.packages("mgcv")           # GAM outcome models
install.packages("randomForest")   # random forest models
install.packages("xgboost")        # gradient boosting
```

## Quick Start

```r
library(incremental)

# Generate example data
sim_data <- simulate_vector_treatment_data(
  n = 1000,
  m = 2,  # 2 treatments
  p = 2,  # 2 covariates
  treatment_effects = c(0.5, 0.3),
  interaction_effects = matrix(c(0, 0.2, 0.2, 0), nrow = 2),
  seed = 123
)
# Treatment and covariate column names follow the "A1..Am" / "X1..Xp" convention
# used throughout the package; estimate_qdelta_mean() requires them explicitly
# (there is no auto-detection).

# Estimate vector incremental treatment effects
result <- estimate_qdelta_mean(
  sim_data$data,
  delta = c(1.5, 2.0),            # Odds scaling factors
  treatment_cols = c("A1", "A2"),
  covariate_cols = c("X1", "X2"),
  divergence = "kl",               # KL-divergence
  estimator = "one_step"           # one-step EIF estimator (recommended)
)

print(result$estimate)    # Point estimate
print(result$se)          # Standard error
print(result$ci_lower)    # Lower CI bound
print(result$ci_upper)    # Upper CI bound
print(result$p_value)     # P-value
```

### Advanced Usage

```r
tcols <- c("A1", "A2")
xcols <- c("X1", "X2")

# Use different divergence measures
kl_result   <- estimate_qdelta_mean(sim_data$data, delta = c(1.5, 2.0), treatment_cols = tcols, covariate_cols = xcols, divergence = "kl")
chi2_result <- estimate_qdelta_mean(sim_data$data, delta = c(1.5, 2.0), treatment_cols = tcols, covariate_cols = xcols, divergence = "chi2")

# Use SuperLearner for nuisance estimation
sl_result <- estimate_qdelta_mean(
  sim_data$data,
  delta = c(1.5, 2.0),
  treatment_cols = tcols, covariate_cols = xcols,
  ps_method = "superlearner",
  outcome_method = "superlearner"
)

# Bootstrap inference
boot_result <- bootstrap_qdelta_mean(
  sim_data$data,
  delta = c(1.5, 2.0),
  treatment_cols = tcols, covariate_cols = xcols,
  n_bootstrap = 1000
)

# Sweep estimates across a grid of delta values
delta_grid <- list(
  c(1.2, 1.2), c(1.5, 1.5), c(2.0, 2.0),
  c(1.2, 1.5), c(1.5, 2.0), c(2.0, 1.2)
)
sens_result <- sensitivity_analysis(sim_data$data, delta_grid, treatment_cols = tcols, covariate_cols = xcols)
```

## Package Structure

### Core Functions

- **`estimate_qdelta_mean()`**: Primary user-facing function for estimation and inference
- **`compute_psi_delta()`**: Internal workhorse implementing outcome-regression, IPW, and EIF-based ("doubly_robust", a legacy name -- see Experimental/robustness notes above) point estimates; usually called via `estimate_qdelta_mean()`
- **`efficient_influence_function()`** / **`efficient_influence_function_active()`**: Efficient influence function computation
- **`compute_h_function_value()`** / internal `divergence_h_pair()` (`R/divergence_h_functions.R`): the four divergence-tilt functions (KL, reverse KL, chi-squared, Hellinger) and their derivatives, dispatched by the `divergence` argument. There is no separate `kl_divergence_qdelta()` / `f_divergence_qdelta()`-style function per divergence -- all four are handled by this single dispatcher.

### Nuisance Estimation

- **`estimate_propensity_scores()`**: Propensity score estimation (glm/nnet, random forest, xgboost, SuperLearner, sequential, mlr3)
- **`estimate_outcome_models()`**: Outcome regression with the same set of base learners

### Inference

- **`inference_qdelta_mean()`**: Wald-type confidence intervals from EIF values
- **`bootstrap_qdelta_mean()`**: Bootstrap inference
- **`hypothesis_test_qdelta_mean()`**: Simple z-test against a null value
- **`sensitivity_analysis()`** / **`estimate_qdelta_sensitivity()`**: Delta-grid sweeps (see "Not yet implemented" above for the risk-ratio bounds, which these do *not* cover)

## Theoretical Background

Vector incremental treatment effects extend incremental propensity score interventions to multiple binary treatments. The method:

1. **Defines causal effects** based on marginal odds scaling: P(A_k = 1 | X) / (1 - P(A_k = 1 | X)) = δ_k × P_0(A_k = 1 | X) / (1 - P_0(A_k = 1 | X))
2. **Finds the Q_delta distribution** that minimizes a divergence D(Q || P_0) subject to the marginal constraints, resolving the non-uniqueness of a joint distribution with specified marginals
3. **Estimates E[Y^{Q_delta}]** with a cross-fitted one-step estimator derived from the efficient influence function
4. **Provides inference** using asymptotic theory or bootstrap methods, and avoids requiring full joint-treatment positivity

### Key Advantages

- **No full joint-treatment positivity assumption required**
- **Handles multiple treatments**: naturally extends to vector-valued treatments
- **Preserves dependence structure**: the minimum-divergence projection departs as little as possible from the observed joint law
- **Interpretable effects**: odds scaling has a clear policy interpretation
- **Efficient and rate-robust**: attains the semiparametric efficiency bound; nuisance error enters only at second order, so the outcome model may converge slowly provided the propensity model converges quickly (both must be consistent -- this is not double robustness)

## Testing

```r
devtools::test()
```

## Documentation

### Vignettes

- **Introduction** (`vignette("introduction", package = "incremental")`): basic usage and examples
- **Feasibility & optimization** (`vignette("feasibility-optimization", package = "incremental")`): theory and diagnostics behind the minimum-divergence projection

### Function Documentation

```r
?estimate_qdelta_mean
?compute_psi_delta
?efficient_influence_function
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes with appropriate tests
4. Submit a pull request

```r
devtools::load_all()
devtools::test()
devtools::check()
devtools::document()
```

## Citation

If you use this package in your research, please cite the accompanying paper:

```
@article{agniel_vector_2026,
  title   = {Vector incremental treatment effects},
  author  = {Agniel, Denis and Rubinstein, Max and Han, Larry and Yu, Anthony
             and Finnerty, Molly T. and Horvitz-Lennon, Marcela V.
             and Normand, Sharon-Lise T.},
  year    = {2026},
  note    = {Preprint and code: \url{https://github.com/denisagniel/incremental}}
}
```

## License

MIT &copy; Denis Agniel. See [LICENSE](LICENSE) / [LICENSE.md](LICENSE.md) for details.

## Acknowledgments

This package implements the vector incremental treatment effects methodology described in the accompanying paper. The method extends the incremental propensity score interventions of Kennedy (2019) and Bonvini et al. (2023) to multiple binary treatments.

## References

- Kennedy, E. H. (2019). Nonparametric causal effects based on incremental propensity score interventions. *Journal of the American Statistical Association*, 114(526), 645-656.
- Bonvini, M., Kennedy, E. H., & Ventura, L. (2023). Incremental propensity score interventions for continuous treatments. *Biometrics*, 79(1), 1-12.
