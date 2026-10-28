# incremental: Vector Incremental Treatment Effects

[![R-CMD-check](https://github.com/denisagniel/incremental/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/denisagniel/incremental/actions/workflows/R-CMD-check.yaml)

An R package implementing vector incremental treatment effects for causal inference with multiple binary treatments. The package provides doubly-robust estimation with efficient influence functions, support for multiple divergence measures (KL, reverse KL, f-divergences), and tools for sensitivity analysis.

This package accompanies the paper "Vector incremental treatment effects" (Agniel, Rubinstein, Han, Yu, Finnerty, Horvitz-Lennon, and Normand). See [Citation](#citation) below.

## Overview

Vector incremental treatment effects extend the concept of incremental propensity score interventions to settings with multiple binary treatments. The method allows researchers to estimate causal effects when realistic interventions may act on many treatments simultaneously, when synergies and interactions among treatments are known to occur, or when the true underlying treatment is difficult to describe comprehensively.

### Key Features

- **Multiple Divergence Measures**: Support for KL-divergence, reverse KL, and general f-divergences (χ², Hellinger, custom)
- **Doubly-Robust Estimation**: Plugin, one-step, and TMLE estimators with efficient influence functions
- **Comprehensive Inference**: Asymptotic and bootstrap confidence intervals, hypothesis testing
- **Flexible Nuisance Estimation**: GLM, SuperLearner, neural networks, and cross-fitting support
- **Sensitivity Analysis**: Bounds under a bounded counterfactual risk-ratio model for violations of ignorability

## Installation

```r
# install.packages("devtools")  # if needed
devtools::install_github("denisagniel/incremental")
```

### Requirements

- R (>= 4.0.0)
- Core dependencies (installed automatically): `dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `R6`, `MASS`, `lpSolve`

### Optional dependencies

For advanced nuisance-estimation methods, install as needed:

```r
install.packages(c("SuperLearner", "sl3"))   # SuperLearner integration
install.packages("nnet")                     # neural network propensity scores
install.packages("mgcv")                     # GAM outcome models
install.packages("randomForest")             # random forest models
install.packages("xgboost")                  # gradient boosting
install.packages("nloptr")                   # advanced optimization
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

# Estimate vector incremental treatment effects
result <- estimate_qdelta_mean(
  sim_data$data,
  delta = c(1.5, 2.0),   # Odds scaling factors
  divergence = "kl",     # KL-divergence
  estimator = "one_step" # One-step estimator
)

print(result$estimate)    # Point estimate
print(result$se)          # Standard error
print(result$ci_lower)    # Lower CI bound
print(result$ci_upper)    # Upper CI bound
print(result$p_value)     # P-value
```

### Advanced Usage

```r
# Use different divergence measures
kl_result <- estimate_qdelta_mean(sim_data$data, delta = c(1.5, 2.0), divergence = "kl")
chi2_result <- estimate_qdelta_mean(sim_data$data, delta = c(1.5, 2.0), divergence = "chi2")

# Use SuperLearner for nuisance estimation
sl_result <- estimate_qdelta_mean(
  sim_data$data,
  delta = c(1.5, 2.0),
  ps_method = "superlearner",
  outcome_method = "superlearner"
)

# Bootstrap inference
boot_result <- bootstrap_qdelta_mean(
  sim_data$data,
  delta = c(1.5, 2.0),
  n_bootstrap = 1000
)

# Sensitivity analysis across a grid of delta values
delta_grid <- list(
  c(1.2, 1.2), c(1.5, 1.5), c(2.0, 2.0),
  c(1.2, 1.5), c(1.5, 2.0), c(2.0, 1.2)
)
sens_result <- sensitivity_analysis(sim_data$data, delta_grid)
```

## Package Structure

### Core Functions

- **`compute_q_delta()`**: Main function to find the Q_delta distribution satisfying marginal odds constraints via minimum-divergence projection
- **`estimate_qdelta_mean()`**: Primary user-facing function for estimation
- **`efficient_influence_function()`**: Efficient influence function computation
- **`compute_psi_delta()`**: Calculate E[Y^{Q_delta}] using identification results

### Divergence Functions

- **`kl_divergence_qdelta()`**: Forward KL D_KL(Q || P_0) implementation
- **`reverse_kl_divergence_qdelta()`**: Reverse KL D_KL(P_0 || Q) implementation
- **`f_divergence_qdelta()`**: General f-divergence framework

### Nuisance Estimation

- **`estimate_propensity_scores()`**: Propensity score estimation with multiple methods
- **`estimate_outcome_models()`**: Outcome regression with various algorithms

### Inference

- **`inference_qdelta_mean()`**: Wald-type confidence intervals
- **`bootstrap_qdelta_mean()`**: Bootstrap inference
- **`hypothesis_test_qdelta_mean()`**: Hypothesis testing
- **`sensitivity_analysis()`**: Sensitivity analysis under a bounded risk-ratio model

## Theoretical Background

Vector incremental treatment effects extend incremental propensity score interventions to multiple binary treatments. The method:

1. **Defines causal effects** based on marginal odds scaling: P(A_k = 1 | X) / (1 - P(A_k = 1 | X)) = δ_k × P_0(A_k = 1 | X) / (1 - P_0(A_k = 1 | X))
2. **Finds the Q_delta distribution** that minimizes a divergence D(Q || P_0) subject to the marginal constraints, resolving the non-uniqueness of a joint distribution with specified marginals
3. **Estimates E[Y^{Q_delta}]** using doubly-robust methods with efficient influence functions
4. **Provides inference** using asymptotic theory or bootstrap methods, and avoids requiring full joint-treatment positivity

### Key Advantages

- **No full joint-treatment positivity assumption required**
- **Handles multiple treatments**: naturally extends to vector-valued treatments
- **Preserves dependence structure**: the minimum-divergence projection departs as little as possible from the observed joint law
- **Interpretable effects**: odds scaling has a clear policy interpretation
- **Doubly-robust estimation** with efficient influence functions

## Examples

### Example 1: Basic Estimation

```r
data <- simulate_vector_treatment_data(
  n = 1000, m = 2, p = 2,
  treatment_effects = c(0.5, 0.3),
  interaction_effects = matrix(c(0, 0.2, 0.2, 0), nrow = 2)
)

result <- estimate_qdelta_mean(data$data, delta = c(1.5, 1.5))
```

### Example 2: Divergence Comparison

```r
kl_result <- estimate_qdelta_mean(data$data, delta = c(1.5, 2.0), divergence = "kl")
chi2_result <- estimate_qdelta_mean(data$data, delta = c(1.5, 2.0), divergence = "chi2")
hellinger_result <- estimate_qdelta_mean(data$data, delta = c(1.5, 2.0), divergence = "hellinger")

estimates <- c(kl_result$estimate, chi2_result$estimate, hellinger_result$estimate)
names(estimates) <- c("KL", "Chi2", "Hellinger")
print(estimates)
```

### Example 3: Sensitivity Analysis

```r
delta_scenarios <- list(
  "moderate" = c(1.2, 1.2),
  "large" = c(1.5, 1.5),
  "mixed" = c(1.2, 1.5),
  "extreme" = c(2.0, 2.0)
)

sens_result <- sensitivity_analysis(data$data, delta_scenarios)
print(sens_result$summary)
```

## Testing

```r
devtools::test()
```

## Documentation

### Vignettes

- **Introduction**: Basic usage and examples (`vignette("introduction", package = "incremental")`)
- **Feasibility & optimization**: theory and diagnostics behind the minimum-divergence projection (`vignette("feasibility-optimization", package = "incremental")`)

### Function Documentation

```r
?estimate_qdelta_mean
?compute_q_delta
?efficient_influence_function
?incremental
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
