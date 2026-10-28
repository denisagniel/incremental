#' @file constants.R
#' @title Package-Wide Constants
#' @description
#' Defines numerical tolerances, thresholds, and default parameters used
#' throughout the incremental package. Centralizes magic numbers for
#' maintainability and consistency.
#' @author Denis Agniel

# Numerical tolerances for convergence
QDELTA_TOLERANCE <- 1e-8
QDELTA_TOLERANCE_RELAXED <- 1e-6
QDELTA_TOLERANCE_CHI2 <- 1e-8  # Chi2 is linear (h(x)=x/2+1): Newton converges in 1 step
CONSTRAINT_TOL <- 1e-8

# Propensity score thresholds
PROPENSITY_ACTIVE_THRESHOLD <- 1e-12  # Minimum for "active" combinations
PROPENSITY_MIN_THRESHOLD <- 1e-8      # Clipping threshold for division safety
PROPENSITY_SUM_TOL <- 1e-6            # Tolerance for sum-to-one check

# Lambda initialization parameters
REVERSE_KL_LAMBDA0_INIT <- -0.1       # Initial lambda_0 for reverse KL
KL_LAMBDA0_INIT <- 0.0                # Initial lambda_0 for KL
CHI2_LAMBDA0_INIT <- 0.0              # Initial lambda_0 for chi-squared
HELLINGER_LAMBDA0_INIT <- 0.0         # Initial lambda_0 for Hellinger

# Chi-squared boundary handling
CHI2_BOUNDARY_LAMBDA0 <- -2           # Special lambda_0 for near-boundary
CHI2_BOUNDARY_THRESHOLD <- 0.01       # Threshold for "near boundary"
CHI2_MIN_ITER_FOR_RELAXATION <- 50    # Min iterations before relaxing tolerance

# Numerical safeguards
REVERSE_KL_LINEAR_TERM_MIN <- -1e12   # Floor for linear terms (h(x) = -1/x, domain x < 0)
OMEGA_REGULARIZATION_EPS <- 1e-4      # Regularization for near-singular Omega

# EIF chunk memory model (see eif_default_chunk_size)
EIF_LIVE_MATRICES <- 14L              # n_chunk x 2^m matrices alive at peak in eif_compute_chunk
EIF_CHUNK_BUDGET_BYTES <- 2e9         # Total working-matrix budget across all workers
EIF_MIN_CHUNK <- 25L                  # Floor on chunk size (the old floor of 1000 broke the budget at m >= 12)
EIF_MAX_CHUNK <- 50000L               # Ceiling on chunk size

# Minimum observations for modeling
MIN_OBS_COMBINATION <- 10             # Min obs per treatment combination for modeling
MIN_OBS_TRAIN <- 5                    # Min obs for training a model

# Convergence warnings
MIN_CONVERGENCE_RATE <- 0.95          # Warn if convergence rate below this

# Adaptive constraint dropping
DEFAULT_MAX_CONSTRAINT_DROPS <- 2     # Max constraints to drop per iteration
DEFAULT_CONDITION_THRESHOLD <- 1e10   # Condition number threshold for singularity

# Default optimization parameters
DEFAULT_MAX_ITER <- 200               # Default maximum Newton-Raphson iterations

# Numerical differentiation
NUMERICAL_DIFF_STEP <- 1e-6           # Step size for numerical Jacobian approximation

# Domain boundary safeguards for h functions
# These prevent division by zero and domain violations in divergence computations
HELLINGER_DOMAIN_BUFFER <- 1e-12      # Buffer for linear_term < 1 check (Hellinger: h(t)=1/(1-t)²)
REVERSE_KL_DOMAIN_BUFFER <- 1e-12     # Buffer for linear_term < 0 check (reverse KL: h(x)=-1/x)
CHI2_DOMAIN_BUFFER <- 1e-12           # Buffer for chi-squared computations

# Near-zero detection thresholds
NEAR_ZERO_THRESHOLD <- 1e-10          # General threshold for "effectively zero"

# Strict domain validation (for reject_if_invalid checks)
DOMAIN_VALIDATION_STRICT_TOL <- 1e-8  # Stricter tolerance for domain checks

# EIF centering check
EIF_CENTERING_THRESHOLD <- 0.1        # Warn if |mean(EIF)| > 0.1 * sd(EIF)

# Minimum safe values (prevent division by zero in normalization)
MIN_SAFE_VALUE <- 1e-8                # Minimum safe value for denominators

# Newton solver perturbation for stalled convergence
NEWTON_PERTURBATION_MAGNITUDE <- 0.01 # Magnitude for escaping local minima

# Chi-squared Levenberg-Marquardt damping parameters
CHI2_LM_DAMPING_STRONG <- 1.0         # Strong damping when constraint norm > 1
CHI2_LM_DAMPING_MEDIUM <- 0.1         # Medium damping when 0.1 < constraint norm < 1
CHI2_LM_DAMPING_WEAK <- 0.01          # Weak damping when constraint norm < 0.1
CHI2_SWITCH_TO_LM_ITER <- 20          # Iteration to switch to Levenberg-Marquardt
CHI2_MAX_STEP_SIZE <- 0.5              # Maximum step size for chi2 (trust region)

# KL divergence overflow protection (Phase 1 fix for 381K convergence)
KL_LINEAR_TERM_MAX <- 700              # exp(700) ≈ 1e304, exp(710) = Inf

# Timeout parameters (Phase 1 fix for 381K convergence)
QDELTA_DEFAULT_TIMEOUT_SECS <- 300    # 5 minutes per observation (conservative)
QDELTA_LARGE_N_TIMEOUT_SECS <- 600    # 10 minutes for large datasets

# Adaptive tolerance for large datasets (Phase 1 fix for 381K convergence)
QDELTA_LARGE_N_THRESHOLD <- 100000    # Threshold for "large dataset" (progress reporting only)
QDELTA_LARGE_N_TOL <- 1e-6            # Reference only: no longer applied automatically; pass qdelta_tol to relax

# Line search parameters (Phase 1 fix for 381K convergence)
LINE_SEARCH_MAX_STEPS <- 20           # Increased from 10
LINE_SEARCH_ARMIJO_C1 <- 1e-4         # Armijo condition parameter
