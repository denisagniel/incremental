# Test helpers for the explicit-column-specification interface.
#
# The package no longer auto-detects columns; every estimation entry point
# requires explicit treatment_cols / covariate_cols. In tests the data are
# constructed with the canonical A1..Am / X1..Xp naming, so these helpers derive
# the column-name vectors for the call sites. (Deriving columns in TEST code is
# fine; the anti-pattern removed from the package was auto-detection inside the
# estimators, not test-side argument construction.)

A_cols <- function(data) grep("^A[0-9]+$", names(data), value = TRUE)
X_cols <- function(data) grep("^X[0-9]+$", names(data), value = TRUE)
elig_cols <- function(data) grep("^elig_A[0-9]+$", names(data), value = TRUE)
