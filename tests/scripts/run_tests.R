#!/usr/bin/env Rscript
#' Run all tests in the incremental package
#'
#' This script runs all testthat tests for the incremental package.
#' You can run it from the command line or source it in R.

# Set the working directory to the package root
package_root <- "package/incremental"

# Check if testthat is available
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat package is required. Install it with: install.packages('testthat')")
}

# Load the package (assuming it's installed or we're in dev mode)
# If using devtools, you might want to use devtools::load_all() instead
if (requireNamespace("devtools", quietly = TRUE)) {
  cat("Using devtools to load package...\n")
  devtools::load_all(package_root)
} else {
  cat("Loading incremental package...\n")
  library(incremental)
}

# Run all tests
cat("\n=== Running all tests ===\n")
test_dir <- file.path(package_root, "tests", "testthat")

if (!dir.exists(test_dir)) {
  stop("Test directory not found: ", test_dir)
}

# Run tests
test_results <- testthat::test_dir(
  test_dir,
  reporter = testthat::default_reporter()
)

# Print summary
cat("\n=== Test Summary ===\n")
print(test_results)

# Return exit code (0 = success, 1 = failure)
if (any(test_results$failed > 0)) {
  quit(status = 1)
} else {
  quit(status = 0)
}

