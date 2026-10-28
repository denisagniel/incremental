#!/usr/bin/env Rscript
#' Debug Test Failures
#'
#' This script helps identify and debug failing tests by:
#' 1. Running tests with detailed output
#' 2. Running each test file individually to isolate failures
#' 3. Providing detailed error information

# Set working directory to package root
package_root <- "package/incremental"
setwd(package_root)

# Load required packages
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat package is required. Install it with: install.packages('testthat')")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools package is required. Install it with: install.packages('devtools')")
}

# Load the package
cat("Loading package in development mode...\n")
devtools::load_all()

# Get list of test files
test_dir <- "tests/testthat"
test_files <- list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE)

cat("\n=== Found", length(test_files), "test files ===\n")
for (f in test_files) {
  cat("  -", basename(f), "\n")
}

# Function to run a single test file with detailed output
run_test_file <- function(test_file) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Running:", basename(test_file), "\n")
  cat(rep("=", 70), "\n", sep = "")
  
  tryCatch({
    # Use test_file for detailed output
    result <- testthat::test_file(test_file, reporter = testthat::default_reporter())
    
    # Print summary
    if (length(result) > 0) {
      passed <- sum(sapply(result, function(x) x$results[[1]]$passed))
      failed <- sum(sapply(result, function(x) x$results[[1]]$failed))
      skipped <- sum(sapply(result, function(x) x$results[[1]]$skipped))
      
      cat("\nSummary for", basename(test_file), ":\n")
      cat("  Passed:", passed, "\n")
      cat("  Failed:", failed, "\n")
      cat("  Skipped:", skipped, "\n")
      
      if (failed > 0) {
        cat("\n⚠️  FAILURES DETECTED ⚠️\n")
        # Print failure details
        for (test_result in result) {
          for (result_item in test_result$results) {
            if (inherits(result_item, "expectation_failure") || 
                inherits(result_item, "expectation_error")) {
              cat("\nTest:", test_result$test, "\n")
              cat("Error:", result_item$message, "\n")
              if (!is.null(result_item$srcref)) {
                cat("Location:", capture.output(print(result_item$srcref)), "\n")
              }
            }
          }
        }
      }
    }
    
    return(list(passed = passed, failed = failed, skipped = skipped))
  }, error = function(e) {
    cat("\n❌ ERROR running test file:", e$message, "\n")
    return(list(passed = 0, failed = 1, skipped = 0))
  })
}

# Run all tests with detailed output
cat("\n\n=== Running all tests with detailed output ===\n")
cat("This will show which specific tests fail and why.\n\n")

all_results <- list()
for (test_file in test_files) {
  results <- run_test_file(test_file)
  all_results[[basename(test_file)]] <- results
}

# Summary
cat("\n\n", rep("=", 70), "\n", sep = "")
cat("=== FINAL SUMMARY ===\n")
cat(rep("=", 70), "\n", sep = "")

total_passed <- sum(sapply(all_results, function(x) x$passed))
total_failed <- sum(sapply(all_results, function(x) x$failed))
total_skipped <- sum(sapply(all_results, function(x) x$skipped))

cat("\nTotal across all test files:\n")
cat("  Passed:", total_passed, "\n")
cat("  Failed:", total_failed, "\n")
cat("  Skipped:", total_skipped, "\n")

if (total_failed > 0) {
  cat("\n⚠️  FAILING TEST FILES:\n")
  for (test_name in names(all_results)) {
    if (all_results[[test_name]]$failed > 0) {
      cat("  -", test_name, ":", all_results[[test_name]]$failed, "failure(s)\n")
    }
  }
  
  cat("\n💡 TIPS FOR DEBUGGING:\n")
  cat("  1. Review the error messages above for each failing test\n")
  cat("  2. Run individual test files: testthat::test_file('tests/testthat/test-XXX.R')\n")
  cat("  3. Run a specific test: testthat::test_that('test name', { ... })\n")
  cat("  4. Check the source code for the functions being tested\n")
  cat("  5. Use debug() to step through failing functions\n")
} else {
  cat("\n✅ All tests passed!\n")
}

