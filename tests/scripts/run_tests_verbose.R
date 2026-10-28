#!/usr/bin/env Rscript
#' Run tests with verbose output to see detailed failure information
#'
#' Usage:
#'   source("run_tests_verbose.R")
#'   # or
#'   Rscript run_tests_verbose.R

# Navigate to package directory
setwd("package/incremental")

# Load devtools and package
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Please install devtools: install.packages('devtools')")
}

library(devtools)
load_all()

# Run tests with detailed reporter
cat("Running tests with detailed output...\n\n")
test_results <- test(
  reporter = testthat::default_reporter()  # Shows detailed output for failures
)

# The test() function will print detailed information about failures
# Check the console output for specific error messages

