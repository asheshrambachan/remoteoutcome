#!/usr/bin/env Rscript
# =============================================================================
# Package Build Script
# =============================================================================
# This script builds the rsv package, generates documentation, and runs checks.
# Run this script from the rsv/ directory or its parent directory.

# Determine package directory
if (basename(getwd()) == "rsv") {
  pkg_dir <- getwd()
} else if (file.exists("rsv")) {
  pkg_dir <- "rsv"
} else {
  stop("Cannot find rsv package directory. Run this from rsv/ or its parent.")
}

cat("Building rsv package at:", pkg_dir, "\n\n")

# Install required packages for building
required_pkgs <- c("devtools", "roxygen2", "testthat")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg)
  }
}

library(devtools)
library(roxygen2)

# Step 1: Generate documentation from roxygen2 comments
cat("Step 1: Generating documentation...\n")
roxygen2::roxygenise(pkg_dir)
cat("Documentation generated in man/\n\n")

# Step 2: Check package
cat("Step 2: Checking package...\n")
check_results <- devtools::check(pkg_dir)
cat("\n")

# Step 3: Build package
cat("Step 3: Building package...\n")
pkg_file <- devtools::build(pkg_dir)
cat("Package built:", pkg_file, "\n\n")

# Step 4: Install package
cat("Step 4: Installing package...\n")
devtools::install(pkg_dir)
cat("Package installed!\n\n")

cat("========================================\n")
cat("Package build complete!\n")
cat("========================================\n")
cat("\nTo load the package, run: library(remoteoutcome)\n")
cat("To view help, run: ?rsv_estimate\n")
