#!/usr/bin/env Rscript
# =============================================================================
# Package Build Script
# =============================================================================
# This script builds the remoteoutcome package, generates documentation, and runs checks.
# Run this script from the remoteoutcome/ directory or its parent directory.

# Determine package directory
if (basename(getwd()) == "remoteoutcome") {
  pkg_dir <- getwd()
} else if (file.exists("remoteoutcome")) {
  pkg_dir <- "remoteoutcome"
} else {
  stop("Cannot find remoteoutcome package directory. Run this from remoteoutcome/ or its parent.")
}

cat("Building remoteoutcome package at:", pkg_dir, "\n\n")

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
cat("Documentation generated in man/\n")

# Step 2: Check package
cat("\nStep 2: Checking package...\n")
check_results <- devtools::check(pkg_dir)

# Step 3: Build package
cat("\nStep 3: Building package and vignettes...\n")
pkg_file <- devtools::build(pkg_dir)
cat("Package built:", pkg_file, "\n")

# Step 4: Install package
cat("\nStep 4: Installing package...\n")
devtools::install(pkg_dir)
cat("Package installed!\n")

# # Step 5: Create site 
create_site <- FALSE
if (create_site){
  required_pkgs <- c("pkgdown", "usethis")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      install.packages(pkg)
    }
  }
  
  cat("\nStep 5: Create site ...\n")
  usethis::use_pkgdown_github_pages()
  pkgdown::build_site(pkg_dir)
}

cat("========================================\n")
cat("Package build complete!\n")
cat("========================================\n")
cat("\nTo load the package, run: library(remoteoutcome)\n")
cat("To view help, run: ?rsv_estimate\n")