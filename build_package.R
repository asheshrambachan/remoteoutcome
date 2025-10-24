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

# Step 2: Build manual (PDF)
cat("\nStep 2: Building manual...\n")
manual_file <- devtools::build_manual(pkg_dir, path = pkg_dir)
cat("Manual built\n")


# Step 3: Check package
cat("\nStep 3: Checking package...\n")
check_results <- devtools::check(pkg_dir)

# Step 4: Build package
cat("\nStep 4: Building package and vignettes...\n")
pkg_file <- devtools::build(pkg_dir)
cat("Package built:", pkg_file, "\n")

# Step 5: Install package
cat("\nStep 5: Installing package...\n")
devtools::install(pkg_dir)
cat("Package installed!\n")


# Step 6: Create site 
create_site <- FALSE
if (create_site){
  required_pkgs <- c("pkgdown", "usethis")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      install.packages(pkg)
    }
  }
  
  cat("\nStep 6: Create site ...\n")
  usethis::use_pkgdown_github_pages()
  pkgdown::build_site(pkg_dir)
}

cat("========================================\n")
cat("Package build complete!\n")
cat("========================================\n")
cat("\nTo load the package, run: library(remoteoutcome)\n")
cat("To view help, run: ?rsv_estimate\n")