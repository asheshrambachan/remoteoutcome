# Installation and Build Instructions

## Quick Installation

From R, run:

```r
# Install devtools if needed
install.packages("devtools")

# Install remoteoutcome package
devtools::install_local("/Users/asheshr/Desktop/rsv_replication-main/rsv")
```

## Building from Source

### Option 1: Using the Build Script (Recommended)

From the terminal, navigate to the package directory and run:

```bash
cd /Users/asheshr/Desktop/rsv_replication-main/rsv
Rscript build_package.R
```

This script will:
1. Generate documentation from roxygen2 comments
2. Check the package for errors
3. Build the package tarball
4. Install the package

### Option 2: Manual Build

From R, run:

```r
library(devtools)
library(roxygen2)

# Set package directory
pkg_dir <- "/Users/asheshr/Desktop/rsv_replication-main/rsv"

# Generate documentation
roxygen2::roxygenise(pkg_dir)

# Check package
devtools::check(pkg_dir)

# Build and install
devtools::install(pkg_dir)
```

### Option 3: Using R CMD

From terminal:

```bash
cd /Users/asheshr/Desktop/rsv_replication-main

# Generate documentation (from R)
R -e "roxygen2::roxygenise('rsv')"

# Build package
R CMD build rsv

# Check package
R CMD check rsv_0.1.0.tar.gz

# Install package
R CMD INSTALL rsv_0.1.0.tar.gz
```

## Verifying Installation

After installation, verify the package works:

```r
library(remoteoutcome)

# View package help
?rsv_estimate

# Run a simple example
set.seed(123)
n <- 500
R <- rnorm(n)
Y <- rbinom(n, 1, plogis(0.5 * R))
D <- rbinom(n, 1, 0.5)
S <- sample(c("e", "o"), n, replace = TRUE)
Y[S == "e"] <- NA
D[S == "o"] <- NA

result <- rsv_estimate(Y = Y, D = D, S = S, R = R,
                       method = "crossfit", se = FALSE)
print(result)
```

## Troubleshooting

### Missing Dependencies

If you get errors about missing packages, install them:

```r
install.packages(c("dplyr", "randomForest", "boot", "fixest"))
```

### Documentation Not Generated

If documentation is not generated, ensure roxygen2 is installed:

```r
install.packages("roxygen2")
roxygen2::roxygenise("/Users/asheshr/Desktop/rsv_replication-main/rsv")
```

### Package Check Warnings

Common warnings and how to fix them:

1. **"Non-standard files/directories found at top level"**: Add files to `.Rbuildignore`
2. **"Undocumented arguments"**: Check that all function parameters have `@param` documentation
3. **"Missing Rd tags"**: Ensure all exported functions have complete roxygen2 documentation

## Development Workflow

For active development:

```r
# Load package for testing without installing
devtools::load_all("/Users/asheshr/Desktop/rsv_replication-main/rsv")

# Run tests (if you add them)
devtools::test("/Users/asheshr/Desktop/rsv_replication-main/rsv")

# Check package
devtools::check("/Users/asheshr/Desktop/rsv_replication-main/rsv")
```
