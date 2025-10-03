# RSV Package Vignette: Poverty Application

## Overview

A comprehensive vignette has been created for the `rsv` R package that demonstrates how to replicate the poverty application figures from the paper. The vignette provides step-by-step guidance on using the package to estimate treatment effects with remotely sensed variables.

## What Was Created

### 1. Vignette File
- **Location**: `rsv/vignettes/poverty-application.Rmd`
- **Format**: R Markdown vignette
- **Length**: ~500 lines of documented code and explanation

### 2. Example Data
- **Location**: `rsv/inst/extdata/poverty_pca.csv`
- **Size**: ~790 KB, 8,312 observations
- **Variables**:
  - `Ycons`, `Ylowinc`, `Ymidinc`: Binary poverty outcomes
  - `D`: Treatment indicator (1 = treated, 0 = control, NA = holdout)
  - `wave`: Experimental design (Treatment, Control, Buffer, Holdout)
  - `R_PC1`: First principal component of satellite features (remotely sensed variable)
  - `clusters`: Cluster identifiers (subdistricts) for bootstrap

### 3. Updated Package Documentation
- Added vignette section to `rsv/README.md`
- Updated `rsv/DESCRIPTION` to include `fixest` in Suggests for benchmark comparisons

## What the Vignette Demonstrates

The vignette walks through a complete analysis workflow:

### 1. Data Preparation
- Loading the poverty dataset from Muralidharan et al. (2017)
- Understanding the experimental design (Treatment, Control, Buffer, Holdout villages)
- Creating sample indicators for RSV estimation

### 2. RSV Estimation - Real Sample Definition
- Experimental sample (S_e): Treatment + Control + Buffer villages
- Observational sample (S_o): Holdout + Buffer villages
- Buffer villages create overlap between samples
- Estimates treatment effects for all three poverty outcomes

### 3. RSV Estimation - Synthetic Sample Definition
- Experimental sample (S_e): All villages with treatment assignment
- Observational sample (S_o): Random 50% of clusters
- Demonstrates methodological flexibility in sample construction

### 4. Benchmark Estimation
- Traditional regression approach using only experimental sample
- Provides comparison point for RSV estimates

### 5. Visualization
- Creates treatment effect plots with 95% confidence intervals
- Replicates the figures from `outputs/poverty/realD/`:
  - `poverty_realD_te_Ycons.jpeg`
  - `poverty_realD_te_Ylowinc.jpeg`
  - `poverty_realD_te_Ymidinc.jpeg`

### 6. Technical Details
- Cluster-robust bootstrap standard errors
- Sample overlap structure analysis
- Interpretation of results

## How to Use the Vignette

### Option 1: Build and View the Vignette (After Installing Package)

```r
# Install the package
devtools::install("rsv", build_vignettes = TRUE)

# View the vignette
vignette("poverty-application", package = "remoteoutcome")

# Or browse all vignettes
browseVignettes("remoteoutcome")
```

### Option 2: Knit the Vignette Directly

```r
# From R, in the rsv_replication-main directory
rmarkdown::render("rsv/vignettes/poverty-application.Rmd")

# This will create an HTML file you can open in a browser
```

### Option 3: Run the Code Chunks Interactively

```r
# Load the package functions
devtools::load_all("rsv")

# Load the data
data_path <- "rsv/inst/extdata/poverty_pca.csv"
poverty_data <- read.csv(data_path)

# Follow the code chunks from the vignette...
```

## Key Functions Demonstrated

### Main RSV Estimation

```r
rsv_estimate(
  Y = poverty_data$Y_obs,          # Outcome (with NAs for missing)
  D = poverty_data$D,              # Treatment indicator
  S = poverty_data$S,              # Sample indicator ("e", "o", or "both")
  R = poverty_data$R_PC1,          # Remotely sensed variable
  method = "none",                 # No sample splitting
  se = TRUE,                       # Compute standard errors
  clusters = poverty_data$clusters, # Cluster identifiers
  B = 500,                         # Bootstrap replications
  alpha = 0.05,                    # 95% confidence interval
  seed = 42                        # Reproducibility
)
```

### Key Parameters Explained

- **Y**: Binary outcome variable (1 = poor, 0 = not poor). Set to NA where not observed.
- **D**: Binary treatment (1 = treated, 0 = control). Set to NA where not observed.
- **S**: Sample membership:
  - `"e"`: Experimental sample only (observe D, not Y)
  - `"o"`: Observational sample only (observe Y, not D)
  - `"both"`: In both samples (observe both Y and D)
- **R**: Remotely sensed variable(s) - here, PC1 of satellite features
- **method**:
  - `"none"`: Use all data for training and testing
  - `"split"`: Random sample split
  - `"crossfit"`: K-fold cross-fitting (recommended for larger datasets)
- **clusters**: Used for cluster-robust bootstrap inference

## Expected Output

The vignette produces:

### Treatment Effect Estimates Table
```
=== Consumption Poverty (Ycons) ===
          estimator    coef      se      lci     uci
         Benchmark -0.0746  0.0447  -0.1482 -0.0010
RSV: Synthetic...  -0.0340  0.0182  -0.0640 -0.0040
   RSV: Real...    -0.0127  0.0099  -0.0291  0.0037
```

### Treatment Effect Plots
- Forest plots showing estimates and confidence intervals
- Three estimators: Benchmark, RSV (Synthetic), RSV (Real)
- One plot for each poverty outcome

## Notes

1. **Bootstrap Replications**: The vignette uses B=500 for speed. For publication-quality results, use B=1000 or more.

2. **Data Source**: The poverty_pca.csv file is derived from:
   - SHRUG village-level data
   - SECC poverty measures
   - Muralidharan et al. (2017) experimental design
   - VIIRS nighttime lights
   - Mosaiks satellite image features
   - PCA transformation to create R_PC1

3. **Comparison with Original Results**: The vignette should produce results similar (but not identical due to bootstrap variation) to:
   - `data/poverty/processed/realD_coefs.csv`
   - Figures in `outputs/poverty/realD/`

4. **Package Installation**: The vignette requires these packages:
   - `rsv` (this package)
   - `dplyr`, `ggplot2` (data manipulation and plotting)
   - `fixest` (benchmark regression)

## Files Modified/Created

```
rsv/
├── inst/
│   └── extdata/
│       └── poverty_pca.csv          # Example data (NEW)
├── vignettes/
│   └── poverty-application.Rmd      # Main vignette (NEW)
├── DESCRIPTION                       # Updated with fixest in Suggests
└── README.md                         # Added vignette section
```

## Citation

If you use this vignette or the rsv package, please cite:

```bibtex
@article{rambachan2025program,
  title={Program Evaluation with Remotely Sensed Outcomes},
  author={Rambachan, Ashesh and Singh, Rahul and Viviano, Davide},
  journal={arXiv preprint arXiv:2411.10959},
  year={2025}
}

@article{muralidharan2017general,
  title={General Equilibrium Effects of (Improving) Public Employment Programs: Experimental Evidence from India},
  author={Muralidharan, Karthik and Niehaus, Paul and Sukhtankar, Sandip},
  journal={NBER Working Paper 23838},
  year={2017}
}
```

## Support

For questions or issues with the vignette:
- Check the package documentation: `?rsv_estimate`
- View the README: `rsv/README.md`
- Report issues: https://github.com/asheshrambachan/remoteoutcome/issues
