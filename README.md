# debiasLLMReporting

R package for comparing methods that correct for measurement error when using LLM-as-a-judge classifiers to estimate population prevalence.

## Overview

When using large language models (LLMs) as classifiers to estimate the prevalence of some outcome in a population, the LLM's imperfect accuracy introduces bias. This package provides simulation tools and estimators to study and correct for this measurement error.

We compare three main approaches:

1. **Measurement Error Correction**: Classical approach using the Rogan-Gladen estimator, which corrects for misclassification using estimated sensitivity (q1) and specificity (q0) from a calibration sample.

2. **Prediction-Powered Inference (PPI)**: Uses a small labeled dataset to debias predictions from a larger unlabeled dataset, combining human labels with LLM predictions.

3. **PPI++**: An extension of PPI that optimizes a tuning parameter (lambda) to minimize variance while maintaining valid confidence intervals.

## Key Features

- Simulation framework for comparing estimator performance across:
  - Varying true prevalence (theta)
  - Different classifier accuracy levels (q0, q1)
  - Multiple labeling budgets (1% to 50% labeled)
  - Different calibration sample allocations (m0:m1 ratios)

- Metrics evaluated:
  - Coverage of confidence intervals
  - CI width
  - Bias and percent bias
  - Calibration estimate quality

- Comprehensive visualization suite with faceted plots

## Installation

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# Install the package and its dependencies
devtools::install_local(".", dependencies = TRUE)
```

## Usage

```r
# Install the package
devtools::install_local(".")

# Load the package
library(debiasLLMReporting)
```

## Reproducibility

The `reproduce/` directory contains scripts to replicate our analyses:

- `simulation_llm_vs_ppi.R` - Monte Carlo simulation comparing estimators across parameter grid
- `apply_estimators.R` - Apply estimators to real LLM judge evaluation data

Results are saved to timestamped directories under `results/`.

## Package Functions

Core estimators:

- `llm_point_and_ci()` - Rogan-Gladen estimator with delta method CI
- `ppi_point_and_ci()` - Standard PPI estimator
- `ppi_pp_point_and_ci_general()` - PPI++ with optimized lambda
- `eif_point_and_ci()` - Efficient influence function estimator
- `fit_misclass_mle()` - Joint MLE for misclassification model

## License

MIT
