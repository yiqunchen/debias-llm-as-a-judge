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
# Install dependencies
install.packages(c("tidyverse", "future", "furrr", "quadprog"))

# Install package
devtools::install_local(".")
```

## Usage

The main simulation script is `simulation_llm_vs_ppi.R`, which runs a full comparison across the parameter grid and generates diagnostic plots.

```r
source("simulation_llm_vs_ppi.R")
```

Results are saved to a timestamped directory under `results/`.

## Package Functions

Core estimators in `R/sim_estimators.R`:

- `llm_point_and_ci()`: Measurement error corrected estimator with finite-sample adjustments
- `ppi_point_and_ci()`: Standard PPI estimator
- `ppi_pp_point_and_ci_general()`: PPI++ with optimized lambda
- `generate_dgp_data()`: Simulate Bernoulli outcomes with controllable prevalence

Additional utilities for multiclass settings and covariate-adjusted estimation are also provided.

## License

MIT
