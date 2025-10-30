# Fact and Friction: Analyses

This repository contains to date code and materials for the simulation study evaluating statistical power to detect cross-lagged effects in continuous-time structural equation models (ctsem).

The project simulates longitudinal data under known parameters, fits a predefined continuous-time model using Stan, and estimates the power to detect the true effect size across varying sample sizes.

---

## Overview Power Analysis

**Goal:**
Estimate power for detecting a known dynamic effect (`gamma_r1 = 0.1`) in a coupled latent-process model using `ctsem`.

**Pipeline summary:**

1. Simulate synthetic continuous-time data (`simulate_ct_data.R`)
2. Fit the ctsem model to each dataset (`fit_once.R`)
3. Repeat across sample sizes and replicates (`main.R`)
4. Aggregate and save results to timestamped `.rds` files in `results/power/`

---

## Requirements

**System prerequisites**

- R (≥ 4.3)
- C++ toolchain for compiling Stan models (e.g., RTools on Windows, `g++` on macOS/Linux)

**Core R packages**

| Package | Purpose |
|----------|----------|
| `ctsem` | Continuous-time SEM estimation |
| `tidyverse` | Data wrangling and plotting |
| `future`, `future.apply`, `progressr` | Parallel processing and progress display |
| `RhpcBLASctl` | BLAS/OMP thread control |
| `tictoc` | Timing messages |
| `rstan` | Backend for model compilation and optimization |

Install all dependencies via:

```r
install.packages(c(
  "ctsem", "tidyverse", "future", "future.apply",
  "progressr", "RhpcBLASctl", "tictoc", "rstan"
))


## Reproducing the Simulation

**Clone the repository**

```r 
git clone https://github.com/da-wi/fact_friction.git
cd fact_friction
```

**Run the main simulation script**
```r 
Rscript scripts/00_power_analysis.R
```

This will:
* compile or load the pre-saved Stan model (`model/ctmodel_compiled.rds`)
* simulate data across all conditions (`Ns`, `reps`, `true_effect`)
* run parallel fits using 4 workers (make sure you have enough memory and cores)
* save a timestamped results file in `results/power/`

Example output file:

`results/power/power_results_raw.20251030_083015.rds`

**Inspect or summarize results**
```r
res <- readRDS("results/power/power_results_raw.20251030_083015.rds")
res |>
  tidyr::unnest(out) |>
  dplyr::group_by(N) |>
  dplyr::summarise(power = mean(detected, na.rm = TRUE))
```

## Notes 
* All random seeds are explicitly set for reproducibility.
* Parallelization is handled through the future package (plan(multisession, workers = 4)).
* Each simulation run automatically saves results with a timestamped filename to prevent overwriting.
* The Stan model is compiled once and reused for efficiency.

## Expected runtime

Approximate runtime depends on system performance:
| Sample size | Replicates | CPU cores | Estimated time |
| ----------- | ---------- | --------- | -------------- |
| 200–400     | 100        | 4         | ~1–2 hours     |
