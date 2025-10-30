# Fact and Friction: Analyses

This repository contains to date code and materials for the Fact and Friction project. 

---

## Overview Power Analysis

The simulation study evaluating statistical power to detect prespecified effects in continuous-time structural equation models (ctsem) (see proposal).

**Goal:**
Estimate power for detecting a known dynamic effect (`gamma_r1 = 0.1`) in a coupled latent-process model using `ctsem`.

### Quick start

**Clone the repository**

```r 
git clone https://github.com/da-wi/fact_friction.git
cd fact_friction
```

**Run the main simulation script**

```r 
Rscript scripts/00_power_analysis.R
```

Example output file:

`results/power/power_results_raw.20251030_083015.rds`

**Inspect or summarize results**
Open `scripts/01_power_assessment.R` in RStudio. 

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

