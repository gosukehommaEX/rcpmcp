# rcpmcp

**Regional Consistency Probability under Multiple Comparison Procedures**

## Overview

`rcpmcp` is an R package for computing regional consistency probabilities (RCPs) and disjunctive power in multi-regional clinical trials (MRCTs) with multiple primary endpoints under multiplicity adjustment procedures.

When multiple primary endpoints are employed, familywise error rate (FWER) control via the Bonferroni procedure leads to monotone inflation of the null RCP with increasing number of endpoints if unadjusted consistency thresholds are used. This package implements:

- **Closed-form RCP computation** via the multivariate normal distribution (Bonferroni)
- **Monte Carlo simulation** for Bonferroni, Holm, Hochberg, and Hommel procedures
- **Adjusted consistency thresholds** (Method 1 and Method 2) that correct for null RCP inflation
- **Sample size determination** under the Bonferroni procedure

Regional consistency is evaluated following the Japanese Ministry of Health, Labour and Welfare (MHLW, 2007) guidance:

- **Method 1**: Region 1 retains at least a fraction `gamma_M1` of the overall treatment effect
- **Method 2**: All regional estimates exceed the threshold `gamma_M2`

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("gosukehommaEX/rcpmcp")
```

## Main Functions

| Function | Description |
|---|---|
| `rcpmcp_single()` | RCP and power for a fixed number of endpoints K. Supports region-specific treatment effects and unknown-variance sensitivity analysis under `approach = "simulation"`. |
| `rcpmcp_multiple()` | RCP and power across K = 1, ..., K_max endpoints. Inherits the simulation extensions of `rcpmcp_single()`. |
| `rcpmcp_get_gamma()` | Adjusted consistency thresholds via root-finding |
| `ssmcp_single()` | Sample size for a fixed K (Bonferroni) |
| `ssmcp_multiple()` | Sample size across K = 1, ..., K_max (Bonferroni) |

## Examples

### 1. RCP for a fixed number of endpoints (closed-form)

```r
library(rcpmcp)

# K = 3 independent endpoints, closed-form solution (Bonferroni)
result <- rcpmcp_single(
  delta    = c(0.2, 0.2, 0.2),
  Sigma    = diag(3),
  N        = 200,
  fs       = c(0.1, 0.45, 0.45),
  K        = 3,
  gamma_M1 = 0.5,
  gamma_M2 = 0,
  alpha    = 0.025
)
print(result)
```

### 2. RCP via Monte Carlo simulation (all four MCP methods)

```r
result_sim <- rcpmcp_single(
  delta    = c(0.2, 0.2, 0.2),
  Sigma    = diag(3),
  N        = 200,
  fs       = c(0.1, 0.45, 0.45),
  K        = 3,
  gamma_M1 = 0.5,
  gamma_M2 = 0,
  alpha    = 0.025,
  approach = "simulation",
  nsim     = 10000,
  seed     = 1
)
print(result_sim)
```

### 3. Sample size determination and adjusted thresholds

```r
# Step 1: Determine sample size for each K = 1, ..., 5 (Bonferroni)
ss <- ssmcp_multiple(
  delta        = rep(0.2, 5),
  Sigma        = diag(5),
  fs           = c(0.1, 0.45, 0.45),
  K_max        = 5,
  alpha        = 0.025,
  target_power = 0.8
)
print(ss)

# Step 2: Compute adjusted consistency thresholds
gamma_res <- rcpmcp_get_gamma(
  Sigma    = diag(5),
  N        = ss$result$N,
  fs       = c(0.1, 0.45, 0.45),
  K_max    = 5,
  gamma_M1 = 0.5,
  gamma_M2 = 0,
  alpha    = 0.025
)
print(gamma_res)

# Step 3: Compare adjusted vs. unadjusted RCP across K
res_adj <- rcpmcp_multiple(
  delta    = rep(0.2, 5),
  Sigma    = diag(5),
  N        = ss$result$N,
  fs       = c(0.1, 0.45, 0.45),
  K_max    = 5,
  gamma_M1 = gamma_res$result$gamma_M1_adj,
  gamma_M2 = gamma_res$result$gamma_M2_adj,
  alpha    = 0.025
)

res_unadj <- rcpmcp_multiple(
  delta    = rep(0.2, 5),
  Sigma    = diag(5),
  N        = ss$result$N,
  fs       = c(0.1, 0.45, 0.45),
  K_max    = 5,
  gamma_M1 = 0.5,
  gamma_M2 = 0,
  alpha    = 0.025
)

# Step 4: Plot
plot(res_adj,
     overlay      = list(res_unadj),
     group_labels = c("Adjusted", "Unadjusted"))
```

### 4. Sensitivity analyses (new in v0.1.1)

The simulation branch of `rcpmcp_single()` (and `rcpmcp_multiple()`) supports two extensions for sensitivity analyses:

- **Region-specific treatment effects**: supply `delta` as a `K`-by-`S` matrix, where element `delta[k, s]` is the true treatment effect for endpoint `k` in region `s`. This allows evaluation of scenarios in which the regional effect departs from the global effect.
- **Unknown variance**: set `variance_known = FALSE` to switch the overall test from a Z-statistic (normal critical value) to a t-statistic with degrees of freedom `nu = N - 2`. The sample variance-covariance matrix is drawn from a Wishart distribution, avoiding the need to generate individual-level data.

The following example combines both extensions. Region 1 has zero treatment effect for all endpoints, while regions 2 and 3 have a treatment effect of 0.2. The variance is treated as unknown.

```r
# Region 1 is null; regions 2-3 have delta = 0.2 for all three endpoints
delta_mat <- matrix(c(0.0, 0.2, 0.2,
                      0.0, 0.2, 0.2,
                      0.0, 0.2, 0.2),
                    nrow = 3, ncol = 3, byrow = TRUE)

result_sens <- rcpmcp_single(
  delta          = delta_mat,
  Sigma          = diag(3),
  N              = 200,
  fs             = c(0.1, 0.45, 0.45),
  K              = 3,
  gamma_M1       = 0.5,
  gamma_M2       = 0,
  alpha          = 0.025,
  approach       = "simulation",
  nsim           = 10000,
  seed           = 1,
  variance_known = FALSE
)
print(result_sens)
```

Both extensions are only supported under `approach = "simulation"`. The closed-form (`approach = "formula"`) solution is preserved unchanged from earlier versions.

## Reproducibility

All figures and tables in the associated manuscript, including the Supporting
Information, can be reproduced from the scripts shipped with the package. To
keep heavy computation separate from figure rendering, each set of scripts is
split into a computation stage that caches intermediate results as `.rds` files
and a rendering stage that reads them.

| Location | Contents |
|---|---|
| `inst/scripts/` | Main-text Figures 1-4 and Table 1 |
| `inst/sup_info/` | Supporting Information Figures S1-S6 |

Within each folder, `data_generate_*.R` runs all sample-size, adjusted-threshold,
closed-form, and Monte Carlo computations and writes the results to `data/*.rds`,
and `table_and_figure_*.R` reads those files and writes the figures (EPS and PDF)
and the LaTeX table to `results/`.

To reproduce the main-text Figures 1-4 and Table 1, set the working directory to
`inst/scripts/` and run:

```r
library(rcpmcp)
source("data_generate_main.R")          # computation -> data/*.rds
source("table_and_figure_manuscript.R") # rendering   -> results/
```

The Supporting Information Figures S1-S6 are reproduced analogously from
`inst/sup_info/`:

```r
library(rcpmcp)
source("data_generate_supplement.R")
source("table_and_figure_supplement.R")
```

All paths are relative to the working directory, so the `data/` and `results/`
subfolders are created alongside the scripts. The number of Monte Carlo
iterations is set by `NSIM` near the top of each computation script (`1000000`
for the published results; reduce it for a quick run). Rendering additionally
requires the `scales`, `patchwork`, and `ggh4x` packages.

## License

GPL (>= 2)
