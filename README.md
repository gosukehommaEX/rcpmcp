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
| `rcpmcp_single()` | RCP and power for a fixed number of endpoints K |
| `rcpmcp_multiple()` | RCP and power across K = 1, ..., K_max endpoints |
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

## Reproducibility

All figures and tables in the associated manuscript can be reproduced by executing:

```r
source(system.file("scripts/table_and_figure_manuscript.R", package = "rcpmcp"))
```

## License

GPL (>= 2)
