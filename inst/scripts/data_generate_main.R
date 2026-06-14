# =============================================================================
# data_generate_main.R
#
# Computation stage for the manuscript figures and tables of:
#   "A Cautionary Note on Regional Consistency Evaluation with Multiple
#    Primary Endpoints in Multi-Regional Clinical Trials"
#
# This script performs all heavy numerical computation (sample sizes, adjusted
# thresholds, closed-form RCP, and Monte Carlo simulation) and writes the
# resulting intermediate data to "data/*.rds". The companion rendering script
# table_and_figure_manuscript.R reads these files and produces the figures and
# the LaTeX table without re-running any computation.
#
# Each Monte Carlo call self-seeds internally (set.seed(seed) inside the rcpmcp
# functions), so separating computation from rendering does not change any
# numerical result.
#
# Usage
#   Set the working directory to this script's location (inst/scripts), then
#   run the whole file. Intermediate results are written to "data/" relative
#   to the working directory. The package must be installed or loaded first:
#     library(rcpmcp)              # installed package, or
#     devtools::load_all()         # from the package root during development
#
# Outputs (relative to the working directory)
#   data/fig1_data.rds
#   data/fig2_data.rds
#   data/fig3_data.rds
#   data/fig4_data.rds
#   data/table1_data.rds
#   data/data_generate_main.log
#   data/sessionInfo_data_generate_main.txt
# =============================================================================

library(rcpmcp)

# --------------------------------------------------------------------------- #
# Output location (relative paths only)
# --------------------------------------------------------------------------- #
data_dir <- "data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_main.log")
cat(sprintf("data_generate_main.R log -- %s\n", format(Sys.time())),
    file = log_file, append = FALSE)

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  message(msg)
}

timed <- function(label, expr) {
  t <- system.time(val <- force(expr))["elapsed"]
  log_msg("%-28s %7.1f sec", label, t)
  val
}

# --------------------------------------------------------------------------- #
# Compute-side helpers
# --------------------------------------------------------------------------- #
make_sigma <- function(rho, k) {
  S <- matrix(rho, k, k)
  diag(S) <- 1
  S
}

.fill_na <- function(x) {
  for (i in seq_along(x)) {
    if (is.na(x[i]) && i > 1L) x[i] <- x[i - 1L]
  }
  x
}

# --------------------------------------------------------------------------- #
# Common parameters
# --------------------------------------------------------------------------- #
K_MAX <- 5
DELTA <- rep(0.2, K_MAX)
ALPHA <- 0.025
GM1_0 <- 0.5     # unadjusted gamma_M1
GM2_0 <- 0.0     # unadjusted gamma_M2
NSIM  <- 1000000 # Monte Carlo iterations (reduce to 10000 for quick checks)
SEED  <- 42

FS_01 <- c(0.1, 0.45, 0.45)   # f1 = 0.1
FS_02 <- c(0.2, 0.40, 0.40)   # f1 = 0.2

RHO_SEQ <- seq(0, 0.8, by = 0.2)

log_msg("Settings: K_MAX = %d, NSIM = %d, SEED = %d, alpha = %.3f",
        K_MAX, NSIM, SEED, ALPHA)

# =============================================================================
# Pre-compute: sample sizes and adjusted thresholds (rho = 0, f1 = 0.1 and 0.2)
# Used by Figures 1, 2, 4 and shared across builders
# =============================================================================

ss_01 <- timed("ss_01 (sample size)", ssmcp_multiple(
  delta = DELTA, Sigma = make_sigma(0, K_MAX),
  fs = FS_01, K_max = K_MAX, alpha = ALPHA, target_power = 0.8))

ss_02 <- timed("ss_02 (sample size)", ssmcp_multiple(
  delta = DELTA, Sigma = make_sigma(0, K_MAX),
  fs = FS_02, K_max = K_MAX, alpha = ALPHA, target_power = 0.8))

gm_01 <- timed("gm_01 (adjusted gamma)", suppressWarnings(
  rcpmcp_get_gamma(Sigma = make_sigma(0, K_MAX), N = ss_01$result$N,
                   fs = FS_01, K_max = K_MAX,
                   gamma_M1 = GM1_0, gamma_M2 = GM2_0, alpha = ALPHA)))

gm_02 <- timed("gm_02 (adjusted gamma)", suppressWarnings(
  rcpmcp_get_gamma(Sigma = make_sigma(0, K_MAX), N = ss_02$result$N,
                   fs = FS_02, K_max = K_MAX,
                   gamma_M1 = GM1_0, gamma_M2 = GM2_0, alpha = ALPHA)))

# Record NA positions before forward-fill (for honest plotting)
na_idx_02 <- which(is.na(gm_02$result$gamma_M2_adj))

# Forward-fill NA so downstream rcpmcp_multiple calls do not error
gm_02$result$gamma_M2_adj <- .fill_na(gm_02$result$gamma_M2_adj)

# =============================================================================
# Figure 1 data: basic properties (N, power, null RCP) vs K
# =============================================================================

.build_fig1_all <- function(ss_obj, gm_obj, fs, f1_label) {
  N_vec <- ss_obj$result$N
  Sigma <- make_sigma(0, K_MAX)

  # Power under H1 (formula + simulation, unadjusted gamma)
  pw_f <- rcpmcp_multiple(
    delta = DELTA, Sigma = Sigma, N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha = ALPHA, approach = "formula"
  )$result

  pw_s <- rcpmcp_multiple(
    delta = DELTA, Sigma = Sigma, N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha = ALPHA, approach = "simulation", nsim = NSIM, seed = SEED
  )$result[["bonferroni"]]

  # RCP0 under H0 (formula + simulation, unadjusted gamma)
  rcp0_f <- rcpmcp_multiple(
    delta = rep(0, K_MAX), Sigma = Sigma, N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha = ALPHA, approach = "formula"
  )$result

  rcp0_s <- rcpmcp_multiple(
    delta = rep(0, K_MAX), Sigma = Sigma, N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha = ALPHA, approach = "simulation", nsim = NSIM, seed = SEED
  )$result[["bonferroni"]]

  # RCP under H1: adjusted vs unadjusted (formula only)
  rcp_adj <- rcpmcp_multiple(
    delta = DELTA, Sigma = Sigma, N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = gm_obj$result$gamma_M1_adj,
    gamma_M2 = gm_obj$result$gamma_M2_adj,
    alpha = ALPHA, approach = "formula"
  )$result

  rcp_unadj <- rcpmcp_multiple(
    delta = DELTA, Sigma = Sigma, N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha = ALPHA, approach = "formula"
  )$result

  list(
    K             = seq_len(K_MAX),
    f1            = f1_label,
    N             = N_vec,
    Power_F       = pw_f$Power,
    Power_S       = pw_s$Power,
    RCP0_M1_F     = rcp0_f$RCP_M1,   RCP0_M1_S   = rcp0_s$RCP_M1,
    RCP0_M2_F     = rcp0_f$RCP_M2,   RCP0_M2_S   = rcp0_s$RCP_M2,
    RCP1_M1_Adj   = rcp_adj$RCP_M1,  RCP1_M1_Unadj = rcp_unadj$RCP_M1,
    RCP1_M2_Adj   = rcp_adj$RCP_M2,  RCP1_M2_Unadj = rcp_unadj$RCP_M2
  )
}

fig1_dat_01 <- timed("Figure 1 (f1 = 0.1)",
                     .build_fig1_all(ss_01, gm_01, FS_01, "f[1]==0.1"))
fig1_dat_02 <- timed("Figure 1 (f1 = 0.2)",
                     .build_fig1_all(ss_02, gm_02, FS_02, "f[1]==0.2"))

saveRDS(list(fig1_dat_01 = fig1_dat_01,
             fig1_dat_02 = fig1_dat_02,
             K_MAX       = K_MAX),
        file.path(data_dir, "fig1_data.rds"))

# =============================================================================
# Figure 2 data: adjusted thresholds and RCP behaviour (closed-form) vs K
# =============================================================================

.build_fig2_df <- function(gm_obj, ss_obj, fs, f1_label, na_idx = integer(0)) {
  N_vec  <- ss_obj$result$N
  Sigma  <- make_sigma(0, K_MAX)
  gm1_v  <- gm_obj$result$gamma_M1_adj
  gm2_v  <- gm_obj$result$gamma_M2_adj   # already forward-filled

  # Restore NA for gamma_M2 for honest display
  gm2_display <- gm2_v
  if (length(na_idx) > 0L) gm2_display[na_idx] <- NA_real_

  # RCP0 under H0 with adjusted gamma
  rcp0_adj <- rcpmcp_multiple(
    delta = rep(0, K_MAX), Sigma = Sigma,
    N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = gm1_v, gamma_M2 = gm2_v,
    alpha = ALPHA, approach = "formula"
  )$result

  # RCP0 under H0 with unadjusted gamma (inflation curve)
  rcp0_unadj <- rcpmcp_multiple(
    delta = rep(0, K_MAX), Sigma = Sigma,
    N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha = ALPHA, approach = "formula"
  )$result

  # RCP under H1: adjusted gamma
  rcp_adj <- rcpmcp_multiple(
    delta = DELTA, Sigma = Sigma,
    N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = gm1_v, gamma_M2 = gm2_v,
    alpha = ALPHA, approach = "formula"
  )$result

  # RCP under H1: unadjusted gamma
  rcp_unadj <- rcpmcp_multiple(
    delta = DELTA, Sigma = Sigma,
    N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha = ALPHA, approach = "formula"
  )$result

  # Mask Method 2 results at positions where gamma_M2_adj was originally NA
  if (length(na_idx) > 0L) {
    rcp0_adj$RCP_M2[na_idx]   <- NA_real_
    rcp_adj$RCP_M2[na_idx]    <- NA_real_
  }

  rbind(
    # Row 1: adjusted gamma
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = NA_character_,
               panel = "gamma[M1]~(adjusted)", value = gm1_v),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = NA_character_,
               panel = "gamma[M2]~(adjusted)", value = gm2_display),
    # Row 2: RCP_{H_0} under H0 -- adjusted vs unadjusted
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Adjusted",
               panel = "RCP[H[0]]~(Method~1)", value = rcp0_adj$RCP_M1),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Unadjusted",
               panel = "RCP[H[0]]~(Method~1)", value = rcp0_unadj$RCP_M1),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Adjusted",
               panel = "RCP[H[0]]~(Method~2)", value = rcp0_adj$RCP_M2),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Unadjusted",
               panel = "RCP[H[0]]~(Method~2)", value = rcp0_unadj$RCP_M2),
    # Row 3: RCP_{H_1} -- adjusted vs unadjusted
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Adjusted",
               panel = "RCP[H[1]]~(Method~1)", value = rcp_adj$RCP_M1),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Unadjusted",
               panel = "RCP[H[1]]~(Method~1)", value = rcp_unadj$RCP_M1),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Adjusted",
               panel = "RCP[H[1]]~(Method~2)", value = rcp_adj$RCP_M2),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = "Unadjusted",
               panel = "RCP[H[1]]~(Method~2)", value = rcp_unadj$RCP_M2)
  )
}

fig2adj_df <- timed("Figure 2 (closed-form)", rbind(
  .build_fig2_df(gm_01, ss_01, FS_01, "f[1]==0.1"),
  .build_fig2_df(gm_02, ss_02, FS_02, "f[1]==0.2", na_idx = na_idx_02)
))

saveRDS(list(fig2adj_df = fig2adj_df,
             K_MAX      = K_MAX),
        file.path(data_dir, "fig2_data.rds"))

# =============================================================================
# Figure 3 data: sensitivity to rho (K = 2..5), f1 = 0.1
# =============================================================================

.build_fig3rho_row <- function(k_val, rho_val, fs) {
  Sigma_k <- make_sigma(rho_val, k_val)

  # Sample size corresponding to this (k, rho) combination
  ss_k <- ssmcp_single(
    delta        = DELTA[seq_len(k_val)],
    Sigma        = Sigma_k,
    fs           = fs,
    K            = k_val,
    alpha        = ALPHA,
    target_power = 0.8
  )
  N_k <- ss_k$N

  # Formula
  res_f <- rcpmcp_single(
    delta    = DELTA[seq_len(k_val)],
    Sigma    = Sigma_k,
    N        = N_k,
    fs       = fs,
    K        = k_val,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha    = ALPHA, approach = "formula"
  )

  # Simulation (Bonferroni only)
  res_s <- rcpmcp_single(
    delta    = DELTA[seq_len(k_val)],
    Sigma    = Sigma_k,
    N        = N_k,
    fs       = fs,
    K        = k_val,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0,
    alpha    = ALPHA, approach = "simulation", nsim = NSIM, seed = SEED
  )
  sr <- res_s$sim_results[res_s$sim_results$mcp_method == "bonferroni", ]

  rbind(
    data.frame(K = k_val, rho = rho_val,
               panel = "Method~1", Approach = "Formula",
               value = res_f$formula_result$RCP_M1),
    data.frame(K = k_val, rho = rho_val,
               panel = "Method~1", Approach = "Simulation",
               value = sr$RCP_M1),
    data.frame(K = k_val, rho = rho_val,
               panel = "Method~2", Approach = "Formula",
               value = res_f$formula_result$RCP_M2),
    data.frame(K = k_val, rho = rho_val,
               panel = "Method~2", Approach = "Simulation",
               value = sr$RCP_M2)
  )
}

fig3rho_df <- timed("Figure 3 (rho sensitivity)", {
  fig3rho_list <- vector("list", length(2:K_MAX) * length(RHO_SEQ))
  idx <- 1L
  for (k_val in 2:K_MAX) {
    for (rho_val in RHO_SEQ) {
      fig3rho_list[[idx]] <- .build_fig3rho_row(k_val, rho_val, FS_01)
      idx <- idx + 1L
    }
  }
  do.call(rbind, fig3rho_list)
})

saveRDS(list(fig3rho_df = fig3rho_df,
             K_MAX      = K_MAX,
             RHO_SEQ    = RHO_SEQ),
        file.path(data_dir, "fig3_data.rds"))

# =============================================================================
# Figure 4 data: null RCP under four MCP procedures (simulation), vs K
# =============================================================================

.build_fig4_df <- function(gm_obj, ss_obj, fs, f1_label) {
  N_vec <- ss_obj$result$N
  gm1_v <- gm_obj$result$gamma_M1_adj
  gm2_v <- gm_obj$result$gamma_M2_adj   # forward-filled

  res_h0 <- rcpmcp_multiple(
    delta = rep(0, K_MAX), Sigma = make_sigma(0, K_MAX),
    N = N_vec, fs = fs, K_max = K_MAX,
    gamma_M1 = gm1_v, gamma_M2 = gm2_v,
    alpha = ALPHA, approach = "simulation", nsim = NSIM, seed = SEED
  )

  all_mth <- c("bonferroni", "holm", "hochberg", "hommel")
  rows <- lapply(all_mth, function(mth) {
    r0 <- res_h0$result[[mth]]
    data.frame(mcp = mth, K = r0$K, f1 = f1_label,
               RCP0_M1 = r0$RCP_M1, RCP0_M2 = r0$RCP_M2)
  })
  do.call(rbind, rows)
}

fig4_df <- timed("Figure 4 (four MCP, sim)", rbind(
  .build_fig4_df(gm_01, ss_01, FS_01, "f[1]==0.1"),
  .build_fig4_df(gm_02, ss_02, FS_02, "f[1]==0.2")
))

# Reference null RCP at K = 1 for the horizontal reference lines
fig4_ref <- list(
  M1_01 = gm_01$RCP0_M1, M2_01 = gm_01$RCP0_M2,
  M1_02 = gm_02$RCP0_M1, M2_02 = gm_02$RCP0_M2
)

saveRDS(list(fig4_df  = fig4_df,
             fig4_ref = fig4_ref,
             K_MAX    = K_MAX),
        file.path(data_dir, "fig4_data.rds"))

# =============================================================================
# Table 1 data: hypothetical MRCT application (up to K = 3)
# =============================================================================

APP_K     <- 3
APP_DELTA <- c(0.30, 0.25, 0.20)
APP_RHO   <- 0.3
APP_SIGMA <- make_sigma(APP_RHO, APP_K)
APP_FS    <- c(0.15, 0.40, 0.45)
APP_ALPHA <- 0.025

ss_app <- timed("Table 1 (sample size)", ssmcp_multiple(
  delta = APP_DELTA, Sigma = APP_SIGMA,
  fs = APP_FS, K_max = APP_K, alpha = APP_ALPHA, target_power = 0.8))

gm_app <- timed("Table 1 (adjusted gamma)", suppressWarnings(
  rcpmcp_get_gamma(
    Sigma = APP_SIGMA, N = ss_app$result$N,
    fs = APP_FS, K_max = APP_K,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0, alpha = APP_ALPHA)))
gm_app$result$gamma_M2_adj <- .fill_na(gm_app$result$gamma_M2_adj)

rcp_unadj_app <- rcpmcp_multiple(
  delta = APP_DELTA, Sigma = APP_SIGMA,
  N = ss_app$result$N, fs = APP_FS, K_max = APP_K,
  gamma_M1 = GM1_0, gamma_M2 = GM2_0,
  alpha = APP_ALPHA, approach = "formula"
)$result

rcp_adj_app <- rcpmcp_multiple(
  delta = APP_DELTA, Sigma = APP_SIGMA,
  N = ss_app$result$N, fs = APP_FS, K_max = APP_K,
  gamma_M1 = gm_app$result$gamma_M1_adj,
  gamma_M2 = gm_app$result$gamma_M2_adj,
  alpha = APP_ALPHA, approach = "formula"
)$result

rcp0_unadj_app <- rcpmcp_multiple(
  delta = rep(0, APP_K), Sigma = APP_SIGMA,
  N = ss_app$result$N, fs = APP_FS, K_max = APP_K,
  gamma_M1 = GM1_0, gamma_M2 = GM2_0,
  alpha = APP_ALPHA, approach = "formula"
)$result

rcp0_adj_app <- rcpmcp_multiple(
  delta = rep(0, APP_K), Sigma = APP_SIGMA,
  N = ss_app$result$N, fs = APP_FS, K_max = APP_K,
  gamma_M1 = gm_app$result$gamma_M1_adj,
  gamma_M2 = gm_app$result$gamma_M2_adj,
  alpha = APP_ALPHA, approach = "formula"
)$result

table1_data <- list(
  APP_K         = APP_K,
  APP_DELTA     = APP_DELTA,
  APP_RHO       = APP_RHO,
  APP_FS        = APP_FS,
  APP_ALPHA     = APP_ALPHA,
  N             = ss_app$result$N,
  gamma_M1_adj  = gm_app$result$gamma_M1_adj,
  gamma_M2_adj  = gm_app$result$gamma_M2_adj,
  null_unadj_M1 = rcp0_unadj_app$RCP_M1,
  null_unadj_M2 = rcp0_unadj_app$RCP_M2,
  null_adj_M1   = rcp0_adj_app$RCP_M1,
  null_adj_M2   = rcp0_adj_app$RCP_M2,
  unadj_M1      = rcp_unadj_app$RCP_M1,
  unadj_M2      = rcp_unadj_app$RCP_M2,
  adj_M1        = rcp_adj_app$RCP_M1,
  adj_M2        = rcp_adj_app$RCP_M2
)

saveRDS(table1_data, file.path(data_dir, "table1_data.rds"))

# =============================================================================
# Record the execution environment for reproducibility
# =============================================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_main.txt")
)

log_msg("All intermediate data written to '%s/'.", data_dir)
