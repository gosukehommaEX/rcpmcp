# =============================================================================
# data_generate_supplement.R
#
# Computation stage for the Supporting Information figures of:
#   "A Cautionary Note on Regional Consistency Evaluation with Multiple
#    Primary Endpoints in Multi-Regional Clinical Trials"
#
# This script performs all heavy numerical computation for the sensitivity
# analyses requested by the editor and reviewers and writes the resulting
# intermediate data to "data/*.rds". The companion rendering script
# table_and_figure_supplement.R reads these files and produces the figures.
#
# Scenarios
#   S1  Unknown variance (Wishart-sampled, t-statistic) vs known variance.
#   S2  Heterogeneous effect in the region of interest (region 1).
#   S3  Mixed endpoints (a subset of endpoints has no treatment effect).
#   S4  Correlation combined with a null region of interest.
#   S5  Unequal endpoint variances.
#   S6  Non-standard (AR(1)) endpoint correlation structure.
#
# Usage
#   Set the working directory to this script's location (inst/sup_info), then
#   run the whole file. Intermediate results are written to "data/" relative
#   to the working directory. The package must be installed or loaded first:
#     library(rcpmcp)              # installed package, or
#     devtools::load_all()         # from the package root during development
#
# Outputs (relative to the working directory)
#   data/siS1_data.rds ... data/siS6_data.rds
#   data/data_generate_supplement.log
#   data/sessionInfo_data_generate_supplement.txt
# =============================================================================

library(rcpmcp)

# --------------------------------------------------------------------------- #
# Output location (relative paths only)
# --------------------------------------------------------------------------- #
data_dir <- "data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_supplement.log")
cat(sprintf("data_generate_supplement.R log -- %s\n", format(Sys.time())),
    file = log_file, append = FALSE)

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  message(msg)
}

timed <- function(label, expr) {
  t <- system.time(val <- force(expr))["elapsed"]
  log_msg("%-34s %7.1f sec", label, t)
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

.ar1_sigma <- function(rho, k) {
  outer(seq_len(k), seq_len(k), function(i, j) rho^abs(i - j))
}

.fill_na <- function(x) {
  for (i in seq_along(x)) {
    if (is.na(x[i]) && i > 1L) x[i] <- x[i - 1L]
  }
  x
}

# --------------------------------------------------------------------------- #
# Common parameters (shared with the main-text baseline, f1 = 0.1)
# --------------------------------------------------------------------------- #
K_MAX     <- 5
DELTA_VAL <- 0.2            # treatment effect where present
ALPHA     <- 0.025
GM1_0     <- 0.5            # unadjusted gamma_M1
GM2_0     <- 0.0            # unadjusted gamma_M2
NSIM      <- 1000000        # Monte Carlo iterations (reduce to 10000 for checks)
SEED      <- 42
FS        <- c(0.1, 0.45, 0.45)
S_REG     <- length(FS)     # number of regions

log_msg("Settings: K_MAX = %d, NSIM = %d, SEED = %d, alpha = %.3f",
        K_MAX, NSIM, SEED, ALPHA)

# Determine the trial sample size (homogeneous alternative) and the adjusted
# thresholds for a given endpoint covariance matrix.
.si_setup <- function(Sigma) {
  ss <- ssmcp_multiple(
    delta = rep(DELTA_VAL, K_MAX), Sigma = Sigma,
    fs = FS, K_max = K_MAX, alpha = ALPHA, target_power = 0.8
  )
  N_vec <- ss$result$N
  gm <- suppressWarnings(rcpmcp_get_gamma(
    Sigma = Sigma, N = N_vec, fs = FS, K_max = K_MAX,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0, alpha = ALPHA
  ))
  na_idx <- which(is.na(gm$result$gamma_M2_adj))
  list(
    N      = N_vec,
    gm1    = gm$result$gamma_M1_adj,
    gm2    = .fill_na(gm$result$gamma_M2_adj),
    na_idx = na_idx
  )
}

# RCP via simulation (Bonferroni), returning a data frame with K, RCP_M1, RCP_M2
.rcp_sim <- function(delta, Sigma, N, gm1, gm2, variance_known = TRUE) {
  rcpmcp_multiple(
    delta = delta, Sigma = Sigma, N = N, fs = FS, K_max = K_MAX,
    gamma_M1 = gm1, gamma_M2 = gm2, alpha = ALPHA,
    approach = "simulation", nsim = NSIM, seed = SEED,
    variance_known = variance_known
  )$result[["bonferroni"]]
}

# RCP via closed-form (Bonferroni)
.rcp_formula <- function(delta, Sigma, N, gm1, gm2) {
  rcpmcp_multiple(
    delta = delta, Sigma = Sigma, N = N, fs = FS, K_max = K_MAX,
    gamma_M1 = gm1, gamma_M2 = gm2, alpha = ALPHA, approach = "formula"
  )$result
}

# Region-heterogeneous delta matrix (K_MAX-by-S): region 1 = r1_val, others base
.delta_region <- function(r1_val, base_val = DELTA_VAL) {
  M <- matrix(base_val, nrow = K_MAX, ncol = S_REG)
  M[, 1] <- r1_val
  M
}

# Baseline setup (independence, equal variances) reused by S1, S2, S3
base <- timed("baseline setup (Sigma = I)", .si_setup(diag(K_MAX)))

# =============================================================================
# S1: Unknown variance (Wishart, t-statistic) vs known variance
# =============================================================================

siS1_df <- timed("S1 unknown variance", {
  d_alt  <- rep(DELTA_VAL, K_MAX)
  d_null <- rep(0, K_MAX)
  Sig    <- diag(K_MAX)

  h1_known   <- .rcp_sim(d_alt,  Sig, base$N, base$gm1, base$gm2, variance_known = TRUE)
  h1_unknown <- .rcp_sim(d_alt,  Sig, base$N, base$gm1, base$gm2, variance_known = FALSE)
  h0_known   <- .rcp_sim(d_null, Sig, base$N, base$gm1, base$gm2, variance_known = TRUE)
  h0_unknown <- .rcp_sim(d_null, Sig, base$N, base$gm1, base$gm2, variance_known = FALSE)

  mk <- function(res, H, variance, panel, value) {
    data.frame(K = res$K, H = H, variance = variance, panel = panel, value = value)
  }
  rbind(
    mk(h0_known,   "H0", "Known",   "Method~1", h0_known$RCP_M1),
    mk(h0_known,   "H0", "Known",   "Method~2", h0_known$RCP_M2),
    mk(h0_unknown, "H0", "Unknown", "Method~1", h0_unknown$RCP_M1),
    mk(h0_unknown, "H0", "Unknown", "Method~2", h0_unknown$RCP_M2),
    mk(h1_known,   "H1", "Known",   "Method~1", h1_known$RCP_M1),
    mk(h1_known,   "H1", "Known",   "Method~2", h1_known$RCP_M2),
    mk(h1_unknown, "H1", "Unknown", "Method~1", h1_unknown$RCP_M1),
    mk(h1_unknown, "H1", "Unknown", "Method~2", h1_unknown$RCP_M2)
  )
})

saveRDS(list(siS1_df = siS1_df, K_MAX = K_MAX),
        file.path(data_dir, "siS1_data.rds"))

# =============================================================================
# S2: Heterogeneous effect in the region of interest (region 1)
#     RCP under H1 with adjusted thresholds, for region-1 effect in {0, 0.1, 0.2}
# =============================================================================

siS2_df <- timed("S2 region-of-interest heterogeneity", {
  r1_vals <- c(0.0, 0.1, 0.2)
  Sig     <- diag(K_MAX)

  rows <- lapply(r1_vals, function(r1) {
    res <- .rcp_sim(.delta_region(r1), Sig, base$N, base$gm1, base$gm2)
    rbind(
      data.frame(K = res$K, r1 = r1, panel = "Method~1", value = res$RCP_M1),
      data.frame(K = res$K, r1 = r1, panel = "Method~2", value = res$RCP_M2)
    )
  })
  do.call(rbind, rows)
})

saveRDS(list(siS2_df = siS2_df, K_MAX = K_MAX),
        file.path(data_dir, "siS2_data.rds"))

# =============================================================================
# S3: Mixed endpoints (odd-indexed endpoints have an effect, even-indexed none)
#     RCP under H1 with adjusted thresholds; homogeneous vs mixed
# =============================================================================

siS3_df <- timed("S3 mixed endpoints", {
  Sig      <- diag(K_MAX)
  d_homog  <- rep(DELTA_VAL, K_MAX)
  d_mixed  <- rep(0, K_MAX)
  d_mixed[seq(1, K_MAX, by = 2)] <- DELTA_VAL   # endpoints 1, 3, 5 active

  configs <- list(Homogeneous = d_homog, Mixed = d_mixed)

  rows <- list()
  for (cfg in names(configs)) {
    d  <- configs[[cfg]]
    rf <- .rcp_formula(d, Sig, base$N, base$gm1, base$gm2)
    rs <- .rcp_sim(d, Sig, base$N, base$gm1, base$gm2)
    rows[[length(rows) + 1L]] <- rbind(
      data.frame(K = rf$K, config = cfg, Approach = "Formula",
                 panel = "Method~1", value = rf$RCP_M1),
      data.frame(K = rf$K, config = cfg, Approach = "Formula",
                 panel = "Method~2", value = rf$RCP_M2),
      data.frame(K = rs$K, config = cfg, Approach = "Simulation",
                 panel = "Method~1", value = rs$RCP_M1),
      data.frame(K = rs$K, config = cfg, Approach = "Simulation",
                 panel = "Method~2", value = rs$RCP_M2)
    )
  }
  do.call(rbind, rows)
})

saveRDS(list(siS3_df = siS3_df, K_MAX = K_MAX),
        file.path(data_dir, "siS3_data.rds"))

# =============================================================================
# S4: Correlation combined with a null region of interest (region 1 = 0)
#     RCP under H1 with adjusted thresholds, for rho in {0, 0.3, 0.5, 0.8}
# =============================================================================

siS4_df <- timed("S4 correlation x null region", {
  rho_vals <- c(0.0, 0.3, 0.5, 0.8)
  d_region <- .delta_region(0.0)   # region 1 null, others = 0.2

  rows <- lapply(rho_vals, function(rho) {
    Sig <- make_sigma(rho, K_MAX)
    st  <- .si_setup(Sig)
    res <- .rcp_sim(d_region, Sig, st$N, st$gm1, st$gm2)
    rbind(
      data.frame(K = res$K, rho = rho, panel = "Method~1", value = res$RCP_M1),
      data.frame(K = res$K, rho = rho, panel = "Method~2", value = res$RCP_M2)
    )
  })
  do.call(rbind, rows)
})

saveRDS(list(siS4_df = siS4_df, K_MAX = K_MAX, RHO_VALS = c(0.0, 0.3, 0.5, 0.8)),
        file.path(data_dir, "siS4_data.rds"))

# =============================================================================
# S5: Unequal endpoint variances (sigma_k linearly spaced in [0.8, 1.2])
#     RCP under H0 and H1 with adjusted thresholds; equal vs unequal variances
# =============================================================================

siS5_df <- timed("S5 unequal variances", {
  sigma_vec  <- seq(0.8, 1.2, length.out = K_MAX)
  Sig_uneq   <- diag(sigma_vec^2)
  st_uneq    <- .si_setup(Sig_uneq)

  d_alt  <- rep(DELTA_VAL, K_MAX)
  d_null <- rep(0, K_MAX)

  build <- function(struct, Sig, st) {
    h0f <- .rcp_formula(d_null, Sig, st$N, st$gm1, st$gm2)
    h1f <- .rcp_formula(d_alt,  Sig, st$N, st$gm1, st$gm2)
    h0s <- .rcp_sim(d_null, Sig, st$N, st$gm1, st$gm2)
    h1s <- .rcp_sim(d_alt,  Sig, st$N, st$gm1, st$gm2)
    rbind(
      data.frame(K = h0f$K, struct = struct, H = "H0", Approach = "Formula",
                 panel = "Method~1", value = h0f$RCP_M1),
      data.frame(K = h0f$K, struct = struct, H = "H0", Approach = "Formula",
                 panel = "Method~2", value = h0f$RCP_M2),
      data.frame(K = h1f$K, struct = struct, H = "H1", Approach = "Formula",
                 panel = "Method~1", value = h1f$RCP_M1),
      data.frame(K = h1f$K, struct = struct, H = "H1", Approach = "Formula",
                 panel = "Method~2", value = h1f$RCP_M2),
      data.frame(K = h0s$K, struct = struct, H = "H0", Approach = "Simulation",
                 panel = "Method~1", value = h0s$RCP_M1),
      data.frame(K = h0s$K, struct = struct, H = "H0", Approach = "Simulation",
                 panel = "Method~2", value = h0s$RCP_M2),
      data.frame(K = h1s$K, struct = struct, H = "H1", Approach = "Simulation",
                 panel = "Method~1", value = h1s$RCP_M1),
      data.frame(K = h1s$K, struct = struct, H = "H1", Approach = "Simulation",
                 panel = "Method~2", value = h1s$RCP_M2)
    )
  }

  rbind(
    build("Equal",   diag(K_MAX), base),
    build("Unequal", Sig_uneq,    st_uneq)
  )
})

saveRDS(list(siS5_df = siS5_df, K_MAX = K_MAX),
        file.path(data_dir, "siS5_data.rds"))

# =============================================================================
# S6: Non-standard correlation structure
#     Independence vs compound symmetry (rho = 0.5) vs AR(1) (rho = 0.5)
#     RCP under H0 and H1 with adjusted thresholds
# =============================================================================

siS6_df <- timed("S6 non-standard correlation", {
  RHO_S6  <- 0.5
  d_alt   <- rep(DELTA_VAL, K_MAX)
  d_null  <- rep(0, K_MAX)

  structs <- list(
    Independence = diag(K_MAX),
    CS           = make_sigma(RHO_S6, K_MAX),
    AR1          = .ar1_sigma(RHO_S6, K_MAX)
  )

  build <- function(struct, Sig, st) {
    h0f <- .rcp_formula(d_null, Sig, st$N, st$gm1, st$gm2)
    h1f <- .rcp_formula(d_alt,  Sig, st$N, st$gm1, st$gm2)
    h0s <- .rcp_sim(d_null, Sig, st$N, st$gm1, st$gm2)
    h1s <- .rcp_sim(d_alt,  Sig, st$N, st$gm1, st$gm2)
    rbind(
      data.frame(K = h0f$K, struct = struct, H = "H0", Approach = "Formula",
                 panel = "Method~1", value = h0f$RCP_M1),
      data.frame(K = h0f$K, struct = struct, H = "H0", Approach = "Formula",
                 panel = "Method~2", value = h0f$RCP_M2),
      data.frame(K = h1f$K, struct = struct, H = "H1", Approach = "Formula",
                 panel = "Method~1", value = h1f$RCP_M1),
      data.frame(K = h1f$K, struct = struct, H = "H1", Approach = "Formula",
                 panel = "Method~2", value = h1f$RCP_M2),
      data.frame(K = h0s$K, struct = struct, H = "H0", Approach = "Simulation",
                 panel = "Method~1", value = h0s$RCP_M1),
      data.frame(K = h0s$K, struct = struct, H = "H0", Approach = "Simulation",
                 panel = "Method~2", value = h0s$RCP_M2),
      data.frame(K = h1s$K, struct = struct, H = "H1", Approach = "Simulation",
                 panel = "Method~1", value = h1s$RCP_M1),
      data.frame(K = h1s$K, struct = struct, H = "H1", Approach = "Simulation",
                 panel = "Method~2", value = h1s$RCP_M2)
    )
  }

  st_cs  <- .si_setup(structs$CS)
  st_ar1 <- .si_setup(structs$AR1)

  rbind(
    build("Independence", structs$Independence, base),
    build("CS",           structs$CS,           st_cs),
    build("AR1",          structs$AR1,          st_ar1)
  )
})

saveRDS(list(siS6_df = siS6_df, K_MAX = K_MAX),
        file.path(data_dir, "siS6_data.rds"))

# =============================================================================
# Record the execution environment for reproducibility
# =============================================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_supplement.txt")
)

log_msg("All Supporting Information data written to '%s/'.", data_dir)
