# =============================================================================
# table_and_figure_manuscript.R
#
# Manuscript figures and tables for:
#   "A Cautionary Note on Regional Consistency Evaluation with Multiple
#    Primary Endpoints in Multi-Regional Clinical Trials"
#
# Prerequisites: source all five function files before running this script.
#   source("rcpmcp_single.R")
#   source("rcpmcp_multiple.R")
#   source("rcpmcp_get_gamma.R")
#   source("ssmcp_single.R")
#   source("ssmcp_multiple.R")
#
# Required packages: mvtnorm, ggplot2, rlang, scales, patchwork, ggh4x
#
# Output files
#   Figure 1 : Figure1.eps  (2 rows x 2 cols)
#     Row 1, Left : Total sample size N_K as a function of K
#     Row 1, Right: Disjunctive power xi_K under H1 as a function of K
#     Row 2, Left : Null RCP phi^(1)_K (Method 1) under H0 with unadjusted
#                   thresholds; formula (solid) vs simulation (dashed)
#     Row 2, Right: Null RCP phi^(2)_K (Method 2) under H0 with unadjusted
#                   thresholds; formula (solid) vs simulation (dashed)
#     Colour: f1 in {0.1, 0.2}; rho = 0
#
#   Figure 2 : Figure2.eps  (3 rows x 2 cols)
#     Cols   : Method 1 (left), Method 2 (right)
#     Row 1  : Adjusted thresholds gamma_{1,K} (left) and gamma_{2,K} (right)
#     Row 2  : Null RCP phi^(m)_K under H0 with adjusted thresholds
#     Row 3  : RCP phi^(m)_K under H1; adjusted (solid) vs unadjusted (dashed)
#     Content: closed-form only; rho = 0; f1 in {0.1, 0.2}
#
#   Figure 3 : Figure3.eps  (4 rows x 2 cols)
#     Rows   : K = 2, 3, 4, 5
#     Cols   : Method 1 (left), Method 2 (right)
#     x-axis : rho in {0, 0.2, 0.4, 0.6, 0.8}
#     Content: RCP phi^(m)_K under H1 with unadjusted thresholds;
#              formula (solid) vs simulation (dashed)
#     Fixed  : f1 = 0.1; N_K computed per (K, rho) at target power 80%
#
#   Figure 4 : Figure4.eps  (2 rows x 2 cols)
#     Rows   : f1 = 0.1 (top), f1 = 0.2 (bottom)
#     Cols   : Method 1 (left), Method 2 (right)
#     Content: Null RCP phi^(m)_K under H0 with adjusted thresholds for four
#              MCP procedures (Bonferroni, Holm, Hochberg, Hommel);
#              Monte Carlo simulation only; rho = 0
#     Reference lines: phi^(m)* (target null RCP at K = 1)
#
#   Table 1  : table1_application.tex
#     Application to a hypothetical MRCT with up to K = 3 primary endpoints
#     (placeholder values; to be replaced with real trial data)
#
# Global settings
#   K_max  = 5
#   delta  = rep(0.2, 5),  sigma = 1  (all endpoints identical)
#   rho    = 0 (Figures 1, 2, 4); rho in {0, 0.2, 0.4, 0.6, 0.8} (Figure 3)
#   f1     in {0.1, 0.2},  S = 3
#             FS_01 = c(0.1, 0.45, 0.45)
#             FS_02 = c(0.2, 0.40, 0.40)
#   alpha  = 0.025
#   gamma_M1 (unadjusted) = 0.5
#   gamma_M2 (unadjusted) = 0
#   target disjunctive power = 0.80
#
# Global plot settings applied to all figures
#   base_size         = 24
#   linewidth         = 1.2
#   point size        = 3
#   legend.key.width  = unit(2.5, "cm")
#   ggsave: dpi = 600, width = 16, height = 12
#   NSIM              = 1000000
# =============================================================================

library(mvtnorm)
library(ggplot2)
library(rlang)
library(scales)
library(patchwork)
library(grid)

# --------------------------------------------------------------------------- #
# Helper: uniform correlation matrix of dimension k with off-diagonal rho
# --------------------------------------------------------------------------- #
make_sigma <- function(rho, k) {
  S <- matrix(rho, k, k)
  diag(S) <- 1
  S
}

# --------------------------------------------------------------------------- #
# Helper: write a character vector to a .tex file
# --------------------------------------------------------------------------- #
write_tex <- function(lines, filename) {
  writeLines(lines, con = filename)
  message("LaTeX table saved: ", filename)
}

# --------------------------------------------------------------------------- #
# Helper: format numeric or NA as LaTeX string
# --------------------------------------------------------------------------- #
.fmt <- function(x, digits = 4, na_str = "---") {
  ifelse(is.na(x), na_str, sprintf(paste0("%.", digits, "f"), x))
}

# --------------------------------------------------------------------------- #
# Helper: forward-fill NA in a numeric vector
# --------------------------------------------------------------------------- #
.fill_na <- function(x) {
  for (i in seq_along(x)) {
    if (is.na(x[i]) && i > 1L) x[i] <- x[i - 1L]
  }
  x
}

# --------------------------------------------------------------------------- #
# Common parameters
# --------------------------------------------------------------------------- #
K_MAX     <- 5
DELTA     <- rep(0.2, K_MAX)
ALPHA     <- 0.025
GM1_0     <- 0.5     # unadjusted gamma_M1
GM2_0     <- 0.0     # unadjusted gamma_M2
NSIM      <- 1000000 # Monte Carlo iterations (reduce to 10000 for quick checks)
SEED      <- 42
BASE_SIZE <- 24      # base font size for all ggplot2 figures
LW        <- 1.2     # linewidth for all figures
PS        <- 3       # point size for all figures
LKW       <- ggplot2::unit(2.5, "cm")   # legend.key.width for all figures

FS_01 <- c(0.1, 0.45, 0.45)   # f1 = 0.1
FS_02 <- c(0.2, 0.40, 0.40)   # f1 = 0.2

# Okabe-Ito colour palette (colour-blind friendly)
OKABE <- c("#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# =============================================================================
# Pre-compute: sample sizes and adjusted thresholds (rho = 0, f1 = 0.1 and 0.2)
# Used by Figures 1, 3, 4, and Table 1
# =============================================================================

ss_01 <- ssmcp_multiple(delta = DELTA, Sigma = make_sigma(0, K_MAX),
                        fs = FS_01, K_max = K_MAX,
                        alpha = ALPHA, target_power = 0.8)

ss_02 <- ssmcp_multiple(delta = DELTA, Sigma = make_sigma(0, K_MAX),
                        fs = FS_02, K_max = K_MAX,
                        alpha = ALPHA, target_power = 0.8)

gm_01 <- suppressWarnings(
  rcpmcp_get_gamma(Sigma = make_sigma(0, K_MAX), N = ss_01$result$N,
                   fs = FS_01, K_max = K_MAX,
                   gamma_M1 = GM1_0, gamma_M2 = GM2_0, alpha = ALPHA))

gm_02 <- suppressWarnings(
  rcpmcp_get_gamma(Sigma = make_sigma(0, K_MAX), N = ss_02$result$N,
                   fs = FS_02, K_max = K_MAX,
                   gamma_M1 = GM1_0, gamma_M2 = GM2_0, alpha = ALPHA))

# Record NA positions before forward-fill (for honest plotting)
na_idx_02 <- which(is.na(gm_02$result$gamma_M2_adj))

# Forward-fill NA so downstream rcpmcp_multiple calls do not error
gm_02$result$gamma_M2_adj <- .fill_na(gm_02$result$gamma_M2_adj)


# =============================================================================
# Figure 1: Basic properties (3 rows x 2 cols) -- built with patchwork
#
# Each panel is a ggplot with facet_wrap(~label) so the strip title matches
# the style of Figures 2-4. Panels are assembled with patchwork.
#
# Row 1: N (left), Power under H1 (right)
# Row 2: RCP_{H_0} M1 (left), RCP_{H_0} M2 (right)  -- formula vs simulation
# Row 3: RCP_{H_1} M1 (left), RCP_{H_1} M2 (right)  -- adjusted vs unadjusted (formula only)
#
# Legend: two rows at bottom
#   Row A (colour x linetype/shape): f1=0.1 / f1=0.2  X  Formula / Simulation
#   Row B (colour x linetype/shape): f1=0.1 / f1=0.2  X  Adjusted / Unadjusted
#
# The colour scale is shared across all panels (f1).
# Rows 1-2 share linetype/shape = Formula/Simulation.
# Row 3 uses linetype/shape = Adjusted/Unadjusted.
# =============================================================================

Y_LIM  <- c(0.6, 1.0)
Y_BRK  <- seq(0.6, 1.0, by = 0.1)

.scale_y_rcp <- function() {
  ggplot2::scale_y_continuous(limits = Y_LIM, breaks = Y_BRK)
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("patchwork is required for Figure 1. ",
       "Install with: install.packages('patchwork')")
}

# --------------------------------------------------------------------------- #
# Compute all data
# --------------------------------------------------------------------------- #

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

fig1_dat_01 <- .build_fig1_all(ss_01, gm_01, FS_01, "f[1]==0.1")
fig1_dat_02 <- .build_fig1_all(ss_02, gm_02, FS_02, "f[1]==0.2")

# Colour values (shared across all panels)
COL_F1 <- OKABE[c(6, 5)]
names(COL_F1) <- c("f[1]==0.1", "f[1]==0.2")

# Common base theme (legend hidden; added via patchwork guide collection)
.theme_fig1 <- function(show_xlab = FALSE) {
  t <- ggplot2::theme_bw(base_size = BASE_SIZE) +
    ggplot2::theme(
      strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text         = ggplot2::element_text(face = "plain"),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "none",
      plot.margin        = ggplot2::margin(3, 5, 3, 5)
    )
  if (!show_xlab)
    t <- t + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  t
}

# --------------------------------------------------------------------------- #
# Panel builders
# --------------------------------------------------------------------------- #

# Panel with a single facet strip title and free y-axis (for N)
.p_N <- function() {
  df <- rbind(
    data.frame(K = fig1_dat_01$K, f1 = fig1_dat_01$f1,
               value = fig1_dat_01$N, label = "N~(required~sample~size)"),
    data.frame(K = fig1_dat_02$K, f1 = fig1_dat_02$f1,
               value = fig1_dat_02$N, label = "N~(required~sample~size)")
  )
  df$f1    <- factor(df$f1,    levels = c("f[1]==0.1", "f[1]==0.2"))
  df$label <- factor(df$label)
  
  ggplot2::ggplot(df, ggplot2::aes(x = K, y = value, colour = f1, group = f1)) +
    ggplot2::geom_line(linewidth = LW) +
    ggplot2::geom_point(size = PS, shape = 16) +
    ggplot2::facet_wrap(~ label, labeller = ggplot2::label_parsed) +
    ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                                breaks = seq_len(K_MAX)) +
    ggplot2::scale_y_continuous(name = NULL) +
    ggplot2::scale_colour_manual(name = NULL, values = COL_F1,
                                 labels = scales::parse_format(),
                                 guide  = "none") +
    .theme_fig1(show_xlab = TRUE)
}

# Formula vs Simulation panels (Power, RCP0_M1, RCP0_M2)
.p_fs <- function(val_F_nm, val_S_nm, label_str, show_xlab = FALSE) {
  df <- rbind(
    data.frame(K = fig1_dat_01$K, f1 = fig1_dat_01$f1,
               Approach = "Formula",    value = fig1_dat_01[[val_F_nm]],
               label = label_str),
    data.frame(K = fig1_dat_01$K, f1 = fig1_dat_01$f1,
               Approach = "Simulation", value = fig1_dat_01[[val_S_nm]],
               label = label_str),
    data.frame(K = fig1_dat_02$K, f1 = fig1_dat_02$f1,
               Approach = "Formula",    value = fig1_dat_02[[val_F_nm]],
               label = label_str),
    data.frame(K = fig1_dat_02$K, f1 = fig1_dat_02$f1,
               Approach = "Simulation", value = fig1_dat_02[[val_S_nm]],
               label = label_str)
  )
  df$f1       <- factor(df$f1,       levels = c("f[1]==0.1", "f[1]==0.2"))
  df$Approach <- factor(df$Approach, levels = c("Formula", "Simulation"))
  df$label    <- factor(df$label)
  
  ggplot2::ggplot(df, ggplot2::aes(x = K, y = value,
                                   colour   = f1,
                                   linetype = Approach,
                                   group    = interaction(f1, Approach))) +
    ggplot2::geom_line(linewidth = LW) +
    ggplot2::geom_point(size = PS, shape = 16) +
    ggplot2::facet_wrap(~ label, labeller = ggplot2::label_parsed) +
    ggplot2::scale_x_continuous(
      name   = if (show_xlab) "Number of endpoints (K)" else NULL,
      breaks = seq_len(K_MAX)
    ) +
    ggplot2::scale_y_continuous(name = NULL, limits = Y_LIM, breaks = Y_BRK) +
    ggplot2::scale_colour_manual(name = NULL, values = COL_F1,
                                 labels = scales::parse_format(),
                                 guide  = "none") +
    ggplot2::scale_linetype_manual(name = NULL,
                                   values = c("Formula"    = "solid",
                                              "Simulation" = "dashed")) +
    .theme_fig1(show_xlab = show_xlab)
}

# Adjusted vs Unadjusted panels (RCP1_M1, RCP1_M2) -- formula only
.p_adj <- function(val_Adj_nm, val_Unadj_nm, label_str, show_xlab = FALSE) {
  df <- rbind(
    data.frame(K = fig1_dat_01$K, f1 = fig1_dat_01$f1,
               Approach = "Adjusted",   value = fig1_dat_01[[val_Adj_nm]],
               label = label_str),
    data.frame(K = fig1_dat_01$K, f1 = fig1_dat_01$f1,
               Approach = "Unadjusted", value = fig1_dat_01[[val_Unadj_nm]],
               label = label_str),
    data.frame(K = fig1_dat_02$K, f1 = fig1_dat_02$f1,
               Approach = "Adjusted",   value = fig1_dat_02[[val_Adj_nm]],
               label = label_str),
    data.frame(K = fig1_dat_02$K, f1 = fig1_dat_02$f1,
               Approach = "Unadjusted", value = fig1_dat_02[[val_Unadj_nm]],
               label = label_str)
  )
  df$f1       <- factor(df$f1,       levels = c("f[1]==0.1", "f[1]==0.2"))
  df$Approach <- factor(df$Approach, levels = c("Adjusted", "Unadjusted"))
  df$label    <- factor(df$label)
  
  ggplot2::ggplot(df, ggplot2::aes(x = K, y = value,
                                   colour   = f1,
                                   linetype = Approach,
                                   shape    = Approach,
                                   group    = interaction(f1, Approach))) +
    ggplot2::geom_line(linewidth = LW) +
    ggplot2::geom_point(size = PS) +
    ggplot2::facet_wrap(~ label, labeller = ggplot2::label_parsed) +
    ggplot2::scale_x_continuous(
      name   = if (show_xlab) "Number of endpoints (K)" else NULL,
      breaks = seq_len(K_MAX)
    ) +
    ggplot2::scale_y_continuous(name = NULL, limits = Y_LIM, breaks = Y_BRK) +
    ggplot2::scale_colour_manual(name = NULL, values = COL_F1,
                                 labels = scales::parse_format()) +
    ggplot2::scale_linetype_manual(name = NULL,
                                   values = c("Adjusted"   = "solid",
                                              "Unadjusted" = "longdash")) +
    ggplot2::scale_shape_manual(name = NULL,
                                values = c("Adjusted"   = 16,
                                           "Unadjusted" = 15)) +
    .theme_fig1(show_xlab = show_xlab)
}

# --------------------------------------------------------------------------- #
# Build four panels (2 rows x 2 cols) -- all legends hidden
# --------------------------------------------------------------------------- #

p1_N    <- .p_N()
p1_Power <- .p_fs("Power_F",   "Power_S",
                  "Power~(H[1])",
                  show_xlab = TRUE)
p1_R0M1 <- .p_fs("RCP0_M1_F", "RCP0_M1_S",
                 "RCP[H[0]]~(Method~1)",
                 show_xlab = TRUE)
p1_R0M2 <- .p_fs("RCP0_M2_F", "RCP0_M2_S",
                 "RCP[H[0]]~(Method~2)",
                 show_xlab = TRUE)

# --------------------------------------------------------------------------- #
# Build legend manually as a separate ggplot
# Encode all legend entries as a dummy data frame:
#   colour entries : f1 = 0.1 (orange), f1 = 0.2 (blue)  -- solid line + point
#   linetype entries: Formula (solid black), Simulation (dashed black)
# Use a single aes mapping that produces exactly the desired legend.
# --------------------------------------------------------------------------- #

leg_df <- data.frame(
  x        = 1,
  y        = 1,
  f1       = factor(c("f[1]==0.1", "f[1]==0.2",
                      "f[1]==0.1", "f[1]==0.2"),
                    levels = c("f[1]==0.1", "f[1]==0.2")),
  Approach = factor(c("Formula", "Formula",
                      "Simulation", "Simulation"),
                    levels = c("Formula", "Simulation")),
  grp      = c("f1==0.1", "f1==0.2", "Sim_01", "Sim_02")
)

p1_leg <- ggplot2::ggplot(
  leg_df,
  ggplot2::aes(x = x, y = y,
               colour   = f1,
               linetype = Approach,
               group    = grp)
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::scale_colour_manual(
    name   = NULL,
    values = COL_F1,
    labels = scales::parse_format()
  ) +
  ggplot2::scale_linetype_manual(
    name   = NULL,
    values = c("Formula" = "solid", "Simulation" = "dashed")
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.key.width = LKW,
    legend.text      = ggplot2::element_text(size = BASE_SIZE * 0.8)
  )

# Extract legend grob
leg_grob <- ggplot2::ggplotGrob(p1_leg)
leg_idx  <- which(sapply(leg_grob$grobs, function(g) g$name) == "guide-box")
leg_only <- leg_grob$grobs[[leg_idx]]

# Wrap legend as patchwork element
p1_leg_wrap <- patchwork::wrap_elements(leg_only)

# --------------------------------------------------------------------------- #
# Assemble Figure 1 with patchwork (2 rows x 2 cols + legend row)
# --------------------------------------------------------------------------- #

Figure1 <-
  (p1_N    | p1_Power) /
  (p1_R0M1 | p1_R0M2 ) /
  p1_leg_wrap +
  patchwork::plot_layout(heights = c(1, 1, 0.15))

print(Figure1)

ggplot2::ggsave(file = paste("Figure1", "eps", sep = "."),
                plot = Figure1, dpi = 600, width = 16, height = 12)


# =============================================================================
# Figure 2: Adjusted gamma and RCP behaviour (3 rows x 2 cols)
#   Row 1: Adjusted gamma (M1 left, M2 right)
#   Row 2: RCP0 under H0 (M1 left, M2 right)  -- with adjusted gamma
#   Row 3: RCP under H1 (M1 left, M2 right)   -- with adjusted gamma
#   Same scenario as Figure 1: rho = 0, f1 in {0.1, 0.2}
#   Colour: f1; linetype/shape: single approach (formula), no sim needed
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
    # Row 2: RCP_{H_0} under H0 with adjusted gamma
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = NA_character_,
               panel = "RCP[H[0]]~(Method~1)", value = rcp0_adj$RCP_M1),
    data.frame(K = seq_len(K_MAX), f1 = f1_label, Approach = NA_character_,
               panel = "RCP[H[0]]~(Method~2)", value = rcp0_adj$RCP_M2),
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

fig2adj_df <- rbind(
  .build_fig2_df(gm_01, ss_01, FS_01, "f[1]==0.1"),
  .build_fig2_df(gm_02, ss_02, FS_02, "f[1]==0.2", na_idx = na_idx_02)
)
fig2adj_df <- fig2adj_df[!is.na(fig2adj_df$value), ]

panel_levels_fig2 <- c(
  "gamma[M1]~(adjusted)",  "gamma[M2]~(adjusted)",
  "RCP[H[0]]~(Method~1)",  "RCP[H[0]]~(Method~2)",
  "RCP[H[1]]~(Method~1)",  "RCP[H[1]]~(Method~2)"
)
fig2adj_df$panel    <- factor(fig2adj_df$panel,    levels = panel_levels_fig2)
fig2adj_df$f1       <- factor(fig2adj_df$f1,       levels = c("f[1]==0.1", "f[1]==0.2"))
fig2adj_df$Approach <- factor(fig2adj_df$Approach, levels = c("Adjusted", "Unadjusted"))

# Row grouping for facet_grid
row_map_fig2 <- c(
  "gamma[M1]~(adjusted)" = "Row1",
  "gamma[M2]~(adjusted)" = "Row1",
  "RCP[H[0]]~(Method~1)" = "Row2",
  "RCP[H[0]]~(Method~2)" = "Row2",
  "RCP[H[1]]~(Method~1)" = "Row3",
  "RCP[H[1]]~(Method~2)" = "Row3"
)
fig2adj_df$row_grp <- factor(row_map_fig2[as.character(fig2adj_df$panel)],
                             levels = c("Row1", "Row2", "Row3"))

col_map_fig2 <- c(
  "gamma[M1]~(adjusted)" = "Method 1",
  "gamma[M2]~(adjusted)" = "Method 2",
  "RCP[H[0]]~(Method~1)" = "Method 1",
  "RCP[H[0]]~(Method~2)" = "Method 2",
  "RCP[H[1]]~(Method~1)" = "Method 1",
  "RCP[H[1]]~(Method~2)" = "Method 2"
)
fig2adj_df$col_grp <- factor(col_map_fig2[as.character(fig2adj_df$panel)],
                             levels = c("Method 1", "Method 2"))

# Row 1 and Row 2 use colour = f1 only (no Approach distinction).
# Row 3 uses colour = f1 AND linetype/shape = Approach (Adjusted/Unadjusted).
# Build two subsets and overlay.

fig2adj_df_row12 <- fig2adj_df[fig2adj_df$row_grp %in% c("Row1", "Row2"), ]
fig2adj_df_row3  <- fig2adj_df[fig2adj_df$row_grp == "Row3", ]

Figure2 <- ggplot2::ggplot(
  mapping = ggplot2::aes(x = K, y = value, colour = f1)
) +
  # Row 1 and Row 2: solid line, shape = 16
  ggplot2::geom_line(
    data      = fig2adj_df_row12,
    mapping   = ggplot2::aes(group = f1),
    linewidth = LW
  ) +
  ggplot2::geom_point(
    data    = fig2adj_df_row12,
    mapping = ggplot2::aes(group = f1),
    size    = PS, shape = 16
  ) +
  # Row 3: linetype = Approach, shape = 16
  ggplot2::geom_line(
    data      = fig2adj_df_row3,
    mapping   = ggplot2::aes(linetype = Approach,
                             group    = interaction(f1, Approach)),
    linewidth = LW
  ) +
  ggplot2::geom_point(
    data    = fig2adj_df_row3,
    mapping = ggplot2::aes(group = interaction(f1, Approach)),
    size    = PS, shape = 16
  ) +
  ggplot2::facet_grid(
    row_grp ~ col_grp,
    scales   = "free_y",
    labeller = ggplot2::labeller(
      row_grp = ggplot2::as_labeller(
        c(Row1 = "Adjusted~threshold",
          Row2 = "RCP[H[0]]",
          Row3 = "RCP[H[1]]"),
        default = ggplot2::label_parsed
      ),
      col_grp = ggplot2::label_value
    )
  ) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(name = NULL) +
  ggplot2::scale_colour_manual(
    name   = NULL,
    values = OKABE[c(6, 5)],
    labels = scales::parse_format()
  ) +
  ggplot2::scale_linetype_manual(
    name   = NULL,
    values = c("Adjusted" = "solid", "Unadjusted" = "longdash"),
    na.translate = FALSE
  ) +
  ggplot2::theme_bw(base_size = BASE_SIZE) +
  ggplot2::theme(
    strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
    strip.text         = ggplot2::element_text(face = "plain"),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position    = "bottom",
    legend.key.width   = LKW,
    plot.margin        = ggplot2::margin(8, 8, 8, 8)
  )

# Apply y-axis limits via ggh4x:
#   Row 1 (gamma panels): free y-axis
#   Row 2 and Row 3 (RCP panels): fixed to [0.6, 1.0] by 0.1
# Panel order in facet_grid(row_grp ~ col_grp): row-major
#   (Row1,M1), (Row1,M2), (Row2,M1), (Row2,M2), (Row3,M1), (Row3,M2)
if (!requireNamespace("ggh4x", quietly = TRUE)) {
  stop("ggh4x is required for Figure 2 y-axis control. ",
       "Install with: install.packages('ggh4x')")
}

Figure2 <- Figure2 +
  ggh4x::facetted_pos_scales(y = list(
    NULL,            # (Row1, Method 1): gamma -- free
    NULL,            # (Row1, Method 2): gamma -- free
    .scale_y_rcp(),  # (Row2, Method 1): RCP_{H_0}
    .scale_y_rcp(),  # (Row2, Method 2): RCP_{H_0}
    .scale_y_rcp(),  # (Row3, Method 1): RCP_{H_1}
    .scale_y_rcp()   # (Row3, Method 2): RCP_{H_1}
  ))

print(Figure2)

ggplot2::ggsave(file = paste("Figure2", "eps", sep = "."),
                plot = Figure2, dpi = 600, width = 16, height = 12)


# =============================================================================
# Figure 3: Sensitivity to rho (4 rows x 2 cols)
#   Rows: k = 2, 3, 4, 5
#   Cols: Method 1 (left), Method 2 (right)
#   x-axis: rho in {0, 0.2, 0.4, 0.6, 0.8}
#   f1 = 0.1 fixed; formula (solid) vs simulation (dashed)
#   N used: N_k computed for each (k, rho) via ssmcp_single (target power 80%)
# =============================================================================

RHO_SEQ <- seq(0, 0.8, by = 0.2)

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

fig3rho_list <- vector("list", length(2:K_MAX) * length(RHO_SEQ))
idx <- 1L
for (k_val in 2:K_MAX) {
  for (rho_val in RHO_SEQ) {
    fig3rho_list[[idx]] <- .build_fig3rho_row(k_val, rho_val, FS_01)
    idx <- idx + 1L
  }
}
fig3rho_df <- do.call(rbind, fig3rho_list)

fig3rho_df$k_label  <- factor(paste0("K == ", fig3rho_df$K),
                              levels = paste0("K == ", 2:K_MAX))
fig3rho_df$panel    <- factor(fig3rho_df$panel,
                              levels = c("Method~1", "Method~2"))
fig3rho_df$Approach <- factor(fig3rho_df$Approach,
                              levels = c("Formula", "Simulation"))

Figure3 <- ggplot2::ggplot(
  fig3rho_df,
  ggplot2::aes(x = rho, y = value,
               colour   = Approach,
               linetype = Approach,
               group    = Approach)
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::facet_grid(
    k_label ~ panel,
    scales   = "free_y",
    labeller = ggplot2::labeller(
      k_label = ggplot2::label_parsed,
      panel   = ggplot2::label_parsed
    )
  ) +
  ggplot2::scale_x_continuous(
    name   = expression(paste("Correlation coefficient (", rho, ")")),
    breaks = RHO_SEQ
  ) +
  ggplot2::scale_y_continuous(name = "RCP", limits = Y_LIM, breaks = Y_BRK) +
  ggplot2::scale_colour_manual(
    name   = NULL,
    values = c("Formula" = OKABE[5], "Simulation" = OKABE[5])
  ) +
  ggplot2::scale_linetype_manual(
    name   = NULL,
    values = c("Formula" = "solid", "Simulation" = "dashed")
  ) +
  ggplot2::theme_bw(base_size = BASE_SIZE) +
  ggplot2::theme(
    strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
    strip.text         = ggplot2::element_text(face = "plain"),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position    = "bottom",
    legend.key.width   = LKW,
    plot.margin        = ggplot2::margin(8, 8, 8, 8)
  )

print(Figure3)

ggplot2::ggsave(file = paste("Figure3", "eps", sep = "."),
                plot = Figure3, dpi = 600, width = 16, height = 12)


# =============================================================================
# Figure 4: Null RCP under four MCP methods (rho = 0 only)
#   Replicates the "conservativeness" check for all four procedures
#   Layout: free (single facet_grid: panel x f1)
#   Colour: MCP method; f1 shown as separate row in facet
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

fig4_df <- rbind(
  .build_fig4_df(gm_01, ss_01, FS_01, "f[1]==0.1"),
  .build_fig4_df(gm_02, ss_02, FS_02, "f[1]==0.2")
)

mth_labels <- c(bonferroni = "Bonferroni", holm = "Holm",
                hochberg = "Hochberg", hommel = "Hommel")

fig4_long <- rbind(
  data.frame(K = fig4_df$K, mcp = fig4_df$mcp, f1 = fig4_df$f1,
             panel = "RCP[H[0]]~(Method~1)", value = fig4_df$RCP0_M1),
  data.frame(K = fig4_df$K, mcp = fig4_df$mcp, f1 = fig4_df$f1,
             panel = "RCP[H[0]]~(Method~2)", value = fig4_df$RCP0_M2)
)
fig4_long$mcp_label <- factor(mth_labels[fig4_long$mcp],
                              levels = c("Bonferroni", "Holm",
                                         "Hochberg", "Hommel"))
fig4_long$panel     <- factor(fig4_long$panel,
                              levels = c("RCP[H[0]]~(Method~1)",
                                         "RCP[H[0]]~(Method~2)"))
fig4_long$f1        <- factor(fig4_long$f1,
                              levels = c("f[1]==0.1", "f[1]==0.2"))

# Reference lines: target RCP0 for each f1 separately
# gm_01$RCP0_M1/M2: reference for f1 = 0.1
# gm_02$RCP0_M1/M2: reference for f1 = 0.2
ref_fig4 <- rbind(
  data.frame(
    f1    = factor("f[1]==0.1", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~1)", levels = levels(fig4_long$panel)),
    yint  = gm_01$RCP0_M1
  ),
  data.frame(
    f1    = factor("f[1]==0.1", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~2)", levels = levels(fig4_long$panel)),
    yint  = gm_01$RCP0_M2
  ),
  data.frame(
    f1    = factor("f[1]==0.2", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~1)", levels = levels(fig4_long$panel)),
    yint  = gm_02$RCP0_M1
  ),
  data.frame(
    f1    = factor("f[1]==0.2", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~2)", levels = levels(fig4_long$panel)),
    yint  = gm_02$RCP0_M2
  )
)

Figure4 <- ggplot2::ggplot(
  fig4_long,
  ggplot2::aes(x = K, y = value,
               colour = mcp_label, shape = mcp_label, group = mcp_label)
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::geom_hline(data = ref_fig4,
                      ggplot2::aes(yintercept = yint),
                      linetype = "dashed", colour = "grey40",
                      linewidth = 0.8, inherit.aes = FALSE) +
  ggplot2::facet_grid(
    f1 ~ panel,
    scales   = "free_y",
    labeller = ggplot2::labeller(
      panel = ggplot2::label_parsed,
      f1    = ggplot2::label_parsed
    )
  ) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(
    name   = expression(RCP[H[0]]),
    limits = Y_LIM, breaks = Y_BRK
  ) +
  ggplot2::scale_colour_manual(name = NULL,
                               values = OKABE[c(6, 1, 5, 3)]) +
  ggplot2::theme_bw(base_size = BASE_SIZE) +
  ggplot2::theme(
    strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
    strip.text         = ggplot2::element_text(face = "plain"),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position    = "bottom",
    legend.key.width   = LKW,
    plot.margin        = ggplot2::margin(8, 8, 8, 8)
  )

print(Figure4)

ggplot2::ggsave(file = paste("Figure4", "eps", sep = "."),
                plot = Figure4, dpi = 600, width = 16, height = 12)


# =============================================================================
# Table 1: Application -- hypothetical MRCT (placeholder values)
#   k = 1, 2, 3 co-primary endpoints
#   Shows: N_k, alpha/k, unadjusted RCP, adjusted threshold, adjusted RCP
# =============================================================================

APP_K     <- 3
APP_DELTA <- c(0.30, 0.25, 0.20)
APP_RHO   <- 0.3
APP_SIGMA <- make_sigma(APP_RHO, APP_K)
APP_FS    <- c(0.15, 0.40, 0.45)
APP_ALPHA <- 0.025

ss_app <- ssmcp_multiple(
  delta = APP_DELTA, Sigma = APP_SIGMA,
  fs = APP_FS, K_max = APP_K,
  alpha = APP_ALPHA, target_power = 0.8
)
gm_app <- suppressWarnings(
  rcpmcp_get_gamma(
    Sigma = APP_SIGMA, N = ss_app$result$N,
    fs = APP_FS, K_max = APP_K,
    gamma_M1 = GM1_0, gamma_M2 = GM2_0, alpha = APP_ALPHA
  )
)
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

tab1_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  paste0("\\caption{Application to a hypothetical multi-regional clinical ",
         "trial with up to $K = 3$ co-primary endpoints. ",
         "True treatment effects: $\\boldsymbol{\\delta} = (0.30, 0.25, 0.20)$, ",
         "$\\sigma_{k} = 1$ ($k = 1, 2, 3$), endpoint correlation $\\rho = 0.3$, ",
         "$f_{1} = 0.15$ (region of interest), $S = 3$, ",
         "$\\alpha = 0.025$, target disjunctive power $= 0.80$. ",
         "Unadjusted thresholds: $\\gamma_{{\\rm M}1} = 0.5$, $\\gamma_{{\\rm M}2} = 0$. ",
         "Adjusted thresholds $\\gamma_{{\\rm M}1,K}$ and $\\gamma_{{\\rm M}2,K}$ are proposed herein.}"),
  "\\label{tab:application}",
  "\\begin{tabular}{r r r cc cc cc}",
  "\\hline",
  paste(" & & &",
        "\\multicolumn{2}{c}{Unadjusted $\\phi^{(m)}_{K}$} &",
        "\\multicolumn{2}{c}{Adjusted} &",
        "\\multicolumn{2}{c}{Adjusted $\\phi^{(m)}_{K}$} \\\\"),
  paste(" & & &",
        "\\multicolumn{2}{c}{under ${\\rm H}_{1}$} &",
        "\\multicolumn{2}{c}{threshold} &",
        "\\multicolumn{2}{c}{under ${\\rm H}_{1}$} \\\\"),
  "\\cline{4-5}\\cline{6-7}\\cline{8-9}",
  paste("$K$ & $N_{K}$ & $\\alpha/K$ &",
        "$\\phi^{(1)}_{K}$ & $\\phi^{(2)}_{K}$ &",
        "$\\gamma_{{\\rm M}1,K}$ & $\\gamma_{{\\rm M}2,K}$ &",
        "$\\phi^{(1)}_{K}$ & $\\phi^{(2)}_{K}$ \\\\"),
  "\\hline"
)

for (k in seq_len(APP_K)) {
  tab1_lines <- c(tab1_lines, sprintf(
    "%d & %d & %.5f & %.4f & %.4f & %.4f & %s & %.4f & %.4f \\\\",
    k,
    ss_app$result$N[k], APP_ALPHA / k,
    rcp_unadj_app$RCP_M1[k], rcp_unadj_app$RCP_M2[k],
    gm_app$result$gamma_M1_adj[k],
    .fmt(gm_app$result$gamma_M2_adj[k]),
    rcp_adj_app$RCP_M1[k], rcp_adj_app$RCP_M2[k]
  ))
}
tab1_lines <- c(tab1_lines, "\\hline", "\\end{tabular}", "\\end{table}")

write_tex(tab1_lines, "table1_application.tex")