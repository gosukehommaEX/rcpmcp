# =============================================================================
# table_and_figure_manuscript.R
#
# Rendering stage for the manuscript figures and tables of:
#   "A Cautionary Note on Regional Consistency Evaluation with Multiple
#    Primary Endpoints in Multi-Regional Clinical Trials"
#
# This script reads the intermediate data written by data_generate_main.R
# (data/*.rds) and produces the figures and the LaTeX table. It performs no
# numerical computation, so figures and tables can be adjusted without
# re-running any simulation.
#
# Prerequisite
#   Run data_generate_main.R first (from this same working directory) so that
#   data/fig1_data.rds, data/fig2_data.rds, data/fig3_data.rds,
#   data/fig4_data.rds, and data/table1_data.rds exist.
#
# Required packages: ggplot2, rlang, scales, patchwork, grid, ggh4x
#
# Outputs (relative to the working directory)
#   results/Figure1.eps, results/Figure1.pdf
#   results/Figure2.eps, results/Figure2.pdf
#   results/Figure3.eps, results/Figure3.pdf
#   results/Figure4.eps, results/Figure4.pdf
#   results/table1_application.tex
# =============================================================================

library(ggplot2)
library(rlang)
library(scales)
library(patchwork)
library(grid)

# --------------------------------------------------------------------------- #
# Read intermediate data produced by data_generate_main.R
# --------------------------------------------------------------------------- #
data_dir <- "data"

fig1_data   <- readRDS(file.path(data_dir, "fig1_data.rds"))
fig2_data   <- readRDS(file.path(data_dir, "fig2_data.rds"))
fig3_data   <- readRDS(file.path(data_dir, "fig3_data.rds"))
fig4_data   <- readRDS(file.path(data_dir, "fig4_data.rds"))
table1_data <- readRDS(file.path(data_dir, "table1_data.rds"))

fig1_dat_01 <- fig1_data$fig1_dat_01
fig1_dat_02 <- fig1_data$fig1_dat_02
fig2adj_df  <- fig2_data$fig2adj_df
fig3rho_df  <- fig3_data$fig3rho_df
fig4_df     <- fig4_data$fig4_df
fig4_ref_v  <- fig4_data$fig4_ref
K_MAX       <- fig1_data$K_MAX
RHO_SEQ     <- fig3_data$RHO_SEQ

# Output location for figures and tables (relative path, separate from data/)
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --------------------------------------------------------------------------- #
# Plotting constants
# --------------------------------------------------------------------------- #
BASE_SIZE <- 24      # base font size for all ggplot2 figures
LW        <- 1.2     # linewidth for all figures
PS        <- 3       # point size for all figures
LKW       <- ggplot2::unit(2.5, "cm")   # legend.key.width for all figures

# Okabe-Ito colour palette (colour-blind friendly)
OKABE <- c("#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

Y_LIM  <- c(0.6, 1.0)
Y_BRK  <- seq(0.6, 1.0, by = 0.1)

# Colour values (shared across all panels)
COL_F1 <- OKABE[c(6, 5)]
names(COL_F1) <- c("f[1]==0.1", "f[1]==0.2")

# --------------------------------------------------------------------------- #
# Rendering helpers
# --------------------------------------------------------------------------- #
.fmt <- function(x, digits = 4, na_str = "---") {
  ifelse(is.na(x), na_str, sprintf(paste0("%.", digits, "f"), x))
}

write_tex <- function(lines, filename) {
  writeLines(lines, con = filename)
  message("LaTeX table saved: ", filename)
}

.scale_y_rcp <- function() {
  ggplot2::scale_y_continuous(limits = Y_LIM, breaks = Y_BRK)
}

# Save a figure to the results directory in both EPS and PDF formats
save_figure <- function(plot, name, width = 16, height = 12, dpi = 600) {
  for (ext in c("eps", "pdf")) {
    ggplot2::ggsave(
      filename = file.path(results_dir, paste(name, ext, sep = ".")),
      plot     = plot,
      device   = ext,
      dpi      = dpi,
      width    = width,
      height   = height
    )
  }
}

# =============================================================================
# Figure 1: Basic properties (2 rows x 2 cols) -- built with patchwork
# =============================================================================

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("patchwork is required for Figure 1. ",
       "Install with: install.packages('patchwork')")
}

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

# Build four panels (2 rows x 2 cols) -- all legends hidden
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

# Build legend manually as a separate ggplot
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

# Assemble Figure 1 with patchwork (2 rows x 2 cols + legend row)
Figure1 <-
  (p1_N    | p1_Power) /
  (p1_R0M1 | p1_R0M2 ) /
  p1_leg_wrap +
  patchwork::plot_layout(heights = c(1, 1, 0.15))

print(Figure1)

save_figure(Figure1, "Figure1")

# =============================================================================
# Figure 2: Adjusted gamma and RCP behaviour (3 rows x 2 cols)
# =============================================================================

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

# Row 1 (adjusted thresholds) uses colour = f1 only (single line per f1).
# Rows 2 and 3 (RCP) use colour = f1 AND linetype = Approach (Adjusted/Unadjusted).
fig2adj_df_row1  <- fig2adj_df[fig2adj_df$row_grp == "Row1", ]
fig2adj_df_row23 <- fig2adj_df[fig2adj_df$row_grp %in% c("Row2", "Row3"), ]

Figure2 <- ggplot2::ggplot(
  mapping = ggplot2::aes(x = K, y = value, colour = f1)
) +
  # Row 1: solid line, shape = 16
  ggplot2::geom_line(
    data      = fig2adj_df_row1,
    mapping   = ggplot2::aes(group = f1),
    linewidth = LW
  ) +
  ggplot2::geom_point(
    data    = fig2adj_df_row1,
    mapping = ggplot2::aes(group = f1),
    size    = PS, shape = 16
  ) +
  # Rows 2 and 3: linetype = Approach, shape = 16
  ggplot2::geom_line(
    data      = fig2adj_df_row23,
    mapping   = ggplot2::aes(linetype = Approach,
                             group    = interaction(f1, Approach)),
    linewidth = LW
  ) +
  ggplot2::geom_point(
    data    = fig2adj_df_row23,
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

save_figure(Figure2, "Figure2")

# =============================================================================
# Figure 3: Sensitivity to rho (4 rows x 2 cols)
# =============================================================================

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

save_figure(Figure3, "Figure3")

# =============================================================================
# Figure 4: Null RCP under four MCP methods (rho = 0 only)
# =============================================================================

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
ref_fig4 <- rbind(
  data.frame(
    f1    = factor("f[1]==0.1", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~1)", levels = levels(fig4_long$panel)),
    yint  = fig4_ref_v$M1_01
  ),
  data.frame(
    f1    = factor("f[1]==0.1", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~2)", levels = levels(fig4_long$panel)),
    yint  = fig4_ref_v$M2_01
  ),
  data.frame(
    f1    = factor("f[1]==0.2", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~1)", levels = levels(fig4_long$panel)),
    yint  = fig4_ref_v$M1_02
  ),
  data.frame(
    f1    = factor("f[1]==0.2", levels = c("f[1]==0.1", "f[1]==0.2")),
    panel = factor("RCP[H[0]]~(Method~2)", levels = levels(fig4_long$panel)),
    yint  = fig4_ref_v$M2_02
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

save_figure(Figure4, "Figure4")

# =============================================================================
# Table 1: Application -- hypothetical MRCT (placeholder values)
# =============================================================================

APP_K     <- table1_data$APP_K
APP_ALPHA <- table1_data$APP_ALPHA

tab1_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\footnotesize",
  paste0("\\caption{Application to a hypothetical multi-regional clinical ",
         "trial with up to $K = 3$ co-primary endpoints. ",
         "True treatment effects: $\\boldsymbol{\\delta} = (0.30, 0.25, 0.20)$, ",
         "$\\sigma_{k} = 1$ ($k = 1, 2, 3$), endpoint correlation $\\rho = 0.3$, ",
         "$f_{1} = 0.15$ (region of interest), $S = 3$, ",
         "$\\alpha = 0.025$, target disjunctive power $= 0.80$. ",
         "Unadjusted thresholds: $\\gamma_{{\\rm M}1} = 0.5$, $\\gamma_{{\\rm M}2} = 0$. ",
         "Adjusted thresholds $\\gamma_{{\\rm M}1,K}$ and $\\gamma_{{\\rm M}2,K}$ are proposed herein.}"),
  "\\label{tab:application}",
  "\\begin{tabular}{r r r cc cccc cccc}",
  "\\hline",
  paste(" & & &",
        "\\multicolumn{2}{c}{Adjusted} &",
        "\\multicolumn{4}{c}{RCP, unadjusted $\\gamma$} &",
        "\\multicolumn{4}{c}{RCP, adjusted $\\gamma$} \\\\"),
  paste(" & & &",
        "\\multicolumn{2}{c}{threshold} &",
        "\\multicolumn{2}{c}{${\\rm H}_{0}$} & \\multicolumn{2}{c}{${\\rm H}_{1}$} &",
        "\\multicolumn{2}{c}{${\\rm H}_{0}$} & \\multicolumn{2}{c}{${\\rm H}_{1}$} \\\\"),
  "\\cline{4-5}\\cline{6-7}\\cline{8-9}\\cline{10-11}\\cline{12-13}",
  paste("$K$ & $N_{K}$ & $\\alpha/K$ &",
        "$\\gamma_{{\\rm M}1,K}$ & $\\gamma_{{\\rm M}2,K}$ &",
        "$\\phi^{(1)}_{K}$ & $\\phi^{(2)}_{K}$ &",
        "$\\phi^{(1)}_{K}$ & $\\phi^{(2)}_{K}$ &",
        "$\\phi^{(1)}_{K}$ & $\\phi^{(2)}_{K}$ &",
        "$\\phi^{(1)}_{K}$ & $\\phi^{(2)}_{K}$ \\\\"),
  "\\hline"
)

for (k in seq_len(APP_K)) {
  tab1_lines <- c(tab1_lines, sprintf(
    "%d & %d & %.5f & %.4f & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
    k,
    table1_data$N[k], APP_ALPHA / k,
    table1_data$gamma_M1_adj[k],
    .fmt(table1_data$gamma_M2_adj[k]),
    table1_data$null_unadj_M1[k], table1_data$null_unadj_M2[k],
    table1_data$unadj_M1[k], table1_data$unadj_M2[k],
    table1_data$null_adj_M1[k], table1_data$null_adj_M2[k],
    table1_data$adj_M1[k], table1_data$adj_M2[k]
  ))
}
tab1_lines <- c(tab1_lines, "\\hline", "\\end{tabular}", "\\end{table}")

write_tex(tab1_lines, file.path(results_dir, "table1_application.tex"))
