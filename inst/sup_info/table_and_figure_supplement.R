# =============================================================================
# table_and_figure_supplement.R
#
# Rendering stage for the Supporting Information figures of:
#   "A Cautionary Note on Regional Consistency Evaluation with Multiple
#    Primary Endpoints in Multi-Regional Clinical Trials"
#
# This script reads the intermediate data written by data_generate_supplement.R
# (data/*.rds) and produces the Supporting Information figures. It performs no
# numerical computation.
#
# Prerequisite
#   Run data_generate_supplement.R first (from this same working directory) so
#   that data/siS1_data.rds ... data/siS6_data.rds exist.
#
# Required packages: ggplot2, scales, grid
#
# Outputs (relative to the working directory)
#   results/FigureS1.eps, results/FigureS1.pdf
#   ...
#   results/FigureS6.eps, results/FigureS6.pdf
# =============================================================================

library(ggplot2)
library(scales)
library(grid)

# --------------------------------------------------------------------------- #
# Read intermediate data produced by data_generate_supplement.R
# --------------------------------------------------------------------------- #
data_dir <- "data"

siS1 <- readRDS(file.path(data_dir, "siS1_data.rds"))
siS2 <- readRDS(file.path(data_dir, "siS2_data.rds"))
siS3 <- readRDS(file.path(data_dir, "siS3_data.rds"))
siS4 <- readRDS(file.path(data_dir, "siS4_data.rds"))
siS5 <- readRDS(file.path(data_dir, "siS5_data.rds"))
siS6 <- readRDS(file.path(data_dir, "siS6_data.rds"))

s1_df <- siS1$siS1_df
s2_df <- siS2$siS2_df
s3_df <- siS3$siS3_df
s4_df <- siS4$siS4_df
s5_df <- siS5$siS5_df
s6_df <- siS6$siS6_df
K_MAX <- siS1$K_MAX

# Output location for figures (relative path, separate from data/)
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --------------------------------------------------------------------------- #
# Plotting constants (shared with the main-text figures)
# --------------------------------------------------------------------------- #
BASE_SIZE <- 24
LW        <- 1.2
PS        <- 3
LKW       <- ggplot2::unit(2.5, "cm")

OKABE <- c("#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# --------------------------------------------------------------------------- #
# Rendering helpers
# --------------------------------------------------------------------------- #
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

.theme_si <- function() {
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
}

H_LABELS <- ggplot2::as_labeller(
  c(H0 = "RCP[H[0]]", H1 = "RCP[H[1]]"),
  default = ggplot2::label_parsed
)

# =============================================================================
# Figure S1: Unknown variance vs known variance (2 rows x 2 cols)
# =============================================================================

s1_df$variance <- factor(s1_df$variance, levels = c("Known", "Unknown"))
s1_df$panel    <- factor(s1_df$panel,    levels = c("Method~1", "Method~2"))
s1_df$H        <- factor(s1_df$H,        levels = c("H0", "H1"))

FigureS1 <- ggplot2::ggplot(
  s1_df,
  ggplot2::aes(x = K, y = value, colour = variance,
               linetype = variance, group = variance)
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::facet_grid(H ~ panel, scales = "free_y",
                      labeller = ggplot2::labeller(H = H_LABELS,
                                                   panel = ggplot2::label_parsed)) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(name = "RCP") +
  ggplot2::scale_colour_manual(name = NULL, values = OKABE[c(5, 6)]) +
  ggplot2::scale_linetype_manual(name = NULL,
                                 values = c("Known" = "solid",
                                            "Unknown" = "dashed")) +
  .theme_si()

print(FigureS1)
save_figure(FigureS1, "FigureS1")

# =============================================================================
# Figure S2: Heterogeneous effect in the region of interest (1 row x 2 cols)
# =============================================================================

s2_df$panel <- factor(s2_df$panel, levels = c("Method~1", "Method~2"))
s2_df$r1f   <- factor(s2_df$r1, levels = c(0.0, 0.1, 0.2),
                      labels = c("0", "0.1", "0.2"))

FigureS2 <- ggplot2::ggplot(
  s2_df,
  ggplot2::aes(x = K, y = value, colour = r1f, group = r1f)
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::facet_wrap(~ panel, scales = "free_y",
                      labeller = ggplot2::label_parsed) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(name = expression(RCP[H[1]])) +
  ggplot2::scale_colour_manual(
    name   = expression(paste("Region-1 effect ", delta[1])),
    values = OKABE[c(6, 1, 5)]
  ) +
  .theme_si()

print(FigureS2)
save_figure(FigureS2, "FigureS2")

# =============================================================================
# Figure S3: Mixed endpoints (1 row x 2 cols)
# =============================================================================

s3_df$panel    <- factor(s3_df$panel,    levels = c("Method~1", "Method~2"))
s3_df$config   <- factor(s3_df$config,   levels = c("Homogeneous", "Mixed"))
s3_df$Approach <- factor(s3_df$Approach, levels = c("Formula", "Simulation"))

FigureS3 <- ggplot2::ggplot(
  s3_df,
  ggplot2::aes(x = K, y = value, colour = config,
               linetype = Approach, group = interaction(config, Approach))
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::facet_wrap(~ panel, scales = "free_y",
                      labeller = ggplot2::label_parsed) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(name = expression(RCP[H[1]])) +
  ggplot2::scale_colour_manual(name = NULL, values = OKABE[c(5, 6)]) +
  ggplot2::scale_linetype_manual(name = NULL,
                                 values = c("Formula" = "solid",
                                            "Simulation" = "dashed")) +
  .theme_si()

print(FigureS3)
save_figure(FigureS3, "FigureS3")

# =============================================================================
# Figure S4: Correlation combined with a null region of interest (1 row x 2 cols)
# =============================================================================

s4_df$panel <- factor(s4_df$panel, levels = c("Method~1", "Method~2"))
s4_df$rhof  <- factor(s4_df$rho, levels = c(0.0, 0.3, 0.5, 0.8),
                      labels = c("0", "0.3", "0.5", "0.8"))

FigureS4 <- ggplot2::ggplot(
  s4_df,
  ggplot2::aes(x = K, y = value, colour = rhof, group = rhof)
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::facet_wrap(~ panel, scales = "free_y",
                      labeller = ggplot2::label_parsed) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(name = expression(RCP[H[1]])) +
  ggplot2::scale_colour_manual(
    name   = expression(paste("Correlation ", rho)),
    values = OKABE[c(6, 1, 5, 3)]
  ) +
  .theme_si()

print(FigureS4)
save_figure(FigureS4, "FigureS4")

# =============================================================================
# Figure S5: Unequal endpoint variances (2 rows x 2 cols)
# =============================================================================

s5_df$panel    <- factor(s5_df$panel,    levels = c("Method~1", "Method~2"))
s5_df$H        <- factor(s5_df$H,        levels = c("H0", "H1"))
s5_df$struct   <- factor(s5_df$struct,   levels = c("Equal", "Unequal"))
s5_df$Approach <- factor(s5_df$Approach, levels = c("Formula", "Simulation"))

FigureS5 <- ggplot2::ggplot(
  s5_df,
  ggplot2::aes(x = K, y = value, colour = struct,
               linetype = Approach, group = interaction(struct, Approach))
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::facet_grid(H ~ panel, scales = "free_y",
                      labeller = ggplot2::labeller(H = H_LABELS,
                                                   panel = ggplot2::label_parsed)) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(name = "RCP") +
  ggplot2::scale_colour_manual(name = "Variances", values = OKABE[c(5, 6)]) +
  ggplot2::scale_linetype_manual(name = NULL,
                                 values = c("Formula" = "solid",
                                            "Simulation" = "dashed")) +
  .theme_si()

print(FigureS5)
save_figure(FigureS5, "FigureS5")

# =============================================================================
# Figure S6: Non-standard correlation structure (2 rows x 2 cols)
# =============================================================================

s6_df$panel    <- factor(s6_df$panel,    levels = c("Method~1", "Method~2"))
s6_df$H        <- factor(s6_df$H,        levels = c("H0", "H1"))
s6_df$struct   <- factor(s6_df$struct,   levels = c("Independence", "CS", "AR1"),
                         labels = c("Independence", "Compound symmetry", "AR(1)"))
s6_df$Approach <- factor(s6_df$Approach, levels = c("Formula", "Simulation"))

FigureS6 <- ggplot2::ggplot(
  s6_df,
  ggplot2::aes(x = K, y = value, colour = struct,
               linetype = Approach, group = interaction(struct, Approach))
) +
  ggplot2::geom_line(linewidth = LW) +
  ggplot2::geom_point(size = PS, shape = 16) +
  ggplot2::facet_grid(H ~ panel, scales = "free_y",
                      labeller = ggplot2::labeller(H = H_LABELS,
                                                   panel = ggplot2::label_parsed)) +
  ggplot2::scale_x_continuous(name = "Number of endpoints (K)",
                              breaks = seq_len(K_MAX)) +
  ggplot2::scale_y_continuous(name = "RCP") +
  ggplot2::scale_colour_manual(name = NULL, values = OKABE[c(6, 5, 3)]) +
  ggplot2::scale_linetype_manual(name = NULL,
                                 values = c("Formula" = "solid",
                                            "Simulation" = "dashed")) +
  .theme_si()

print(FigureS6)
save_figure(FigureS6, "FigureS6")
