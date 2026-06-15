#' Adjusted Consistency Thresholds under Multiple Comparison Procedures
#'
#' @description
#' Compute multiplicity-adjusted consistency thresholds \eqn{\gamma_{M1,k}}
#' and \eqn{\gamma_{M2,K}} for \eqn{K = 1, 2, \ldots, K_{\max}} endpoints,
#' such that the regional consistency probability (RCP) under the null
#' hypothesis (\eqn{\delta = 0}) remains constant across all values of
#' \eqn{K}.
#'
#' The reference RCP for each method is defined as the RCP under
#' \eqn{\delta = 0} with \eqn{K = 1} using the unadjusted thresholds
#' \code{gamma_M1} and \code{gamma_M2}. For \eqn{K \geq 2}, the adjusted
#' threshold \eqn{\gamma_{M1,K}} (resp. \eqn{\gamma_{M2,K}}) is obtained by
#' solving:
#' \deqn{
#'   \mathrm{RCP}^{(m)}_{\delta=0}(K,\, \gamma_K) =
#'   \mathrm{RCP}^{(m)}_{\delta=0}(1,\, \gamma)
#' }
#' via \code{stats::uniroot}. The adjustment for Method 1 and Method 2 is
#' performed independently.
#'
#' The adjusted thresholds can then be passed to \code{\link{rcpmcp_single}}
#' or \code{\link{rcpmcp_multiple}} to evaluate RCPs under the alternative
#' hypothesis with multiplicity-controlled consistency criteria.
#'
#' \strong{Scalar vs. vector \code{N}}:
#' \code{N} accepts either a single positive integer (the same total sample
#' size is used for all values of \eqn{K}) or an integer vector of length
#' \code{K_max} (a separate sample size \eqn{N_k} is used when computing
#' the null RCP for each \eqn{K}). The vector form is intended for use with
#' the output of \code{\link{ssmcp_multiple}}, where each \eqn{N_k} is the
#' sample size required to achieve a target disjunctive power for exactly
#' \eqn{K} endpoints. A typical workflow is:
#' \preformatted{
#' ss    <- ssmcp_multiple(delta = rep(0.2, 5), ...)
#' gamma <- rcpmcp_get_gamma(N = ss$result$N, ...)
#' res   <- rcpmcp_multiple(delta = rep(0.2, 5), N = ss$result$N,
#'                          gamma_M1 = gamma$result$gamma_M1_adj, ...)
#' }
#' Note that when \code{N} is a vector, the reference RCP at \eqn{K = 1}
#' is computed using \eqn{N_1} (the first element of \code{N}).
#'
#' @param Sigma \code{K_max}-by-\code{K_max} numeric variance-covariance
#'   matrix for the \code{K_max} endpoints. Must be symmetric and positive
#'   definite. \code{NULL} (default) sets \code{Sigma = diag(K_max)}, i.e.,
#'   unit variance and independence.
#' @param N Positive integer or integer vector of length \code{K_max}. Total
#'   sample size across all regions. When a scalar is supplied the same
#'   sample size is used for all \eqn{K}. When a vector of length
#'   \code{K_max} is supplied, \eqn{N_k} is used when solving for the
#'   adjusted threshold at \eqn{K} endpoints. All elements must be positive
#'   integers.
#' @param fs Numeric vector. Proportion of patients in each region
#'   \eqn{s = 1, \ldots, S}. Must sum to 1 and all elements must be positive.
#'   Default is \code{c(0.1, 0.45, 0.45)}.
#' @param K_max Positive integer. Maximum number of endpoints. Adjusted
#'   thresholds are computed for \eqn{K = 1, 2, \ldots, K_{\max}}.
#'   Must match the dimension of \code{Sigma}. Default is \code{5}.
#' @param gamma_M1 Numeric scalar. Unadjusted effect retention threshold for
#'   Method 1 (used as the reference at \eqn{K = 1}). Must be in
#'   \eqn{[0, 1)}. Default is \code{0.5}.
#' @param gamma_M2 Numeric scalar. Unadjusted consistency threshold for
#'   Method 2 (used as the reference at \eqn{K = 1}). Must be non-negative.
#'   Default is \code{0}.
#' @param alpha Numeric scalar. One-sided familywise significance level before
#'   multiplicity adjustment. Default is \code{0.025}.
#' @param tol Numeric scalar. Tolerance passed to \code{stats::uniroot}.
#'   Default is \code{1e-6}.
#'
#' @return An object of class \code{"rcpmcp_gamma"}, which is a list
#'   containing:
#' \describe{
#'   \item{\code{Sigma}}{Variance-covariance matrix used.}
#'   \item{\code{N}}{Total sample size(s) as supplied (scalar or vector).}
#'   \item{\code{fs}}{Regional patient proportions.}
#'   \item{\code{K_max}}{Maximum number of endpoints evaluated.}
#'   \item{\code{gamma_M1}}{Unadjusted threshold for Method 1.}
#'   \item{\code{gamma_M2}}{Unadjusted threshold for Method 2.}
#'   \item{\code{alpha}}{One-sided familywise significance level.}
#'   \item{\code{RCP0_M1}}{Reference RCP for Method 1 (\eqn{K = 1},
#'     \eqn{\delta = 0}, computed at \eqn{N_1}).}
#'   \item{\code{RCP0_M2}}{Reference RCP for Method 2 (\eqn{K = 1},
#'     \eqn{\delta = 0}, computed at \eqn{N_1}).}
#'   \item{\code{result}}{A data frame with \code{K_max} rows and columns:
#'     \code{k}, \code{N}, \code{alpha_adj}, \code{RCP0_M1} (RCP under
#'     \eqn{\delta = 0} after adjustment), \code{RCP0_M2},
#'     \code{gamma_M1_adj}, \code{gamma_M2_adj}.}
#' }
#'
#' @seealso \code{\link{rcpmcp_single}}, \code{\link{rcpmcp_multiple}},
#'   \code{\link{ssmcp_multiple}}
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats pnorm qnorm uniroot cov2cor
#'
#' @examples
#' \donttest{
#' # Example 1: Scalar N, independent endpoints
#' gamma_result <- rcpmcp_get_gamma(
#'   Sigma    = diag(5),
#'   N        = 200,
#'   fs       = c(0.1, 0.45, 0.45),
#'   K_max    = 5,
#'   gamma_M1 = 0.5,
#'   gamma_M2 = 0,
#'   alpha    = 0.025
#' )
#' print(gamma_result)
#' }
#'
#' \donttest{
#' # Example 2: Vector N from ssmcp_multiple
#' ss <- ssmcp_multiple(
#'   delta        = rep(0.2, 5),
#'   Sigma        = diag(5),
#'   fs           = c(0.1, 0.45, 0.45),
#'   K_max        = 5,
#'   alpha        = 0.025,
#'   target_power = 0.8
#' )
#' gamma_vec <- rcpmcp_get_gamma(
#'   Sigma    = diag(5),
#'   N        = ss$result$N,
#'   fs       = c(0.1, 0.45, 0.45),
#'   K_max    = 5,
#'   gamma_M1 = 0.5,
#'   gamma_M2 = 0,
#'   alpha    = 0.025
#' )
#' print(gamma_vec)
#' }
#'
#' @export
#' @param r Numeric scalar. Positive allocation ratio of the experimental
#'   group to the control group. Default is \code{1} (equal allocation).
rcpmcp_get_gamma <- function(Sigma    = NULL,
                             N,
                             fs       = c(0.1, 0.45, 0.45),
                             K_max    = 5,
                             gamma_M1 = 0.5,
                             gamma_M2 = 0,
                             alpha    = 0.025,
                             r        = 1,
                             tol      = 1e-6) {

  # ========== Input Validation ==========
  if (!is.null(Sigma)) {
    if (!is.matrix(Sigma) || nrow(Sigma) != K_max || ncol(Sigma) != K_max) {
      stop("Sigma must be a K_max-by-K_max numeric matrix")
    }
    if (!isSymmetric(Sigma)) {
      stop("Sigma must be symmetric")
    }
    if (any(eigen(Sigma, symmetric = TRUE)$values <= 0)) {
      stop("Sigma must be positive definite")
    }
    if (any(diag(Sigma) <= 0)) {
      stop("diagonal elements of Sigma must be positive")
    }
  }
  # N: scalar or vector of length K_max, all positive integers
  if (!is.numeric(N) || any(N <= 0) || any(N != as.integer(N))) {
    stop("N must be a positive integer or a vector of positive integers")
  }
  if (!length(N) %in% c(1L, K_max)) {
    stop("N must be a scalar or a vector of length K_max")
  }
  if (!is.numeric(fs) || any(fs <= 0) || abs(sum(fs) - 1) > 1e-8) {
    stop("fs must be a numeric vector of positive values summing to 1")
  }
  if (!is.numeric(K_max) || length(K_max) != 1 ||
      K_max < 1 || K_max != as.integer(K_max)) {
    stop("K_max must be a single positive integer")
  }
  if (!is.numeric(gamma_M1) || length(gamma_M1) != 1 ||
      gamma_M1 < 0 || gamma_M1 >= 1) {
    stop("gamma_M1 must be a single numeric value in [0, 1)")
  }
  if (!is.numeric(gamma_M2) || length(gamma_M2) != 1 || gamma_M2 < 0) {
    stop("gamma_M2 must be a single non-negative numeric value")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(r) || length(r) != 1 || r <= 0) {
    stop("r must be a single positive numeric value")
  }
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("tol must be a single positive numeric value")
  }

  # ========== Setup ==========
  Sigma_use <- if (is.null(Sigma)) diag(K_max) else Sigma
  c_alloc   <- (r + 1)^2 / r
  Sigma_var <- c_alloc * Sigma_use

  # Expand scalar N to a vector of length K_max for uniform indexing
  N_vec <- if (length(N) == 1L) rep(N, K_max) else N

  # ========== Reference RCP at K = 1 (uses N_vec[1]) ==========
  # Pre-compute fixed quantities for K = 1, delta = 0
  pre1 <- rcpmcp_single_precompute(
    delta     = rep(0, 1),
    Sigma_use = Sigma_var[1L, 1L, drop = FALSE],
    N         = N_vec[1L],
    fs        = fs,
    K         = 1L,
    alpha     = alpha
  )
  # Power at K = 1, delta = 0 (= alpha, but computed consistently)
  power1 <- .rcpmcp_formula_power(pre1)

  # RCP at K = 1 with unadjusted thresholds
  ref     <- .rcpmcp_formula_rcp_only(
    pre       = pre1,
    delta     = rep(0, 1),
    Sigma_use = Sigma_var[1L, 1L, drop = FALSE],
    gamma_M1  = gamma_M1,
    gamma_M2  = gamma_M2,
    power_out = power1
  )
  rcp0_m1   <- ref$RCP_M1
  rcp0_m2   <- ref$RCP_M2
  alpha_adj1 <- pre1$alpha_adj

  # ========== Solve for adjusted thresholds at K = 1, ..., K_max ==========
  result <- data.frame(
    K            = integer(K_max),
    N            = integer(K_max),
    alpha_adj    = numeric(K_max),
    RCP0_M1      = numeric(K_max),
    RCP0_M2      = numeric(K_max),
    gamma_M1_adj = numeric(K_max),
    gamma_M2_adj = numeric(K_max)
  )

  # K = 1: no adjustment needed
  result[1L, ] <- list(
    K            = 1L,
    N            = N_vec[1L],
    alpha_adj    = alpha_adj1,
    RCP0_M1      = rcp0_m1,
    RCP0_M2      = rcp0_m2,
    gamma_M1_adj = gamma_M1,
    gamma_M2_adj = gamma_M2
  )

  for (k in seq(2L, K_max)) {

    idx_k     <- seq_len(k)
    Sigma_k   <- Sigma_var[idx_k, idx_k, drop = FALSE]
    delta_k   <- rep(0, k)

    # ------------------------------------------------------------------
    # Pre-compute all gamma-independent quantities for this k once.
    # Power (= denominator of RCP) depends only on k, N, fs, alpha, Sigma,
    # and delta -- not on gamma_M1 / gamma_M2 -- so it is computed here
    # once and cached, eliminating redundant pmvnorm calls inside uniroot.
    # ------------------------------------------------------------------
    pre_k <- rcpmcp_single_precompute(
      delta     = delta_k,
      Sigma_use = Sigma_k,
      N         = N_vec[k],
      fs        = fs,
      K         = k,
      alpha     = alpha
    )
    power_k <- .rcpmcp_formula_power(pre_k)

    # ---- Method 1: solve for gamma_M1_adj in [gamma_M1, 1 - tol] ----
    # f(g) = RCP_M1_H0(k, g) - rcp0_m1 = 0
    # RCP_M1 decreases monotonically as g increases.
    # f(gamma_M1) > 0  (RCP inflated at k > 1 with unadjusted threshold)
    # f(1 - tol)  < 0  (RCP near 0 as g approaches 1)
    f_m1 <- function(g) {
      .rcpmcp_formula_rcp_only(
        pre       = pre_k,
        delta     = delta_k,
        Sigma_use = Sigma_k,
        gamma_M1  = g,
        gamma_M2  = gamma_M2,
        power_out = power_k,
        which     = "m1"
      )$RCP_M1 - rcp0_m1
    }

    gm1_adj <- tryCatch(
      stats::uniroot(f_m1, lower = gamma_M1, upper = 1 - tol,
                     tol = tol)$root,
      error = function(e) {
        warning(sprintf(
          "uniroot failed for Method 1 at K = %d: %s",
          k, conditionMessage(e)
        ))
        NA_real_
      }
    )

    # ---- Method 2: solve for gamma_M2_adj in [gamma_M2, upper_M2] ----
    # Determine a finite upper bound where RCP_M2 drops below rcp0_m2.
    # Start at max(gamma_M2 + 1, 1) and double until f < 0.
    # A generous initial value avoids excessive bracket-expansion iterations.
    upper_m2 <- max(gamma_M2 + 1, 1)
    for (i in seq_len(30L)) {
      rcp_upper <- tryCatch(
        .rcpmcp_formula_rcp_only(
          pre       = pre_k,
          delta     = delta_k,
          Sigma_use = Sigma_k,
          gamma_M1  = gamma_M1,
          gamma_M2  = upper_m2,
          power_out = power_k,
          which     = "m2"
        )$RCP_M2,
        error = function(e) NA_real_
      )
      if (!is.na(rcp_upper) && !is.nan(rcp_upper) && rcp_upper < rcp0_m2) break
      upper_m2 <- upper_m2 * 2
    }

    f_m2 <- function(g) {
      .rcpmcp_formula_rcp_only(
        pre       = pre_k,
        delta     = delta_k,
        Sigma_use = Sigma_k,
        gamma_M1  = gamma_M1,
        gamma_M2  = g,
        power_out = power_k,
        which     = "m2"
      )$RCP_M2 - rcp0_m2
    }

    gm2_adj <- tryCatch(
      stats::uniroot(f_m2, lower = gamma_M2, upper = upper_m2,
                     tol = tol)$root,
      error = function(e) {
        warning(sprintf(
          "uniroot failed for Method 2 at K = %d: %s",
          k, conditionMessage(e)
        ))
        NA_real_
      }
    )

    # Record adjusted RCPs under H0 for verification (one call each)
    rcp_adj_m1 <- if (!is.na(gm1_adj)) {
      .rcpmcp_formula_rcp_only(
        pre       = pre_k,
        delta     = delta_k,
        Sigma_use = Sigma_k,
        gamma_M1  = gm1_adj,
        gamma_M2  = gamma_M2,
        power_out = power_k,
        which     = "m1"
      )$RCP_M1
    } else NA_real_

    rcp_adj_m2 <- if (!is.na(gm2_adj)) {
      .rcpmcp_formula_rcp_only(
        pre       = pre_k,
        delta     = delta_k,
        Sigma_use = Sigma_k,
        gamma_M1  = gamma_M1,
        gamma_M2  = gm2_adj,
        power_out = power_k,
        which     = "m2"
      )$RCP_M2
    } else NA_real_

    result[k, ] <- list(
      K            = k,
      N            = N_vec[k],
      alpha_adj    = alpha / k,
      RCP0_M1      = rcp_adj_m1,
      RCP0_M2      = rcp_adj_m2,
      gamma_M1_adj = gm1_adj,
      gamma_M2_adj = gm2_adj
    )
  }

  # ========== Output ==========
  out <- list(
    Sigma    = Sigma_use,
    N        = N,
    fs       = fs,
    K_max    = K_max,
    gamma_M1 = gamma_M1,
    gamma_M2 = gamma_M2,
    alpha    = alpha,
    r        = r,
    RCP0_M1  = rcp0_m1,
    RCP0_M2  = rcp0_m2,
    result   = result
  )
  class(out) <- "rcpmcp_gamma"
  return(out)
}


#' Print Method for rcpmcp_gamma Objects
#'
#' @param x An object of class \code{"rcpmcp_gamma"}.
#' @param digits A non-negative integer specifying the number of decimal
#'   places. Default is \code{4}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @rdname rcpmcp_get_gamma
#' @export
print.rcpmcp_gamma <- function(x, digits = 4, ...) {
  cat("\nAdjusted Consistency Thresholds under Multiple Comparison Procedures\n")
  cat(strrep("-", 70), "\n")

  # Display N: scalar or vector
  if (length(x$N) == 1L) {
    cat(sprintf("   Total sample size : N        = %d\n", x$N))
  } else {
    cat(sprintf("   Total sample size : N        = (%s)\n",
                paste(x$N, collapse = ", ")))
  }

  cat(sprintf("   Regional props.   : fs       = (%s)\n",
              paste(x$fs, collapse = ", ")))
  cat(sprintf("   Max. endpoints    : K_max    = %d\n",   x$K_max))
  cat(sprintf("   Base threshold    : gamma_M1 = %.4f  (Method 1)\n",
              x$gamma_M1))
  cat(sprintf("   Base threshold    : gamma_M2 = %.4f  (Method 2)\n",
              x$gamma_M2))
  cat(sprintf("   Significance lvl  : alpha    = %.4f\n", x$alpha))
  cat(sprintf("   Allocation ratio  : r        = %s\n", if (is.null(x$r)) 1 else x$r))

  corr_use <- cov2cor(x$Sigma)
  is_indep <- all(abs(corr_use - diag(x$K_max)) < 1e-8)
  if (is_indep) {
    cat("   Endpoint corr.    : Independent (identity matrix)\n")
  } else {
    cat("   Endpoint corr.    : User-specified (see $Sigma)\n")
  }

  cat(sprintf(
    "\n   Reference RCP (K=1, delta=0) : RCP_M1 = %.4f,  RCP_M2 = %.4f\n",
    x$RCP0_M1, x$RCP0_M2
  ))

  cat(strrep("-", 70), "\n")
  cat("Adjusted thresholds (RCP0 = RCP under delta=0):\n\n")

  df <- x$result
  df$alpha_adj    <- formatC(df$alpha_adj,    format = "f", digits = digits + 1)
  df$RCP0_M1      <- formatC(df$RCP0_M1,      format = "f", digits = digits)
  df$RCP0_M2      <- formatC(df$RCP0_M2,      format = "f", digits = digits)
  df$gamma_M1_adj <- formatC(df$gamma_M1_adj, format = "f", digits = digits)
  df$gamma_M2_adj <- formatC(df$gamma_M2_adj, format = "f", digits = digits)

  colnames(df) <- c("K", "N", "alpha/K",
                    "RCP0 (M1)", "RCP0 (M2)",
                    "gamma_M1_adj", "gamma_M2_adj")
  print(df, row.names = FALSE, right = TRUE)

  cat("\n")
  invisible(x)
}


#' @rdname rcpmcp_get_gamma
#'
#' @description
#' \code{plot.rcpmcp_gamma}: Visualise multiplicity-adjusted consistency
#' thresholds (\eqn{\gamma_{M1,K}}, \eqn{\gamma_{M2,K}}) and the
#' corresponding null RCPs (\eqn{\mathrm{RCP0}_{M1}},
#' \eqn{\mathrm{RCP0}_{M2}}) as a function of the number of endpoints
#' \eqn{K}. Strip labels use \code{ggplot2::label_parsed} so that
#' mathematical expressions such as \eqn{\gamma_{M1}} render correctly.
#' When additional \code{rcpmcp_gamma} objects are supplied via
#' \code{overlay}, all curves are overlaid on a single figure with distinct
#' colours.
#'
#' Available panels (controlled by the \code{panels} argument):
#' \itemize{
#'   \item \code{"gamma_M1"}: Adjusted effect-retention threshold (Method 1).
#'   \item \code{"gamma_M2"}: Adjusted consistency threshold (Method 2).
#'   \item \code{"RCP0_M1"}:  Null RCP after adjustment (Method 1); should
#'     be flat at the reference value.
#'   \item \code{"RCP0_M2"}:  Null RCP after adjustment (Method 2); should
#'     be flat at the reference value.
#' }
#'
#' @param overlay An optional \strong{named list} of additional
#'   \code{"rcpmcp_gamma"} objects to overlay. All objects must share the
#'   same \code{K_max} as \code{x}. Default is \code{NULL}.
#' @param group_labels Character vector or \code{expression}. Legend labels
#'   for all conditions in the order \code{c(x, overlay[[1]], \ldots)}.
#'   Supports \code{expression()} for mathematical annotation, e.g.
#'   \code{expression(rho == 0, rho == 0.5)}.
#'   Default is \code{NULL} (auto-generated labels).
#' @param panels Character vector selecting which panels to display. Any
#'   non-empty subset of \code{c("gamma_M1", "gamma_M2", "RCP0_M1",
#'   "RCP0_M2")}. The display order is fixed as listed. Default shows all
#'   four panels.
#' @param show_reference Logical. If \code{TRUE} (default), horizontal
#'   reference lines are drawn at the unadjusted base thresholds (threshold
#'   panels) and the reference null RCPs (null-RCP panels).
#' @param layout Integer vector of length 2 specifying panel arrangement as
#'   \code{c(nrow, ncol)}. For example, \code{c(1, 4)} produces 1 row and
#'   4 columns, \code{c(2, 2)} produces a 2-by-2 grid, and \code{c(4, 1)}
#'   (the default) produces 4 rows and 1 column. When \code{ncol > 1}, strip
#'   labels are automatically moved to the top of each panel.
#' @param fixed_rcp_ylim Logical. If \code{TRUE}, the y-axis of
#'   \code{"RCP0_M1"} and \code{"RCP0_M2"} panels is fixed to \eqn{[0, 1]},
#'   while threshold panels retain free scaling. Requires the
#'   \pkg{ggh4x} package. Default is \code{FALSE}.
#' @param colours Character vector of colour codes. Recycled if shorter than
#'   the number of conditions. Default uses the Okabe-Ito palette
#'   (colour-blind friendly).
#' @param point_size Numeric scalar. Size of the points. Default \code{2.5}.
#' @param line_width Numeric scalar. Width of the lines. Default \code{0.8}.
#' @param base_size Numeric scalar. Base font size for
#'   \code{ggplot2::theme_bw}. Default \code{11}.
#'
#' @return \code{plot.rcpmcp_gamma} returns a \code{ggplot} object.
#'   Save with \code{ggplot2::ggsave()}.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline facet_wrap label_parsed scale_x_continuous scale_y_continuous scale_colour_manual theme_bw theme element_text element_rect element_blank margin unit
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#' # plot(): Single object, all four panels (default)
#' gamma_res <- rcpmcp_get_gamma(
#'   Sigma = diag(5), N = 200,
#'   fs = c(0.1, 0.45, 0.45), K_max = 5,
#'   gamma_M1 = 0.5, gamma_M2 = 0, alpha = 0.025
#' )
#' p <- plot(gamma_res)
#' print(p)
#'
#' # plot(): Threshold panels only
#' p2 <- plot(gamma_res, panels = c("gamma_M1", "gamma_M2"))
#' print(p2)
#' }
#'
#' \donttest{
#' # plot(): Overlay multiple conditions with Greek rho in legend
#' make_sigma <- function(rho, k) {
#'   S <- matrix(rho, k, k); diag(S) <- 1; S
#' }
#' gamma_res_overlay <- rcpmcp_get_gamma(
#'   Sigma = diag(5), N = 200,
#'   fs = c(0.1, 0.45, 0.45), K_max = 5,
#'   gamma_M1 = 0.5, gamma_M2 = 0, alpha = 0.025
#' )
#' gamma05 <- rcpmcp_get_gamma(
#'   Sigma = make_sigma(0.5, 5), N = 200,
#'   fs = c(0.1, 0.45, 0.45), K_max = 5
#' )
#' p3 <- plot(gamma_res_overlay,
#'            overlay      = list(gamma05),
#'            group_labels = expression(rho == 0, rho == 0.5))
#' print(p3)
#' }
#'
#' @export
plot.rcpmcp_gamma <- function(x,
                              ...,
                              overlay        = NULL,
                              group_labels   = NULL,
                              panels         = c("gamma_M1", "gamma_M2",
                                                 "RCP0_M1",  "RCP0_M2"),
                              show_reference = TRUE,
                              layout         = c(length(panels), 1L),
                              fixed_rcp_ylim = FALSE,
                              colours        = NULL,
                              point_size     = 2.5,
                              line_width     = 0.8,
                              base_size      = 11) {

  # ========== Build ordered object list ==========
  if (!is.null(overlay)) {
    if (!is.list(overlay) ||
        !all(vapply(overlay, inherits, logical(1L), "rcpmcp_gamma"))) {
      stop("overlay must be a list of 'rcpmcp_gamma' objects.")
    }
    obj_list <- c(list(x), overlay)
  } else {
    obj_list <- list(x)
  }
  n_cond       <- length(obj_list)
  single_input <- (n_cond == 1L)

  # ========== Internal keys and legend labels ==========
  internal_keys <- paste0(".cond", seq_len(n_cond))

  if (!is.null(group_labels)) {
    n_lab <- if (is.expression(group_labels)) length(group_labels) else
      length(group_labels)
    if (n_lab != n_cond) {
      stop(paste0(
        "length(group_labels) must equal 1 + length(overlay) = ", n_cond, "."
      ))
    }
    legend_labels <- group_labels
  } else {
    legend_labels <- paste0("Condition ", seq_len(n_cond))
  }
  names(obj_list) <- internal_keys

  # ========== Validate panels argument ==========
  valid_panels <- c("gamma_M1", "gamma_M2", "RCP0_M1", "RCP0_M2")
  panels       <- match.arg(panels, valid_panels, several.ok = TRUE)
  if (length(panels) == 0L) {
    stop("panels must contain at least one valid panel name.")
  }

  # ========== Validate layout argument ==========
  if (!is.numeric(layout) || length(layout) != 2L ||
      any(layout < 1L) || any(layout != as.integer(layout))) {
    stop("layout must be an integer vector of length 2 with positive values.")
  }
  nrow_layout <- as.integer(layout[1L])
  ncol_layout <- as.integer(layout[2L])

  # ========== Validate fixed_rcp_ylim argument ==========
  if (!is.logical(fixed_rcp_ylim) || length(fixed_rcp_ylim) != 1L) {
    stop("fixed_rcp_ylim must be a single logical value.")
  }
  if (fixed_rcp_ylim && !requireNamespace("ggh4x", quietly = TRUE)) {
    stop("fixed_rcp_ylim = TRUE requires the 'ggh4x' package. ",
         "Install it with: install.packages(\"ggh4x\")")
  }

  # ========== Validate K_max consistency ==========
  K_vals <- vapply(obj_list, function(o) as.integer(o$K_max), integer(1L))
  if (length(unique(K_vals)) > 1L) {
    stop("All rcpmcp_gamma objects must have the same K_max.")
  }
  K_max <- unique(K_vals)

  # ========== Panel strip labels (parsed for math rendering) ==========
  # label_parsed() renders strip text as R plotmath expressions.
  # The factor levels must therefore be valid plotmath strings.
  panel_plotmath <- c(
    gamma_M1 = "gamma[M1]~(adjusted)",
    gamma_M2 = "gamma[M2]~(adjusted)",
    RCP0_M1  = "RCP[0]~(Method~1)",
    RCP0_M2  = "RCP[0]~(Method~2)"
  )
  panel_col_map <- c(
    gamma_M1 = "gamma_M1_adj",
    gamma_M2 = "gamma_M2_adj",
    RCP0_M1  = "RCP0_M1",
    RCP0_M2  = "RCP0_M2"
  )
  active_panels <- valid_panels[valid_panels %in% panels]
  # Use plotmath strings as factor levels so label_parsed can parse them
  panel_levels  <- unname(panel_plotmath[active_panels])

  # ========== Build tidy long data frame ==========
  df_list <- lapply(internal_keys, function(key) {
    res <- obj_list[[key]]$result
    do.call(rbind, lapply(active_panels, function(pnl) {
      data.frame(
        condition = key,
        K         = res$K,
        panel     = unname(panel_plotmath[pnl]),
        value     = res[[panel_col_map[pnl]]],
        stringsAsFactors = FALSE
      )
    }))
  })
  df_long           <- do.call(rbind, df_list)
  df_long$panel     <- factor(df_long$panel,     levels = panel_levels)
  df_long$condition <- factor(df_long$condition, levels = internal_keys)

  # ========== Colours ==========
  # Okabe-Ito palette (colour-blind friendly, first 8 colours)
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#000000"
  )
  col_use <- if (!is.null(colours)) {
    rep_len(colours, n_cond)
  } else {
    okabe_ito[seq_len(n_cond)]
  }
  names(col_use) <- internal_keys

  # ========== Reference lines ==========
  # threshold panels -> unadjusted base threshold from the first object
  # RCP0 panels      -> stored reference RCP scalar from the first object
  ref_df <- NULL
  if (show_reference) {
    o_ref <- obj_list[[internal_keys[1L]]]
    ref_map <- c(
      gamma_M1 = o_ref$gamma_M1,
      gamma_M2 = o_ref$gamma_M2,
      RCP0_M1  = o_ref$RCP0_M1,
      RCP0_M2  = o_ref$RCP0_M2
    )
    ref_rows <- lapply(active_panels, function(pnl) {
      data.frame(
        panel = unname(panel_plotmath[pnl]),
        yint  = unname(ref_map[pnl]),
        stringsAsFactors = FALSE
      )
    })
    ref_df       <- do.call(rbind, ref_rows)
    ref_df$panel <- factor(ref_df$panel, levels = panel_levels)
  }

  # ========== Strip position: left when single column, top otherwise ==========
  strip_pos <- if (ncol_layout == 1L) "left" else "top"

  # ========== Build base plot ==========
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x      = .data$K,
      y      = .data$value,
      colour = .data$condition,
      group  = .data$condition
    )
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(
      ~ panel,
      nrow           = nrow_layout,
      ncol           = ncol_layout,
      scales         = "free_y",
      strip.position = strip_pos,
      labeller       = ggplot2::label_parsed
    ) +
    ggplot2::scale_x_continuous(
      name   = "Number of endpoints (K)",
      breaks = seq_len(K_max)
    ) +
    ggplot2::scale_y_continuous(name = NULL) +
    ggplot2::scale_colour_manual(
      name   = if (single_input) NULL else "Condition",
      values = col_use,
      labels = if (single_input) ggplot2::waiver() else legend_labels
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text         = ggplot2::element_text(face = "bold",
                                                 size = base_size),
      strip.placement    = "outside",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.spacing      = ggplot2::unit(0.8, "lines"),
      legend.position    = if (single_input) "none" else "bottom",
      legend.title       = ggplot2::element_text(face = "bold"),
      axis.title.x       = ggplot2::element_text(
        margin = ggplot2::margin(t = 6)
      ),
      plot.margin        = ggplot2::margin(8, 8, 8, 8)
    )

  # ========== Reference lines layer ==========
  if (!is.null(ref_df)) {
    p <- p + ggplot2::geom_hline(
      data        = ref_df,
      ggplot2::aes(yintercept = .data$yint),
      linetype    = "dashed",
      colour      = "grey40",
      linewidth   = 0.6,
      inherit.aes = FALSE
    )
  }

  # ========== Per-panel y-axis limits via ggh4x ==========
  # When fixed_rcp_ylim = TRUE, RCP0 panels are fixed to [0, 1] while
  # threshold panels retain free_y scaling.
  if (fixed_rcp_ylim) {
    rcp_plotmath <- unname(panel_plotmath[c("RCP0_M1", "RCP0_M2")])
    scales_list  <- lapply(panel_levels, function(lv) {
      if (lv %in% rcp_plotmath) {
        ggplot2::scale_y_continuous(limits = c(0, 1))
      } else {
        ggplot2::scale_y_continuous(limits = NULL)
      }
    })
    p <- p + ggh4x::facetted_pos_scales(y = scales_list)
  }

  p
}
