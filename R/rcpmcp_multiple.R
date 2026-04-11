#' Regional Consistency Probability under Multiple Comparison Procedures
#' (Multiple Evaluation)
#'
#' @description
#' Compute power and regional consistency probabilities (RCPs) for
#' \eqn{K = 1, 2, \ldots, K_{\max}} endpoints by repeatedly calling
#' \code{\link{rcpmcp_single}}. This function is intended to investigate
#' how power and RCPs change as the number of endpoints increases under a
#' specified multiplicity adjustment procedure.
#'
#' For each value of \eqn{K}, \code{\link{rcpmcp_single}} is called with
#' the first \eqn{K} elements of \code{delta} and the upper-left
#' \eqn{K}-by-\eqn{K} submatrix of \code{Sigma}. Results are collected
#' into a single data frame (formula) or a named list of data frames
#' (simulation) for easy inspection and plotting.
#'
#' When \code{approach = "formula"}, the Bonferroni closed-form solution is
#' always used. When \code{approach = "simulation"}, all four multiplicity
#' adjustment procedures (\code{"bonferroni"}, \code{"holm"},
#' \code{"hochberg"}, \code{"hommel"}) are computed simultaneously from the
#' same random draws at each \eqn{K}, and results are stored in a named list
#' \code{$result} with one data frame per procedure.
#'
#' \strong{Scalar vs. vector \code{N}}:
#' \code{N} accepts either a single positive integer (the same total sample
#' size is used for all values of \eqn{K}) or an integer vector of length
#' \code{K_max} (a separate sample size \eqn{N_k} is used for each
#' \eqn{K = 1, \ldots, K_{\max}}). The vector form is intended for use with
#' the output of \code{\link{ssmcp_multiple}}, where each \eqn{N_k} is the
#' sample size required to achieve a target disjunctive power for exactly
#' \eqn{K} endpoints. A typical workflow is:
#' \preformatted{
#' ss  <- ssmcp_multiple(delta = rep(0.2, 5), ...)
#' res <- rcpmcp_multiple(delta = rep(0.2, 5), N = ss$result$N, ...)
#' }
#'
#' \strong{Scalar vs. vector \code{gamma_M1} / \code{gamma_M2}}:
#' Both arguments accept either a single scalar (the same threshold is used
#' for all values of \eqn{K}) or a numeric vector of length \code{K_max}
#' (a separate threshold \eqn{\gamma_{M1,K}} or \eqn{\gamma_{M2,K}} is used
#' for each \eqn{K}). The vector form is intended for use with the output of
#' \code{\link{rcpmcp_get_gamma}}, enabling direct comparison of
#' multiplicity-adjusted versus unadjusted thresholds:
#' \preformatted{
#' gamma_res <- rcpmcp_get_gamma(...)
#' res_adj   <- rcpmcp_multiple(...,
#'                gamma_M1 = gamma_res$result$gamma_M1_adj,
#'                gamma_M2 = gamma_res$result$gamma_M2_adj)
#' res_unadj <- rcpmcp_multiple(..., gamma_M1 = 0.5, gamma_M2 = 0)
#' plot(res_adj, overlay      = list(res_unadj),
#'               group_labels = c("Adjusted", "Unadjusted"))
#' }
#'
#' @param delta Numeric vector of length \code{K_max}. True treatment effects
#'   (mean differences) for each endpoint. All elements must be non-negative.
#'   Zero values correspond to the null hypothesis and can be used to evaluate
#'   type I error.
#' @param Sigma \code{K_max}-by-\code{K_max} numeric variance-covariance
#'   matrix for the \code{K_max} endpoints. Must be symmetric and positive
#'   definite. \code{NULL} (default) sets \code{Sigma = diag(K_max)}, i.e.,
#'   unit variance and independence.
#' @param N Positive integer or integer vector of length \code{K_max}. Total
#'   sample size across all regions. When a scalar is supplied the same
#'   sample size is used for all \eqn{K}. When a vector of length
#'   \code{K_max} is supplied, \eqn{N_k} is used for the evaluation with
#'   \eqn{K} endpoints. All elements must be positive integers.
#' @param fs Numeric vector. Proportion of patients in each region
#'   \eqn{s = 1, \ldots, S}. Must sum to 1 and all elements must be positive.
#'   Default is \code{c(0.1, 0.45, 0.45)}.
#' @param K_max Positive integer. Maximum number of endpoints. Results are
#'   computed for \eqn{K = 1, 2, \ldots, K_{\max}}. Must match the length of
#'   \code{delta} and the dimension of \code{Sigma}. Default is \code{5}.
#' @param gamma_M1 Numeric scalar or vector of length \code{K_max}. Effect
#'   retention threshold for Method 1. Region 1 must satisfy
#'   \eqn{\hat{\delta}_{1j} > \gamma_{M1} \hat{\delta}_j}.
#'   When a scalar is supplied the same threshold is used for all \eqn{K}.
#'   When a vector of length \code{K_max} is supplied, \eqn{\gamma_{M1,K}}
#'   is used for the evaluation with \eqn{K} endpoints, enabling
#'   multiplicity-adjusted thresholds from \code{\link{rcpmcp_get_gamma}}.
#'   All elements must be in \eqn{[0, 1]}. Default is \code{0.5}.
#' @param gamma_M2 Numeric scalar or vector of length \code{K_max}.
#'   Consistency threshold for Method 2. All regional estimates must satisfy
#'   \eqn{\hat{\delta}_{sj} > \gamma_{M2}}.
#'   When a scalar is supplied the same threshold is used for all \eqn{K}.
#'   When a vector of length \code{K_max} is supplied, \eqn{\gamma_{M2,K}}
#'   is used for the evaluation with \eqn{K} endpoints.
#'   All elements must be non-negative. Default is \code{0}.
#' @param alpha Numeric scalar. One-sided familywise significance level before
#'   multiplicity adjustment. Default is \code{0.025}.
#' @param approach Character scalar. Calculation approach: \code{"formula"}
#'   for the closed-form Bonferroni solution, or \code{"simulation"} for
#'   Monte Carlo simulation. Default is \code{"formula"}.
#' @param nsim Positive integer. Number of Monte Carlo iterations. Used only
#'   when \code{approach = "simulation"}. Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for reproducibility. Used
#'   only when \code{approach = "simulation"}. Default is \code{1}.
#'
#' @return An object of class \code{"rcpmcp_multiple"}, which is a list
#'   containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used.}
#'   \item{\code{nsim}}{Number of Monte Carlo iterations (\code{NULL} for
#'     the \code{"formula"} approach).}
#'   \item{\code{delta}}{True treatment effects.}
#'   \item{\code{Sigma}}{Variance-covariance matrix used.}
#'   \item{\code{N}}{Total sample size(s) as supplied (scalar or vector).}
#'   \item{\code{fs}}{Regional patient proportions.}
#'   \item{\code{K_max}}{Maximum number of endpoints evaluated.}
#'   \item{\code{gamma_M1}}{Effect retention threshold(s) for Method 1
#'     as supplied (scalar or vector of length \code{K_max}).}
#'   \item{\code{gamma_M2}}{Consistency threshold(s) for Method 2
#'     as supplied (scalar or vector of length \code{K_max}).}
#'   \item{\code{alpha}}{One-sided familywise significance level.}
#'   \item{\code{result}}{
#'     For \code{approach = "formula"}: a data frame with \code{K_max} rows
#'     and columns \code{k}, \code{N}, \code{alpha_adj}, \code{gamma_M1},
#'     \code{gamma_M2}, \code{Power}, \code{RCP_M1}, \code{RCP_M2}.
#'     For \code{approach = "simulation"}: a named list with four elements
#'     (\code{"bonferroni"}, \code{"holm"}, \code{"hochberg"},
#'     \code{"hommel"}), each a data frame of the same structure computed
#'     from the same random draws.}
#' }
#'
#' @seealso \code{\link{rcpmcp_single}}, \code{\link{ssmcp_multiple}},
#'   \code{\link{rcpmcp_get_gamma}}
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats p.adjust pnorm qnorm rnorm
#'
#' @examples
#' # Example 1: Scalar N, closed-form solution
#' result_formula <- rcpmcp_multiple(
#'   delta    = rep(0.2, 5),
#'   Sigma    = diag(5),
#'   N        = 200,
#'   fs       = c(0.1, 0.45, 0.45),
#'   K_max    = 5,
#'   gamma_M1 = 0.5,
#'   gamma_M2 = 0,
#'   alpha    = 0.025,
#'   approach = "formula"
#' )
#' print(result_formula)
#'
#' # Example 2: Simulation (all four MCP methods computed simultaneously)
#' result_sim <- rcpmcp_multiple(
#'   delta    = rep(0.2, 5),
#'   Sigma    = diag(5),
#'   N        = 200,
#'   fs       = c(0.1, 0.45, 0.45),
#'   K_max    = 5,
#'   gamma_M1 = 0.5,
#'   gamma_M2 = 0,
#'   alpha    = 0.025,
#'   approach = "simulation",
#'   nsim     = 10000,
#'   seed     = 1
#' )
#' print(result_sim)
#'
#' # Example 3: Print selected MCP methods only
#' print(result_sim, mcp_method = c("bonferroni", "holm"))
#'
#' @export
rcpmcp_multiple <- function(delta,
                            Sigma    = NULL,
                            N,
                            fs       = c(0.1, 0.45, 0.45),
                            K_max    = 5,
                            gamma_M1 = 0.5,
                            gamma_M2 = 0,
                            alpha    = 0.025,
                            approach = "formula",
                            nsim     = 1e4,
                            seed     = 1) {

  # ========== Input Validation ==========
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("delta must be a numeric vector of non-negative values")
  }
  if (length(delta) != K_max) {
    stop("length of delta must equal K_max")
  }
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
  # gamma_M1: scalar or vector of length K_max, values in [0, 1]
  if (!is.numeric(gamma_M1) || any(gamma_M1 < 0) || any(gamma_M1 > 1)) {
    stop("gamma_M1 must be numeric with all values in [0, 1]")
  }
  if (!length(gamma_M1) %in% c(1L, K_max)) {
    stop("gamma_M1 must be a scalar or a vector of length K_max")
  }
  # gamma_M2: scalar or vector of length K_max, non-negative
  if (!is.numeric(gamma_M2) || any(gamma_M2 < 0)) {
    stop("gamma_M2 must be numeric with all non-negative values")
  }
  if (!length(gamma_M2) %in% c(1L, K_max)) {
    stop("gamma_M2 must be a scalar or a vector of length K_max")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1)")
  }
  if (!approach %in% c("formula", "simulation")) {
    stop('approach must be either "formula" or "simulation"')
  }
  if (approach == "simulation") {
    if (!is.numeric(nsim) || length(nsim) != 1 ||
        nsim <= 0 || nsim != as.integer(nsim)) {
      stop("nsim must be a single positive integer")
    }
    if (!is.numeric(seed) || length(seed) != 1 ||
        seed < 0 || seed != as.integer(seed)) {
      stop("seed must be a single non-negative integer")
    }
  }

  # ========== Setup ==========
  Sigma_use <- if (is.null(Sigma)) diag(K_max) else Sigma

  # Expand scalar inputs to vectors of length K_max for uniform indexing
  N_vec        <- if (length(N)        == 1L) rep(N,        K_max) else N
  gamma_M1_vec <- if (length(gamma_M1) == 1L) rep(gamma_M1, K_max) else gamma_M1
  gamma_M2_vec <- if (length(gamma_M2) == 1L) rep(gamma_M2, K_max) else gamma_M2

  # ========== Helper: build an empty result data frame ==========
  .empty_result_df <- function() {
    data.frame(
      K         = integer(K_max),
      N         = integer(K_max),
      alpha_adj = numeric(K_max),
      gamma_M1  = numeric(K_max),
      gamma_M2  = numeric(K_max),
      Power     = numeric(K_max),
      RCP_M1    = numeric(K_max),
      RCP_M2    = numeric(K_max)
    )
  }

  # ========== Formula Approach ==========
  if (approach == "formula") {

    result <- .empty_result_df()

    for (k in seq_len(K_max)) {
      idx_k     <- seq_len(k)
      delta_k   <- delta[idx_k]
      Sigma_k   <- Sigma_use[idx_k, idx_k, drop = FALSE]
      Sigma_arg <- if (all(Sigma_k == diag(k))) NULL else Sigma_k

      res_k <- rcpmcp_single(
        delta    = delta_k,
        Sigma    = Sigma_arg,
        N        = N_vec[k],
        fs       = fs,
        K        = k,
        gamma_M1 = gamma_M1_vec[k],
        gamma_M2 = gamma_M2_vec[k],
        alpha    = alpha,
        approach = "formula"
      )

      result[k, ] <- list(
        K         = k,
        N         = N_vec[k],
        alpha_adj = res_k$alpha_adj,
        gamma_M1  = gamma_M1_vec[k],
        gamma_M2  = gamma_M2_vec[k],
        Power     = res_k$formula_result$Power,
        RCP_M1    = res_k$formula_result$RCP_M1,
        RCP_M2    = res_k$formula_result$RCP_M2
      )
    }

    # ========== Simulation Approach ==========
  } else {

    all_methods <- c("bonferroni", "holm", "hochberg", "hommel")

    # Initialise one result data frame per MCP method
    result <- stats::setNames(
      lapply(all_methods, function(m) .empty_result_df()),
      all_methods
    )

    for (k in seq_len(K_max)) {
      idx_k     <- seq_len(k)
      delta_k   <- delta[idx_k]
      Sigma_k   <- Sigma_use[idx_k, idx_k, drop = FALSE]
      Sigma_arg <- if (all(Sigma_k == diag(k))) NULL else Sigma_k

      # rcpmcp_single computes all four MCP methods from a single random draw
      res_k <- rcpmcp_single(
        delta    = delta_k,
        Sigma    = Sigma_arg,
        N        = N_vec[k],
        fs       = fs,
        K        = k,
        gamma_M1 = gamma_M1_vec[k],
        gamma_M2 = gamma_M2_vec[k],
        alpha    = alpha,
        approach = "simulation",
        nsim     = nsim,
        seed     = seed
      )

      # Store results for each MCP method
      for (mth in all_methods) {
        sr_row <- res_k$sim_results[res_k$sim_results$mcp_method == mth, ]
        result[[mth]][k, ] <- list(
          K         = k,
          N         = N_vec[k],
          alpha_adj = res_k$alpha_adj,
          gamma_M1  = gamma_M1_vec[k],
          gamma_M2  = gamma_M2_vec[k],
          Power     = sr_row$Power,
          RCP_M1    = sr_row$RCP_M1,
          RCP_M2    = sr_row$RCP_M2
        )
      }
    }
  }

  # ========== Output ==========
  out <- list(
    approach = approach,
    nsim     = if (approach == "simulation") nsim else NULL,
    delta    = delta,
    Sigma    = Sigma_use,
    N        = N,
    fs       = fs,
    K_max    = K_max,
    gamma_M1 = gamma_M1,
    gamma_M2 = gamma_M2,
    alpha    = alpha,
    result   = result
  )
  class(out) <- "rcpmcp_multiple"
  return(out)
}


#' Print Method for rcpmcp_multiple Objects
#'
#' @description
#' Print a summary of an \code{"rcpmcp_multiple"} object.
#'
#' For \code{approach = "formula"}, a single results table is shown
#' (Bonferroni closed-form). For \code{approach = "simulation"}, results
#' for all four multiplicity adjustment procedures are stored in
#' \code{$result} and displayed in separate blocks. The \code{mcp_method}
#' argument selects which procedure(s) to include and in what order.
#'
#' @param x An object of class \code{"rcpmcp_multiple"}.
#' @param digits A non-negative integer specifying the number of decimal
#'   places for probability values. Default is \code{4}.
#' @param mcp_method Character vector. One or more of \code{"bonferroni"},
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"} specifying which
#'   procedure(s) to display. The display order follows the order supplied.
#'   Only used when \code{x$approach = "simulation"}; ignored otherwise.
#'   Default is \code{c("bonferroni", "holm", "hochberg", "hommel")}
#'   (all four).
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @rdname rcpmcp_multiple
#' @export
print.rcpmcp_multiple <- function(x,
                                  digits     = 4,
                                  mcp_method = c("bonferroni", "holm",
                                                 "hochberg",  "hommel"),
                                  ...) {
  cat("\nRegional Consistency Probability under Multiple Comparison Procedures\n")
  cat(strrep("-", 70), "\n")

  if (x$approach == "simulation") {
    cat(sprintf("   Approach          : Simulation-Based (nsim = %d)\n",
                x$nsim))
  } else {
    cat("   Approach          : Closed-Form Solution (Bonferroni)\n")
  }

  cat(sprintf("   Treatment effect  : delta    = (%s)\n",
              paste(formatC(x$delta, format = "f", digits = digits),
                    collapse = ", ")))
  cat(sprintf("   Std. deviation    : sd       = (%s)\n",
              paste(formatC(sqrt(diag(x$Sigma)), format = "f", digits = digits),
                    collapse = ", ")))

  # Display N: scalar or vector
  if (length(x$N) == 1L) {
    cat(sprintf("   Total sample size : N        = %d\n", x$N))
  } else {
    cat(sprintf("   Total sample size : N        = (%s)\n",
                paste(x$N, collapse = ", ")))
  }

  cat(sprintf("   Regional props.   : fs       = (%s)\n",
              paste(x$fs, collapse = ", ")))
  cat(sprintf("   Max. endpoints    : K_max    = %d\n", x$K_max))

  # Display gamma_M1 / gamma_M2: scalar or vector
  if (length(x$gamma_M1) == 1L) {
    cat(sprintf("   Threshold (M1)    : gamma_M1 = %.4f\n", x$gamma_M1))
  } else {
    cat(sprintf("   Threshold (M1)    : gamma_M1 = (%s)\n",
                paste(formatC(x$gamma_M1, format = "f", digits = digits),
                      collapse = ", ")))
  }
  if (length(x$gamma_M2) == 1L) {
    cat(sprintf("   Threshold (M2)    : gamma_M2 = %.4f\n", x$gamma_M2))
  } else {
    cat(sprintf("   Threshold (M2)    : gamma_M2 = (%s)\n",
                paste(formatC(x$gamma_M2, format = "f", digits = digits),
                      collapse = ", ")))
  }

  cat(sprintf("   Significance lvl  : alpha    = %.4f\n", x$alpha))

  corr_use <- cov2cor(x$Sigma)
  is_indep <- all(abs(corr_use - diag(x$K_max)) < 1e-8)
  if (is_indep) {
    cat("   Endpoint corr.    : Independent (identity matrix)\n")
  } else {
    cat("   Endpoint corr.    : User-specified (see $Sigma)\n")
  }

  cat(strrep("-", 70), "\n")

  # ---- Helper: format and print one result data frame ----
  .print_result_df <- function(df) {
    df_fmt            <- df
    df_fmt$alpha_adj  <- formatC(df$alpha_adj, format = "f",
                                 digits = digits + 1)
    df_fmt$gamma_M1   <- formatC(df$gamma_M1,  format = "f", digits = digits)
    df_fmt$gamma_M2   <- formatC(df$gamma_M2,  format = "f", digits = digits)
    for (col in c("Power", "RCP_M1", "RCP_M2")) {
      df_fmt[[col]] <- formatC(df[[col]], format = "f", digits = digits)
    }
    colnames(df_fmt) <- c("K", "N", "alpha_adj", "gamma_M1", "gamma_M2",
                          "Power", "RCP (Method 1)", "RCP (Method 2)")
    print(df_fmt, row.names = FALSE, right = TRUE)
  }

  # ---- Formula: single block ----
  if (x$approach == "formula") {
    cat("Results:\n\n")
    .print_result_df(x$result)

    # ---- Simulation: one block per selected MCP method ----
  } else {
    valid_methods <- c("bonferroni", "holm", "hochberg", "hommel")
    mcp_method    <- mcp_method[mcp_method %in% valid_methods]
    if (length(mcp_method) == 0L) {
      warning("No valid mcp_method specified; displaying all four methods.")
      mcp_method <- valid_methods
    }

    cap <- function(s) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))

    for (mth in mcp_method) {
      cat(sprintf("Results [%s]:\n\n", cap(mth)))
      .print_result_df(x$result[[mth]])
      cat("\n")
    }
  }

  cat("\n")
  invisible(x)
}


#' @rdname rcpmcp_multiple
#'
#' @description
#' \code{plot.rcpmcp_multiple}: Visualise disjunctive power, RCP (Method 1),
#' and RCP (Method 2) as a function of the number of endpoints \eqn{K}.
#' Selected panels are arranged with \code{ggplot2::facet_wrap(scales =
#' "free_y")}. When additional \code{rcpmcp_multiple} objects are supplied
#' via \code{overlay}, all curves are overlaid on a single figure with
#' distinct colours.
#'
#' Available panels (controlled by the \code{panels} argument):
#' \itemize{
#'   \item \code{"Power"}:  Disjunctive power.
#'   \item \code{"RCP_M1"}: RCP under Method 1 (effect-retention).
#'   \item \code{"RCP_M2"}: RCP under Method 2 (all-regions-positive).
#' }
#'
#' When \code{approach = "simulation"}, the \code{mcp_method} argument
#' (scalar) selects which of the four procedures to plot. To compare
#' procedures visually, create separate \code{rcpmcp_multiple} objects and
#' use \code{overlay}.
#'
#' A horizontal reference line is drawn at \code{target_power} in the Power
#' panel when \code{"Power"} is included. The x-axis represents the number
#' of endpoints \eqn{K = 1, \ldots, K_{\max}}.
#'
#' @param overlay An optional \strong{named list} of additional
#'   \code{"rcpmcp_multiple"} objects to overlay on the same figure
#'   (e.g. results for different correlation levels). All objects must share
#'   the same \code{K_max} as \code{x}. Default is \code{NULL} (plot
#'   \code{x} alone).
#' @param group_labels Character vector or \code{expression}. Legend labels
#'   for all conditions in the order \code{c(x, overlay[[1]], overlay[[2]],
#'   \ldots)}. Supports \code{expression()} for mathematical annotation,
#'   e.g. \code{expression(rho == 0, rho == 0.5)}.
#'   Default is \code{NULL} (auto-generated labels).
#' @param panels Character vector selecting which panels to display. Any
#'   non-empty subset of \code{c("Power", "RCP_M1", "RCP_M2")}. The display
#'   order is always \code{Power -> RCP_M1 -> RCP_M2}. Default shows all
#'   three panels.
#' @param mcp_method Character scalar. Which multiplicity adjustment procedure
#'   to plot when \code{approach = "simulation"}. One of
#'   \code{"bonferroni"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"}.
#'   Ignored when \code{approach = "formula"}. Default is
#'   \code{"bonferroni"}.
#' @param target_power Numeric scalar. A horizontal reference line drawn in
#'   the Power panel when \code{"Power"} is included in \code{panels}. Set
#'   to \code{NULL} to suppress. Default is \code{0.8}.
#' @param layout Integer vector of length 2 specifying panel arrangement as
#'   \code{c(nrow, ncol)}. For example, \code{c(1, 3)} produces 1 row and
#'   3 columns, \code{c(3, 1)} (the default) produces 3 rows and 1 column.
#'   When \code{ncol > 1}, strip labels are automatically moved to the top
#'   of each panel.
#' @param fixed_rcp_ylim Logical. If \code{TRUE}, the y-axis of
#'   \code{"RCP_M1"} and \code{"RCP_M2"} panels is fixed to \eqn{[0, 1]},
#'   while the Power panel retains free scaling. Requires the \pkg{ggh4x}
#'   package. Default is \code{FALSE}.
#' @param colours Character vector of colour codes. Recycled if shorter than
#'   the number of conditions. Default uses the Okabe-Ito palette
#'   (colour-blind friendly).
#' @param point_size Numeric scalar. Size of the points. Default \code{2.5}.
#' @param line_width Numeric scalar. Width of the lines. Default \code{0.8}.
#' @param base_size Numeric scalar. Base font size for
#'   \code{ggplot2::theme_bw}. Default \code{11}.
#'
#' @return \code{plot.rcpmcp_multiple} returns a \code{ggplot} object.
#'   Save with \code{ggplot2::ggsave()}.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline
#'   facet_wrap scale_x_continuous scale_y_continuous scale_colour_manual
#'   theme_bw theme element_text element_rect element_blank margin unit
#' @importFrom rlang .data
#'
#' @examples
#' # plot(): Single object, all three panels (default)
#' res <- rcpmcp_multiple(
#'   delta = rep(0.3, 5), Sigma = diag(5),
#'   N = 200, fs = c(0.1, 0.45, 0.45), K_max = 5
#' )
#' p <- plot(res)
#' print(p)
#'
#' # plot(): RCP panels only
#' p2 <- plot(res, panels = c("RCP_M1", "RCP_M2"))
#' print(p2)
#'
#' # plot(): Overlay multiple conditions (varying correlation)
#' make_sigma <- function(rho, k) {
#'   S <- matrix(rho, k, k); diag(S) <- 1; S
#' }
#' res05 <- rcpmcp_multiple(
#'   delta = rep(0.3, 5), Sigma = make_sigma(0.5, 5),
#'   N = 200, fs = c(0.1, 0.45, 0.45), K_max = 5
#' )
#' p3 <- plot(res,
#'            overlay      = list(res05),
#'            group_labels = expression(rho == 0, rho == 0.5))
#' print(p3)
#'
#' @export
plot.rcpmcp_multiple <- function(x,
                                 ...,
                                 overlay        = NULL,
                                 group_labels   = NULL,
                                 panels         = c("Power", "RCP_M1",
                                                    "RCP_M2"),
                                 mcp_method     = "bonferroni",
                                 target_power   = 0.8,
                                 layout         = c(length(panels), 1L),
                                 fixed_rcp_ylim = FALSE,
                                 colours        = NULL,
                                 point_size     = 2.5,
                                 line_width     = 0.8,
                                 base_size      = 11) {

  # ========== Validate mcp_method ==========
  valid_methods <- c("bonferroni", "holm", "hochberg", "hommel")
  if (length(mcp_method) != 1L || !mcp_method %in% valid_methods) {
    stop('mcp_method must be a single value: "bonferroni", "holm", ',
         '"hochberg", or "hommel".')
  }

  # ========== Helper: extract result data frame ==========
  # For formula objects: return $result directly.
  # For simulation objects: return $result[[mcp_method]].
  .get_result_df <- function(obj) {
    if (obj$approach == "formula") {
      obj$result
    } else {
      obj$result[[mcp_method]]
    }
  }

  # ========== Build ordered object list: x first, then overlay ==========
  if (!is.null(overlay)) {
    if (!is.list(overlay) ||
        !all(vapply(overlay, inherits, logical(1L), "rcpmcp_multiple"))) {
      stop("overlay must be a list of 'rcpmcp_multiple' objects.")
    }
    obj_list <- c(list(x), overlay)
  } else {
    obj_list <- list(x)
  }
  n_cond       <- length(obj_list)
  single_input <- (n_cond == 1L)

  # ========== Assign internal keys and legend labels ==========
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
  valid_panels <- c("Power", "RCP_M1", "RCP_M2")
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
    stop("All rcpmcp_multiple objects must have the same K_max.")
  }
  K_max <- unique(K_vals)

  # ========== Build tidy long data frame ==========
  panel_label_map <- c(
    Power  = "Power",
    RCP_M1 = "RCP (Method 1)",
    RCP_M2 = "RCP (Method 2)"
  )
  panel_col_map <- c(
    Power  = "Power",
    RCP_M1 = "RCP_M1",
    RCP_M2 = "RCP_M2"
  )
  active_panels <- valid_panels[valid_panels %in% panels]
  panel_levels  <- unname(panel_label_map[active_panels])

  df_list <- lapply(internal_keys, function(key) {
    res <- .get_result_df(obj_list[[key]])
    do.call(rbind, lapply(active_panels, function(pnl) {
      data.frame(
        condition = key,
        K         = res$K,
        panel     = unname(panel_label_map[pnl]),
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
      strip.position = strip_pos
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

  # ========== Reference line in Power panel ==========
  if (!is.null(target_power) && "Power" %in% panels) {
    if (!is.numeric(target_power) || length(target_power) != 1L ||
        target_power <= 0 || target_power >= 1) {
      stop("target_power must be a single numeric value in (0, 1) or NULL.")
    }
    ref_df <- data.frame(
      panel = factor("Power", levels = panel_levels),
      yint  = target_power
    )
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
  # When fixed_rcp_ylim = TRUE, RCP panels are fixed to [0, 1] while
  # the Power panel retains free_y scaling.
  if (fixed_rcp_ylim) {
    rcp_labels  <- unname(panel_label_map[c("RCP_M1", "RCP_M2")])
    scales_list <- lapply(panel_levels, function(lv) {
      if (lv %in% rcp_labels) {
        ggplot2::scale_y_continuous(limits = c(0, 1))
      } else {
        ggplot2::scale_y_continuous(limits = NULL)
      }
    })
    p <- p + ggh4x::facetted_pos_scales(y = scales_list)
  }

  p
}
