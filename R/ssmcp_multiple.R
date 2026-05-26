#' Sample Size Calculation under the Bonferroni Procedure
#' (Multiple Evaluation)
#'
#' @description
#' Compute the required total sample size for \eqn{K = 1, 2, \ldots, K_{\max}}
#' endpoints by repeatedly calling \code{\link{ssmcp_single}}. This function
#' is intended to investigate how the required sample size changes as the
#' number of endpoints increases under the Bonferroni multiplicity adjustment
#' procedure.
#'
#' For each value of \eqn{K}, \code{\link{ssmcp_single}} is called using the
#' first \eqn{K} elements of \code{delta} and the upper-left
#' \eqn{K}-by-\eqn{K} submatrix of \code{Sigma}. The target is the
#' disjunctive power (probability that at least one endpoint is significant).
#'
#' Since the Bonferroni procedure is the most conservative among common
#' familywise error rate (FWER)-controlling procedures, the sample sizes
#' computed here provide conservative upper bounds that are also valid for
#' step-down procedures such as Holm, Hochberg, and Hommel. To verify the
#' power achieved under a step-down procedure at a given \eqn{N}, use
#' \code{\link{rcpmcp_multiple}} with \code{approach = "simulation"} and
#' the corresponding \code{mcp_method}.
#'
#' @param delta Numeric vector of length \code{K_max}. True treatment effects
#'   (mean differences) for each endpoint. All elements must be positive.
#' @param Sigma \code{K_max}-by-\code{K_max} numeric variance-covariance
#'   matrix for the \code{K_max} endpoints. Must be symmetric and positive
#'   definite. \code{NULL} (default) sets \code{Sigma = diag(K_max)}, i.e.,
#'   unit variance and independence.
#' @param fs Numeric vector. Proportion of patients in each region
#'   \eqn{s = 1, \ldots, S}. Must sum to 1 and all elements must be positive.
#'   Default is \code{c(0.1, 0.45, 0.45)}.
#' @param K_max Positive integer. Maximum number of endpoints. Results are
#'   computed for \eqn{K = 1, 2, \ldots, K_{\max}}. Must match the length of
#'   \code{delta} and the dimension of \code{Sigma}. Default is \code{5}.
#' @param alpha Numeric scalar. One-sided familywise significance level before
#'   multiplicity adjustment. Default is \code{0.025}.
#' @param target_power Numeric scalar. Target disjunctive power. Must be in
#'   \eqn{(0, 1)}. Default is \code{0.8}.
#' @param N_min Positive integer. Minimum sample size passed to
#'   \code{\link{ssmcp_single}}. Default is \code{10}.
#' @param N_max Positive integer. Maximum sample size passed to
#'   \code{\link{ssmcp_single}}. Default is \code{1000000}.
#' @param tol Numeric scalar. Tolerance passed to \code{stats::uniroot}.
#'   Default is \code{1e-6}.
#'
#' @return An object of class \code{"ssmcp_multiple"}, which is a list
#'   containing:
#' \describe{
#'   \item{\code{delta}}{True treatment effects.}
#'   \item{\code{Sigma}}{Variance-covariance matrix used.}
#'   \item{\code{fs}}{Regional patient proportions.}
#'   \item{\code{K_max}}{Maximum number of endpoints evaluated.}
#'   \item{\code{alpha}}{One-sided familywise significance level.}
#'   \item{\code{target_power}}{Target disjunctive power.}
#'   \item{\code{result}}{A data frame with \code{K_max} rows and columns:
#'     \code{k}, \code{alpha_adj}, \code{N}, \code{power_achieved}.}
#' }
#'
#' @seealso \code{\link{ssmcp_single}}, \code{\link{rcpmcp_multiple}}
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats pnorm qnorm uniroot cov2cor
#' @importFrom utils combn
#'
#' @examples
#' \donttest{
#' # Example 1: K_max = 5, independent endpoints
#' ss_mult <- ssmcp_multiple(
#'   delta        = rep(0.2, 5),
#'   Sigma        = diag(5),
#'   fs           = c(0.1, 0.45, 0.45),
#'   K_max        = 5,
#'   alpha        = 0.025,
#'   target_power = 0.8
#' )
#' print(ss_mult)
#' }
#'
#' \donttest{
#' # Example 2: K_max = 5, uniform correlation rho = 0.5
#' rho              <- 0.5
#' Sigma_corr       <- matrix(rho, nrow = 5, ncol = 5)
#' diag(Sigma_corr) <- 1
#' ss_corr <- ssmcp_multiple(
#'   delta        = rep(0.2, 5),
#'   Sigma        = Sigma_corr,
#'   fs           = c(0.1, 0.45, 0.45),
#'   K_max        = 5,
#'   alpha        = 0.025,
#'   target_power = 0.8
#' )
#' print(ss_corr)
#' }
#'
#' @export
ssmcp_multiple <- function(delta,
                           Sigma        = NULL,
                           fs           = c(0.1, 0.45, 0.45),
                           K_max        = 5,
                           alpha        = 0.025,
                           target_power = 0.8,
                           N_min        = 10,
                           N_max        = 1000000,
                           tol          = 1e-6) {

  # ========== Input Validation ==========
  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a numeric vector of positive values")
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
  if (!is.numeric(fs) || any(fs <= 0) || abs(sum(fs) - 1) > 1e-8) {
    stop("fs must be a numeric vector of positive values summing to 1")
  }
  if (!is.numeric(K_max) || length(K_max) != 1 ||
      K_max < 1 || K_max != as.integer(K_max)) {
    stop("K_max must be a single positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(target_power) || length(target_power) != 1 ||
      target_power <= 0 || target_power >= 1) {
    stop("target_power must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(N_min) || length(N_min) != 1 ||
      N_min < 1 || N_min != as.integer(N_min)) {
    stop("N_min must be a single positive integer")
  }
  if (!is.numeric(N_max) || length(N_max) != 1 ||
      N_max <= N_min || N_max != as.integer(N_max)) {
    stop("N_max must be a single positive integer greater than N_min")
  }
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("tol must be a single positive numeric value")
  }

  # Use Sigma = diag(K_max) if NULL
  Sigma_use <- if (is.null(Sigma)) diag(K_max) else Sigma

  # ========== Loop over k = 1, ..., K_max ==========
  result <- data.frame(
    K             = integer(K_max),
    alpha_adj     = numeric(K_max),
    N             = integer(K_max),
    power_achieved = numeric(K_max)
  )

  for (k in seq_len(K_max)) {
    delta_k   <- delta[seq_len(k)]
    Sigma_k   <- Sigma_use[seq_len(k), seq_len(k), drop = FALSE]
    Sigma_arg <- if (all(Sigma_k == diag(k))) NULL else Sigma_k

    res_k <- ssmcp_single(
      delta        = delta_k,
      Sigma        = Sigma_arg,
      fs           = fs,
      K            = k,
      alpha        = alpha,
      target_power = target_power,
      N_min        = N_min,
      N_max        = N_max,
      tol          = tol
    )

    result[k, ] <- list(
      K             = k,
      alpha_adj     = res_k$alpha_adj,
      N             = res_k$N,
      power_achieved = res_k$power_achieved
    )
  }

  # ========== Output ==========
  out <- list(
    delta        = delta,
    Sigma        = Sigma_use,
    fs           = fs,
    K_max        = K_max,
    alpha        = alpha,
    target_power = target_power,
    result       = result
  )
  class(out) <- "ssmcp_multiple"
  return(out)
}


#' Print Method for ssmcp_multiple Objects
#'
#' @param x An object of class \code{"ssmcp_multiple"}.
#' @param digits A non-negative integer specifying the number of decimal
#'   places for probability values. Default is \code{4}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @rdname ssmcp_multiple
#' @export
print.ssmcp_multiple <- function(x, digits = 4, ...) {
  cat("\nSample Size under Bonferroni Procedure (Multiple Evaluation)\n")
  cat(strrep("-", 70), "\n")

  cat("   MCP method        : Bonferroni\n")
  cat("   Power definition  : Disjunctive (>= 1 endpoint significant)\n")
  cat(sprintf("   Treatment effect  : delta    = (%s)\n",
              paste(formatC(x$delta, format = "f", digits = digits),
                    collapse = ", ")))
  cat(sprintf("   Std. deviation    : sd       = (%s)\n",
              paste(formatC(sqrt(diag(x$Sigma)), format = "f", digits = digits),
                    collapse = ", ")))
  cat(sprintf("   Regional props.   : fs       = (%s)\n",
              paste(x$fs, collapse = ", ")))
  cat(sprintf("   Max. endpoints    : K_max    = %d\n",   x$K_max))
  cat(sprintf("   Significance lvl  : alpha    = %.4f\n", x$alpha))
  cat(sprintf("   Target power      : 1-beta   = %.4f\n", x$target_power))

  corr_use <- cov2cor(x$Sigma)
  is_indep <- all(abs(corr_use - diag(x$K_max)) < 1e-8)
  if (is_indep) {
    cat("   Endpoint corr.    : Independent (identity matrix)\n")
  } else {
    cat("   Endpoint corr.    : User-specified (see $Sigma)\n")
  }

  cat(strrep("-", 70), "\n")
  cat("Results:\n\n")

  df <- x$result
  df$alpha_adj      <- formatC(df$alpha_adj,      format = "f",
                               digits = digits + 1)
  df$power_achieved <- formatC(df$power_achieved, format = "f",
                               digits = digits)
  colnames(df) <- c("K", "alpha_adj", "N", "Power achieved")
  print(df, row.names = FALSE, right = TRUE)

  cat("\n")
  invisible(x)
}
