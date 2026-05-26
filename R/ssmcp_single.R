#' Sample Size Calculation under the Bonferroni Procedure
#' (Single Evaluation)
#'
#' @description
#' Calculate the required total sample size \eqn{N} for a multi-regional
#' clinical trial (MRCT) with \code{K} endpoints such that the disjunctive
#' power (probability that at least one endpoint achieves statistical
#' significance) meets or exceeds \code{target_power} under the Bonferroni
#' multiplicity adjustment procedure.
#'
#' Under the Bonferroni procedure, each endpoint is tested at the adjusted
#' significance level \eqn{\alpha_K = \alpha / K}. The disjunctive power is
#' computed via the inclusion-exclusion principle using
#' \code{mvtnorm::pmvnorm}, following the same approach as
#' \code{\link{rcpmcp_single}} with \code{approach = "formula"}.
#'
#' Sample size determination for step-down procedures such as Holm,
#' Hochberg, and Hommel is not supported because these methods are data
#' dependent and do not admit a fixed per-endpoint significance level prior
#' to observing data. Since the Bonferroni procedure is the most conservative
#' among common familywise error rate (FWER)-controlling procedures, the
#' sample size obtained here provides a conservative upper bound that is also
#' valid for Holm, Hochberg, and Hommel. Users wishing to verify the
#' achieved power under a step-down procedure for a given \eqn{N} can do so
#' via \code{\link{rcpmcp_single}} with \code{approach = "simulation"} and
#' the corresponding \code{mcp_method}.
#'
#' The per-endpoint closed-form solution
#' \eqn{N_0 = \{(z_{1-\alpha/K} + z_{1-\beta}) \sigma_j / \delta_j\}^2}
#' is used as the initial value for \code{stats::uniroot}, taking the maximum
#' over all endpoints.
#'
#' @param delta Numeric vector of length \code{k}. True treatment effects
#'   (mean differences) for each endpoint. All elements must be positive.
#' @param Sigma \code{k}-by-\code{k} numeric variance-covariance matrix for
#'   the \code{K} endpoints. Must be symmetric and positive definite, with
#'   diagonal elements representing the variance of each endpoint's treatment
#'   effect estimate (i.e., \eqn{\sigma_j^2}). Off-diagonal elements represent
#'   covariances between endpoints. \code{NULL} (default) sets
#'   \code{Sigma = diag(k)}, i.e., unit variance and independence.
#' @param fs Numeric vector. Proportion of patients in each region
#'   \eqn{s = 1, \ldots, S}. Must sum to 1 and all elements must be positive.
#'   Default is \code{c(0.1, 0.45, 0.45)}.
#' @param K Positive integer. Number of endpoints. Must match the length of
#'   \code{delta} and the dimension of \code{Sigma}. Default is \code{1}.
#' @param alpha Numeric scalar. One-sided familywise significance level before
#'   multiplicity adjustment. Default is \code{0.025}.
#' @param target_power Numeric scalar. Target disjunctive power. Must be in
#'   \eqn{(0, 1)}. Default is \code{0.8}.
#' @param N_min Positive integer. Minimum sample size for the
#'   \code{stats::uniroot} search. Default is \code{10}.
#' @param N_max Positive integer. Maximum sample size for the
#'   \code{stats::uniroot} search. Default is \code{1000000}.
#' @param tol Numeric scalar. Tolerance passed to \code{stats::uniroot}.
#'   Default is \code{1e-6}.
#'
#' @return An object of class \code{"ssmcp_single"}, which is a list
#'   containing:
#' \describe{
#'   \item{\code{delta}}{True treatment effects.}
#'   \item{\code{Sigma}}{Variance-covariance matrix used.}
#'   \item{\code{fs}}{Regional patient proportions.}
#'   \item{\code{K}}{Number of endpoints (\eqn{K}).}
#'   \item{\code{alpha}}{One-sided familywise significance level.}
#'   \item{\code{alpha_adj}}{Bonferroni-adjusted significance level
#'     \eqn{\alpha / k}.}
#'   \item{\code{target_power}}{Target disjunctive power.}
#'   \item{\code{N}}{Required total sample size (rounded up to integer).}
#'   \item{\code{power_achieved}}{Disjunctive power achieved at \code{N}.}
#' }
#'
#' @seealso \code{\link{rcpmcp_single}}, \code{\link{ssmcp_multiple}}
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats pnorm qnorm uniroot cov2cor
#' @importFrom utils combn
#'
#' @examples
#' # Example 1: K = 1 (single endpoint, reduces to standard formula)
#' ss1 <- ssmcp_single(
#'   delta        = 0.2,
#'   Sigma        = as.matrix(1),
#'   fs           = c(0.1, 0.45, 0.45),
#'   K            = 1,
#'   alpha        = 0.025,
#'   target_power = 0.8
#' )
#' print(ss1)
#'
#' \donttest{
#' # Example 2: K = 3, independent endpoints
#' ss3 <- ssmcp_single(
#'   delta        = c(0.2, 0.2, 0.2),
#'   Sigma        = diag(3),
#'   fs           = c(0.1, 0.45, 0.45),
#'   K            = 3,
#'   alpha        = 0.025,
#'   target_power = 0.8
#' )
#' print(ss3)
#' }
#'
#' @export
ssmcp_single <- function(delta,
                         Sigma        = NULL,
                         fs           = c(0.1, 0.45, 0.45),
                         K            = 1,
                         alpha        = 0.025,
                         target_power = 0.8,
                         N_min        = 10,
                         N_max        = 1000000,
                         tol          = 1e-6) {

  # ========== Input Validation ==========
  if (!is.numeric(delta) || any(delta <= 0)) {
    stop("delta must be a numeric vector of positive values")
  }
  if (length(delta) != K) {
    stop("length of delta must equal K")
  }
  if (!is.null(Sigma)) {
    if (!is.matrix(Sigma) || nrow(Sigma) != K || ncol(Sigma) != K) {
      stop("Sigma must be a K-by-K numeric matrix")
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
  if (!is.numeric(K) || length(K) != 1 || K < 1 || K != as.integer(K)) {
    stop("K must be a single positive integer")
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

  # ========== Common Setup ==========
  Sigma_use <- if (is.null(Sigma)) diag(K) else Sigma
  sd_j      <- sqrt(diag(Sigma_use))
  R_Z       <- cov2cor(Sigma_use)
  alpha_adj <- alpha / K
  z_crit    <- stats::qnorm(1 - alpha_adj)

  # ------------------------------------------------------------------
  # Internal: compute Bonferroni disjunctive power for a given N
  # Power = Pr(union_j {Z_j > z_{1-alpha/K}})
  # via inclusion-exclusion over all 2^k - 1 non-empty subsets
  # ------------------------------------------------------------------
  .disjunctive_power <- function(N_val) {
    lambda <- delta / (sd_j / sqrt(N_val))   # non-centrality parameters

    all_subsets <- unlist(
      lapply(seq_len(K), function(m) combn(seq_len(K), m, simplify = FALSE)),
      recursive = FALSE
    )

    pw <- 0
    for (sub in all_subsets) {
      m    <- length(sub)
      sign <- (-1)^(m + 1)
      if (m == 1) {
        pw <- pw + sign * stats::pnorm(lambda[sub] - z_crit)
      } else {
        pw <- pw + sign * mvtnorm::pmvnorm(
          lower = rep(z_crit, m),
          upper = rep(Inf, m),
          mean  = lambda[sub],
          corr  = R_Z[sub, sub, drop = FALSE],
          seed  = 1L
        )[1]
      }
    }
    pw
  }

  # ------------------------------------------------------------------
  # Initial value for uniroot
  # Per-endpoint closed-form (most conservative endpoint):
  #   N_init = max_j { ((z_crit + z_{1-beta}) * sd_j / delta_j)^2 }
  # ------------------------------------------------------------------
  z_beta <- stats::qnorm(target_power)
  N_init <- max(((z_crit + z_beta) * sd_j / delta)^2)
  N_init <- max(N_min + 1L, min(N_max - 1L, ceiling(N_init)))

  # Objective: disjunctive_power(N) - target_power = 0
  .obj <- function(N_val) .disjunctive_power(N_val) - target_power

  # Expand upper bracket until objective turns non-negative
  upper <- N_init
  for (i in seq_len(60)) {
    if (.obj(upper) >= 0) break
    upper <- min(N_max, upper * 2L)
  }

  if (.obj(N_max) < 0) {
    warning(sprintf(
      "Target power %.4f cannot be achieved within N_max = %d.",
      target_power, N_max
    ))
    N_sol <- N_max
  } else if (.obj(N_min) >= 0) {
    N_sol <- N_min
  } else {
    N_sol <- ceiling(
      stats::uniroot(.obj, lower = N_min, upper = upper, tol = tol)$root
    )
  }

  power_achieved <- .disjunctive_power(N_sol)

  # ========== Output ==========
  out <- list(
    delta          = delta,
    Sigma          = Sigma_use,
    fs             = fs,
    K              = K,
    alpha          = alpha,
    alpha_adj      = alpha_adj,
    target_power   = target_power,
    N              = N_sol,
    power_achieved = power_achieved
  )
  class(out) <- "ssmcp_single"
  return(out)
}


#' Print Method for ssmcp_single Objects
#'
#' @param x An object of class \code{"ssmcp_single"}.
#' @param digits A non-negative integer specifying the number of decimal
#'   places for probability values. Default is \code{4}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @rdname ssmcp_single
#' @export
print.ssmcp_single <- function(x, digits = 4, ...) {
  cat("\nSample Size under Bonferroni Procedure (Single Evaluation)\n")
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
  cat(sprintf("   No. of endpoints  : K        = %d\n",   x$K))
  cat(sprintf("   Significance lvl  : alpha    = %.4f\n", x$alpha))
  cat(sprintf("   Adjusted level    : alpha/K  = %.4f\n", x$alpha_adj))
  cat(sprintf("   Target power      : 1-beta   = %.4f\n", x$target_power))

  corr_use <- cov2cor(x$Sigma)
  is_indep <- all(abs(corr_use - diag(x$K)) < 1e-8)
  if (is_indep) {
    cat("   Endpoint corr.    : Independent (identity matrix)\n")
  } else {
    cat("   Endpoint corr.    : User-specified (see $Sigma)\n")
  }

  cat(strrep("-", 70), "\n")
  cat("Results:\n\n")
  cat(sprintf("   Required N        : %d\n", x$N))
  cat(sprintf("   Power achieved    : %s\n",
              formatC(x$power_achieved, format = "f", digits = digits)))
  cat("\n")

  invisible(x)
}
