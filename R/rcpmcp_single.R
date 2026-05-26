# =============================================================================
# Internal helper: pre-compute quantities shared across formula calls
# =============================================================================

#' Pre-compute Shared Quantities for Formula-Based RCP Calculations
#'
#' @description
#' An internal helper that computes and caches all quantities derived from
#' \code{delta}, \code{Sigma}, \code{N}, \code{fs}, \code{k}, and
#' \code{alpha} that do not depend on \code{gamma_M1} or \code{gamma_M2}.
#' The returned list is consumed by \code{\link{rcpmcp_single}} (formula
#' branch) and \code{\link{.rcpmcp_formula_rcp_only}}.
#'
#' Separating pre-computation from threshold-dependent calculation allows
#' \code{\link{rcpmcp_get_gamma}} to call \code{.rcpmcp_formula_rcp_only}
#' repeatedly (once per \code{uniroot} iteration) without redundant
#' recomputation of fixed quantities.
#'
#' @param delta     Numeric vector of length \code{K}. True treatment effects.
#' @param Sigma_use \code{K}-by-\code{K} variance-covariance matrix
#'   (already resolved from user input; never \code{NULL}).
#' @param N         Positive integer. Total sample size.
#' @param fs        Numeric vector summing to 1. Regional patient proportions.
#' @param K         Positive integer. Number of endpoints (\eqn{K}).
#' @param alpha     Numeric scalar. One-sided familywise significance level.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{S}}{Number of regions.}
#'   \item{\code{Ns}}{Integer vector of regional sample sizes (sum = N).}
#'   \item{\code{f1}}{Fraction of patients in region 1.}
#'   \item{\code{alpha_adj}}{Bonferroni-adjusted level \eqn{\alpha/K}.}
#'   \item{\code{z_crit}}{Critical value \eqn{z_{1-\alpha/K}}.}
#'   \item{\code{sd_j}}{Per-endpoint standard deviations (length k).}
#'   \item{\code{se_overall}}{Overall standard errors (length k).}
#'   \item{\code{se_s}}{\code{K}-by-\code{S} matrix of regional SE.}
#'   \item{\code{lambda}}{Non-centrality parameters (length K).}
#'   \item{\code{R_Z}}{Correlation matrix of \eqn{Z_j}'s
#'     (\code{K}-by-\code{K}).}
#'   \item{\code{all_subsets}}{All \eqn{2^K - 1} non-empty subsets of
#'     \eqn{\{1,\ldots,K\}}.}
#'   \item{\code{sqrt_fs}}{Pre-computed \eqn{\sqrt{f_s}} (length S).}
#' }
#'
#' @keywords internal
rcpmcp_single_precompute <- function(delta, Sigma_use, N, fs, K, alpha) {
  S         <- length(fs)
  alpha_adj <- alpha / K
  z_crit    <- stats::qnorm(1 - alpha_adj)
  sd_j      <- sqrt(diag(Sigma_use))

  # Regional sample sizes: use floor() for all regions, then assign the
  # remainder to the last region so that sum(Ns) == N exactly.
  Ns    <- floor(fs * N)
  Ns[S] <- N - sum(Ns[-S])

  f1         <- fs[1]
  se_overall <- sd_j / sqrt(N)
  se_s       <- outer(sd_j, sqrt(Ns), "/")   # k-by-S
  lambda     <- delta / se_overall
  R_Z        <- cov2cor(Sigma_use)

  # Pre-compute all 2^K - 1 non-empty subsets (fixed for a given K)
  all_subsets <- unlist(
    lapply(seq_len(K), function(m) combn(seq_len(K), m, simplify = FALSE)),
    recursive = FALSE
  )

  list(
    S           = S,
    Ns          = Ns,
    f1          = f1,
    alpha_adj   = alpha_adj,
    z_crit      = z_crit,
    sd_j        = sd_j,
    se_overall  = se_overall,
    se_s        = se_s,
    lambda      = lambda,
    R_Z         = R_Z,
    all_subsets = all_subsets,
    sqrt_fs     = sqrt(fs)
  )
}


# =============================================================================
# Internal helper: formula RCP for M1 and M2 only (Power skipped)
# =============================================================================

#' Compute Formula-Based RCP (Method 1 and Method 2) Without Power
#'
#' @description
#' An internal function used by \code{\link{rcpmcp_get_gamma}} during
#' \code{uniroot} iterations. Given pre-computed fixed quantities (from
#' \code{\link{rcpmcp_single_precompute}}) and threshold values
#' \code{gamma_M1} / \code{gamma_M2}, it computes only the numerators of
#' RCP (Method 1) and RCP (Method 2) via the inclusion-exclusion principle,
#' without recomputing the Power denominator.
#'
#' The caller (\code{rcpmcp_get_gamma}) supplies a cached \code{power_out}
#' value so that \code{RCP_M1 = numer_m1 / power_out} and
#' \code{RCP_M2 = numer_m2 / power_out} can be returned.
#'
#' @param pre       List returned by \code{\link{rcpmcp_single_precompute}}.
#' @param delta     Numeric vector of length \code{K}. True treatment effects.
#' @param Sigma_use \code{K}-by-\code{K} variance-covariance matrix.
#' @param gamma_M1  Numeric scalar. Effect retention threshold for Method 1.
#' @param gamma_M2  Numeric scalar. Consistency threshold for Method 2.
#' @param power_out Numeric scalar. Pre-computed Power (denominator for RCP).
#'
#' @return A named list with \code{RCP_M1} and \code{RCP_M2}.
#'
#' @keywords internal
.rcpmcp_formula_rcp_only <- function(pre, delta, Sigma_use,
                                     gamma_M1, gamma_M2, power_out) {
  S           <- pre$S
  Ns          <- pre$Ns
  f1          <- pre$f1
  z_crit      <- pre$z_crit
  sd_j        <- pre$sd_j
  se_overall  <- pre$se_overall
  se_s        <- pre$se_s
  lambda      <- pre$lambda
  R_Z         <- pre$R_Z
  all_subsets <- pre$all_subsets
  sqrt_fs     <- pre$sqrt_fs
  K           <- length(delta)

  # ---------- Method 1 pre-computation (depends on gamma_M1) ----------
  # D_j = hat_delta_{1j} - gamma_M1 * hat_delta_j
  # Var(D_j) = Sigma[j,j] * c_DD
  # c_DD = (1 - 2*gamma_M1*f1) / Ns[1] + gamma_M1^2 / N
  N     <- sum(Ns)   # recover N (exact, since sum(Ns) == N by construction)
  c_DD  <- (1 - 2 * gamma_M1 * f1) / Ns[1] + gamma_M1^2 / N
  var_d <- diag(Sigma_use) * c_DD
  sd_d  <- sqrt(var_d)
  mu_V  <- (1 - gamma_M1) * delta / sd_d

  # Corr(V_j, Z_j) within endpoint j
  rho_VZ_diag <- diag(Sigma_use) * (f1 / Ns[1] - gamma_M1 / N) /
    (se_overall * sd_d)

  # Common cross-endpoint factor (depends on gamma_M1)
  cross_factor <- f1 / Ns[1] - gamma_M1 / N

  # Helper: Pr(all A_j^{M1} for j in idx)
  .prob_m1_subset <- function(idx) {
    m   <- length(idx)
    dim <- 2L * m

    lower_vec <- numeric(dim)
    mean_vec  <- numeric(dim)
    for (ii in seq_len(m)) {
      j <- idx[ii]
      lower_vec[2L * ii - 1L] <- 0
      lower_vec[2L * ii]      <- z_crit
      mean_vec[2L * ii - 1L]  <- mu_V[j]
      mean_vec[2L * ii]       <- lambda[j]
    }

    corr_mat <- diag(dim)
    for (ii in seq_len(m)) {
      j <- idx[ii]
      corr_mat[2L * ii - 1L, 2L * ii] <- rho_VZ_diag[j]
      corr_mat[2L * ii, 2L * ii - 1L] <- rho_VZ_diag[j]

      for (jj in seq_len(m)) {
        if (ii == jj) next
        l      <- idx[jj]
        sig_jl <- Sigma_use[j, l]

        corr_mat[2L * ii - 1L, 2L * jj - 1L] <-
          sig_jl * c_DD / (sd_d[j] * sd_d[l])
        corr_mat[2L * ii - 1L, 2L * jj] <-
          sig_jl * cross_factor / (sd_d[j] * se_overall[l])
        corr_mat[2L * ii, 2L * jj - 1L] <-
          sig_jl * cross_factor / (se_overall[j] * sd_d[l])
        corr_mat[2L * ii, 2L * jj] <- R_Z[j, l]
      }
    }

    mvtnorm::pmvnorm(
      lower = lower_vec,
      upper = rep(Inf, dim),
      mean  = mean_vec,
      corr  = corr_mat,
      seed  = 1L
    )[1]
  }

  # Helper: Pr(all A_j^{M2} for j in idx) (depends on gamma_M2)
  .prob_m2_subset <- function(idx) {
    m   <- length(idx)
    blk <- S + 1L
    dim <- blk * m

    lower_vec <- numeric(dim)
    mean_vec  <- numeric(dim)
    for (ii in seq_len(m)) {
      j    <- idx[ii]
      base <- (ii - 1L) * blk
      lower_vec[base + 1L] <- z_crit
      mean_vec[base + 1L]  <- lambda[j]
      for (s in seq_len(S)) {
        lower_vec[base + 1L + s] <- gamma_M2 / se_s[j, s]
        mean_vec[base + 1L + s]  <- delta[j] / se_s[j, s]
      }
    }

    corr_mat <- diag(dim)
    for (ii in seq_len(m)) {
      j       <- idx[ii]
      base_ii <- (ii - 1L) * blk
      for (s in seq_len(S)) {
        corr_mat[base_ii + 1L, base_ii + 1L + s] <- sqrt_fs[s]
        corr_mat[base_ii + 1L + s, base_ii + 1L] <- sqrt_fs[s]
      }
      for (jj in seq_len(m)) {
        if (ii == jj) next
        l       <- idx[jj]
        base_jj <- (jj - 1L) * blk
        r_jl    <- R_Z[j, l]
        corr_mat[base_ii + 1L, base_jj + 1L] <- r_jl
        for (s in seq_len(S)) {
          corr_mat[base_ii + 1L,     base_jj + 1L + s] <- sqrt_fs[s] * r_jl
          corr_mat[base_ii + 1L + s, base_jj + 1L]     <- sqrt_fs[s] * r_jl
          corr_mat[base_ii + 1L + s, base_jj + 1L + s] <- r_jl
        }
      }
    }

    mvtnorm::pmvnorm(
      lower = lower_vec,
      upper = rep(Inf, dim),
      mean  = mean_vec,
      corr  = corr_mat,
      seed  = 1L
    )[1]
  }

  # Inclusion-exclusion for M1 and M2 numerators only
  numer_m1 <- 0
  numer_m2 <- 0
  for (sub in all_subsets) {
    sign     <- (-1)^(length(sub) + 1L)
    numer_m1 <- numer_m1 + sign * .prob_m1_subset(sub)
    numer_m2 <- numer_m2 + sign * .prob_m2_subset(sub)
  }

  list(
    RCP_M1 = if (power_out > 0) numer_m1 / power_out else NA_real_,
    RCP_M2 = if (power_out > 0) numer_m2 / power_out else NA_real_
  )
}


# =============================================================================
# Internal helper: formula Power only
# =============================================================================

#' Compute Formula-Based Disjunctive Power
#'
#' @description
#' An internal function that computes the Bonferroni disjunctive power
#' using pre-computed quantities from
#' \code{\link{rcpmcp_single_precompute}}. Used by
#' \code{\link{rcpmcp_get_gamma}} to cache Power once per \eqn{K}, avoiding
#' recomputation during \code{uniroot} iterations.
#'
#' @param pre   List returned by \code{\link{rcpmcp_single_precompute}}.
#'
#' @return Numeric scalar. Disjunctive power.
#'
#' @keywords internal
.rcpmcp_formula_power <- function(pre) {
  z_crit      <- pre$z_crit
  lambda      <- pre$lambda
  R_Z         <- pre$R_Z
  all_subsets <- pre$all_subsets

  .prob_power_subset <- function(idx) {
    m <- length(idx)
    if (m == 1L) {
      return(stats::pnorm(lambda[idx] - z_crit))
    }
    mvtnorm::pmvnorm(
      lower = rep(z_crit, m),
      upper = rep(Inf, m),
      mean  = lambda[idx],
      corr  = R_Z[idx, idx, drop = FALSE],
      seed  = 1L
    )[1]
  }

  pw <- 0
  for (sub in all_subsets) {
    pw <- pw + (-1)^(length(sub) + 1L) * .prob_power_subset(sub)
  }
  pw
}


# =============================================================================
# Main exported function
# =============================================================================

#' Regional Consistency Probability under Multiple Comparison Procedures
#' (Single Evaluation)
#'
#' @description
#' Calculate the power and regional consistency probabilities (RCPs) for
#' multi-regional clinical trials (MRCTs) with \code{K} endpoints subject
#' to a multiplicity adjustment procedure.
#'
#' \strong{Closed-form solution} (\code{approach = "formula"}):
#' Power and RCP are computed under the Bonferroni procedure, with the
#' adjusted significance level \eqn{\alpha_K = \alpha / K}, using the
#' inclusion-exclusion principle over all \eqn{2^K - 1} non-empty subsets
#' of endpoints. Each term is evaluated via multivariate normal integration
#' (\code{mvtnorm::pmvnorm}). The endpoint covariance structure specified in
#' \code{Sigma} is fully accounted for in all calculations. All integration
#' variables are standardized so that the covariance matrix reduces to a
#' correlation matrix, following the approach in Homma (2021).
#'
#' \strong{Monte Carlo simulation} (\code{approach = "simulation"}):
#' Power and RCP are estimated via Monte Carlo simulation using a single set
#' of random draws shared across all four multiplicity adjustment procedures
#' (\code{"bonferroni"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"}).
#' Results for all four procedures are stored simultaneously in
#' \code{$sim_results}. All four methods control the familywise error rate
#' (FWER). Adjusted p-values are computed via \code{stats::p.adjust}.
#'
#' RCP is defined as:
#' \deqn{
#'   \mathrm{RCP}^{(m)}(K) = \frac{
#'     \Pr\!\left(\bigcup_{j=1}^{K}
#'       \bigl\{ C_j^{(m)} \cap \{ Z_j > z_{1-\alpha/K} \} \bigr\}
#'     \right)
#'   }{
#'     \Pr\!\left(\bigcup_{j=1}^{K}
#'       \{ Z_j > z_{1-\alpha/K} \}
#'     \right)
#'   }
#' }
#' where \eqn{C_j^{(m)}} denotes the regional consistency event for
#' endpoint \eqn{j} under Method \eqn{m} (m = 1 or 2).
#'
#' Two consistency evaluation methods follow the Japanese MHLW (2007)
#' guidance:
#' \itemize{
#'   \item \strong{Method 1}: Region 1 retains at least a fraction
#'     \code{gamma_M1} of the overall treatment effect:
#'     \eqn{\hat{\delta}_{1j} > \gamma_{M1} \hat{\delta}_j}.
#'   \item \strong{Method 2}: All regional estimates exceed \code{gamma_M2}:
#'     \eqn{\hat{\delta}_{sj} > \gamma_{M2}} for all \eqn{s = 1, \ldots, S}.
#' }
#'
#' @param delta Numeric vector of length \code{K} or a numeric matrix of
#'   dimension \code{K}-by-\code{S}. When a vector is supplied, the same
#'   treatment effect is assumed for all regions, that is,
#'   \eqn{\delta_{k,s} = \delta_{k}} for all \eqn{s}. When a matrix is
#'   supplied, element \code{delta[k, s]} represents the true treatment
#'   effect for endpoint \eqn{k} in region \eqn{s}, allowing
#'   region-specific treatment effects. Region-specific treatment effects
#'   are only supported under \code{approach = "simulation"}. All elements
#'   must be non-negative. Zero values correspond to the null hypothesis
#'   and can be used to evaluate type I error.
#' @param Sigma \code{K}-by-\code{K} numeric variance-covariance matrix for
#'   the \code{K} endpoints. Must be symmetric and positive definite, with
#'   diagonal elements representing the variance of each endpoint's treatment
#'   effect estimate (i.e., \eqn{\sigma_j^2}). Off-diagonal elements represent
#'   covariances between endpoints. \code{NULL} (default) sets
#'   \code{Sigma = diag(k)}, i.e., unit variance and independence.
#' @param N Positive integer. Total sample size across all regions.
#' @param fs Numeric vector. Proportion of patients in each region
#'   \eqn{s = 1, \ldots, S}. Must sum to 1 and all elements must be positive.
#'   Default is \code{c(0.1, 0.45, 0.45)}.
#' @param K Positive integer. Number of endpoints. Must match the length of
#'   \code{delta} and the dimension of \code{Sigma}. Default is \code{1}.
#' @param gamma_M1 Numeric scalar. Effect retention threshold for Method 1.
#'   Region 1 must satisfy
#'   \eqn{\hat{\delta}_{1j} > \gamma_{M1} \hat{\delta}_j}.
#'   Must be in \eqn{[0, 1]}. Default is \code{0.5}.
#' @param gamma_M2 Numeric scalar. Consistency threshold for Method 2.
#'   All regional estimates must satisfy
#'   \eqn{\hat{\delta}_{sj} > \gamma_{M2}}.
#'   Must be non-negative. Default is \code{0}.
#' @param alpha Numeric scalar. One-sided familywise significance level before
#'   multiplicity adjustment. Default is \code{0.025}.
#' @param approach Character scalar. Calculation approach: \code{"formula"}
#'   for the closed-form Bonferroni solution, or \code{"simulation"} for
#'   Monte Carlo simulation. Default is \code{"formula"}.
#' @param nsim Positive integer. Number of Monte Carlo iterations. Used only
#'   when \code{approach = "simulation"}. Default is \code{10000}.
#' @param seed Non-negative integer. Random seed for reproducibility. Used
#'   only when \code{approach = "simulation"}. Default is \code{1}.
#' @param variance_known Logical scalar. If \code{TRUE} (default), the
#'   common variance is treated as known and Wald-type Z-statistics with
#'   normal critical values are used. If \code{FALSE}, the variance is
#'   treated as unknown and Welch-pool t-statistics with degrees of freedom
#'   \eqn{\nu = N - 2} are used; the sample variance-covariance matrix is
#'   sampled from a Wishart distribution to avoid generating
#'   individual-level data. Only supported under
#'   \code{approach = "simulation"}.
#'
#' @return An object of class \code{"rcpmcp_single"}, which is a list
#'   containing:
#' \describe{
#'   \item{\code{approach}}{Calculation approach used.}
#'   \item{\code{nsim}}{Number of Monte Carlo iterations (\code{NULL} for
#'     the \code{"formula"} approach).}
#'   \item{\code{delta}}{True treatment effects.}
#'   \item{\code{Sigma}}{Variance-covariance matrix used.}
#'   \item{\code{N}}{Total sample size.}
#'   \item{\code{fs}}{Regional patient proportions.}
#'   \item{\code{k}}{Number of endpoints (\eqn{K}).}
#'   \item{\code{gamma_M1}}{Effect retention threshold for Method 1.}
#'   \item{\code{gamma_M2}}{Consistency threshold for Method 2.}
#'   \item{\code{alpha}}{One-sided familywise significance level.}
#'   \item{\code{alpha_adj}}{Bonferroni-adjusted significance level
#'     \eqn{\alpha / K}.}
#'   \item{\code{formula_result}}{(\code{approach = "formula"} only.)
#'     A named list with \code{Power}, \code{RCP_M1}, \code{RCP_M2}.
#'     \code{NULL} for \code{approach = "simulation"}.}
#'   \item{\code{sim_results}}{(\code{approach = "simulation"} only.)
#'     A data frame with four rows (one per procedure: \code{"bonferroni"},
#'     \code{"holm"}, \code{"hochberg"}, \code{"hommel"}) and columns
#'     \code{mcp_method}, \code{Power}, \code{RCP_M1}, \code{RCP_M2},
#'     all computed from the same random draws. \code{NULL} for
#'     \code{approach = "formula"}.}
#' }
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats p.adjust pnorm pt qnorm qt rnorm rWishart cov2cor
#' @importFrom utils combn
#'
#' @examples
#' # Example 1: Closed-form solution with K = 3 endpoints (Bonferroni)
#' result1 <- rcpmcp_single(
#'   delta    = c(0.2, 0.2, 0.2),
#'   Sigma    = diag(3),
#'   N        = 200,
#'   fs       = c(0.1, 0.45, 0.45),
#'   K        = 3,
#'   gamma_M1 = 0.5,
#'   gamma_M2 = 0,
#'   alpha    = 0.025
#' )
#' print(result1)
#'
#' \donttest{
#' # Example 2: Simulation (all four MCP methods from a single random draw)
#' result2 <- rcpmcp_single(
#'   delta    = c(0.2, 0.2, 0.2),
#'   Sigma    = diag(3),
#'   N        = 200,
#'   fs       = c(0.1, 0.45, 0.45),
#'   K        = 3,
#'   gamma_M1 = 0.5,
#'   gamma_M2 = 0,
#'   alpha    = 0.025,
#'   approach = "simulation",
#'   nsim     = 10000,
#'   seed     = 1
#' )
#' print(result2)
#'
#' # Example 3: Print selected MCP methods only
#' print(result2, mcp_method = c("bonferroni", "holm"))
#' }
#'
#' \donttest{
#' # Example 4: Region-specific treatment effects (e.g., null effect in region 1)
#' delta_mat <- matrix(c(0.0, 0.0, 0.0,
#'                       0.2, 0.2, 0.2,
#'                       0.2, 0.2, 0.2), nrow = 3, ncol = 3, byrow = FALSE)
#' result4 <- rcpmcp_single(
#'   delta    = delta_mat,
#'   Sigma    = diag(3),
#'   N        = 200,
#'   fs       = c(0.1, 0.45, 0.45),
#'   K        = 3,
#'   gamma_M1 = 0.5,
#'   gamma_M2 = 0,
#'   alpha    = 0.025,
#'   approach = "simulation",
#'   nsim     = 10000,
#'   seed     = 1
#' )
#' print(result4)
#' }
#'
#' \donttest{
#' # Example 5: Unknown variance (Wishart-sampled sample covariance)
#' result5 <- rcpmcp_single(
#'   delta          = c(0.2, 0.2, 0.2),
#'   Sigma          = diag(3),
#'   N              = 200,
#'   fs             = c(0.1, 0.45, 0.45),
#'   K              = 3,
#'   gamma_M1       = 0.5,
#'   gamma_M2       = 0,
#'   alpha          = 0.025,
#'   approach       = "simulation",
#'   nsim           = 10000,
#'   seed           = 1,
#'   variance_known = FALSE
#' )
#' print(result5)
#' }
#'
#' @export
rcpmcp_single <- function(delta,
                          Sigma          = NULL,
                          N,
                          fs             = c(0.1, 0.45, 0.45),
                          K              = 1,
                          gamma_M1       = 0.5,
                          gamma_M2       = 0,
                          alpha          = 0.025,
                          approach       = "formula",
                          nsim           = 1e4,
                          seed           = 1,
                          variance_known = TRUE) {

  # ========== Input Validation ==========
  # delta may be either a numeric vector of length K (region-homogeneous
  # treatment effects) or a numeric K-by-S matrix (region-specific
  # treatment effects, only allowed under approach = "simulation").
  if (is.matrix(delta)) {
    if (!is.numeric(delta) || any(delta < 0)) {
      stop("delta must contain only non-negative numeric values")
    }
    if (nrow(delta) != K) {
      stop("nrow(delta) must equal K when delta is supplied as a matrix")
    }
    if (ncol(delta) != length(fs)) {
      stop("ncol(delta) must equal length(fs) when delta is a matrix")
    }
    delta_is_matrix <- TRUE
  } else {
    if (!is.numeric(delta) || any(delta < 0)) {
      stop("delta must be a numeric vector of non-negative values")
    }
    if (length(delta) != K) {
      stop("length of delta must equal K")
    }
    delta_is_matrix <- FALSE
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
  if (!is.numeric(N) || length(N) != 1 || N <= 0 || N != as.integer(N)) {
    stop("N must be a single positive integer")
  }
  if (!is.numeric(fs) || any(fs <= 0) || abs(sum(fs) - 1) > 1e-8) {
    stop("fs must be a numeric vector of positive values summing to 1")
  }
  if (!is.numeric(K) || length(K) != 1 || K < 1 || K != as.integer(K)) {
    stop("K must be a single positive integer")
  }
  if (!is.numeric(gamma_M1) || length(gamma_M1) != 1 ||
      gamma_M1 < 0 || gamma_M1 > 1) {
    stop("gamma_M1 must be a single numeric value in [0, 1]")
  }
  if (!is.numeric(gamma_M2) || length(gamma_M2) != 1 || gamma_M2 < 0) {
    stop("gamma_M2 must be a single non-negative numeric value")
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
  if (!is.logical(variance_known) || length(variance_known) != 1 ||
      is.na(variance_known)) {
    stop("variance_known must be a single logical value (TRUE or FALSE)")
  }

  # Formula branch only supports region-homogeneous delta and known variance.
  if (approach == "formula") {
    if (delta_is_matrix) {
      stop('region-specific delta (matrix input) is only supported under approach = "simulation"')
    }
    if (!variance_known) {
      stop('variance_known = FALSE is only supported under approach = "simulation"')
    }
  }

  # ========== Common Setup ==========
  Sigma_use <- if (is.null(Sigma)) diag(K) else Sigma

  # ========== Formula Approach (Bonferroni closed-form) ==========
  if (approach == "formula") {

    # Pre-compute all gamma-independent quantities once
    pre <- rcpmcp_single_precompute(
      delta     = delta,
      Sigma_use = Sigma_use,
      N         = N,
      fs        = fs,
      K         = K,
      alpha     = alpha
    )

    # Compute Power (independent of gamma_M1 / gamma_M2)
    power_out <- .rcpmcp_formula_power(pre)

    # Compute RCP for M1 and M2
    rcp_res <- .rcpmcp_formula_rcp_only(
      pre       = pre,
      delta     = delta,
      Sigma_use = Sigma_use,
      gamma_M1  = gamma_M1,
      gamma_M2  = gamma_M2,
      power_out = power_out
    )

    formula_result <- list(
      Power  = power_out,
      RCP_M1 = rcp_res$RCP_M1,
      RCP_M2 = rcp_res$RCP_M2
    )
    sim_results <- NULL
    alpha_adj   <- pre$alpha_adj

    # ========== Simulation Approach ==========
  } else {

    # Derive Ns with floor() + remainder correction (sum(Ns) == N exactly)
    S     <- length(fs)
    Ns    <- floor(fs * N)
    Ns[S] <- N - sum(Ns[-S])

    sd_j       <- sqrt(diag(Sigma_use))
    se_overall <- sd_j / sqrt(N)
    alpha_adj  <- alpha / K

    # Build the K-by-S matrix of region-specific true treatment effects.
    # When delta is supplied as a vector, broadcast it across all S regions.
    if (delta_is_matrix) {
      delta_mat <- delta
    } else {
      delta_mat <- matrix(delta, nrow = K, ncol = S, byrow = FALSE)
    }

    set.seed(seed)

    # ------------------------------------------------------------------
    # Generate a single set of random draws shared across all MCP methods.
    # hat_delta_s_arr: nsim-by-K-by-S array of regional treatment effect
    # estimates. For region s, the K-dimensional vector of regional means
    # follows a multivariate normal distribution with mean delta_mat[, s]
    # and variance-covariance matrix Sigma_use / Ns[s]. Regions are
    # independent.
    # ------------------------------------------------------------------
    hat_delta_s_arr <- array(NA_real_, dim = c(nsim, K, S))
    for (s in seq_len(S)) {
      Sigma_s <- Sigma_use / Ns[s]
      L_s     <- t(chol(Sigma_s))
      Z_raw   <- matrix(stats::rnorm(nsim * K), nrow = K, ncol = nsim)
      hat_s   <- L_s %*% Z_raw + delta_mat[, s]
      hat_delta_s_arr[, , s] <- t(hat_s)
    }

    # Overall estimate: hat_delta_j = sum_s fs[s] * hat_delta_{s,j}
    hat_delta_overall <- matrix(0, nrow = nsim, ncol = K)
    for (s in seq_len(S)) {
      hat_delta_overall <- hat_delta_overall +
        fs[s] * hat_delta_s_arr[, , s]
    }

    # ------------------------------------------------------------------
    # Compute one-sided p-values for the overall test of each endpoint.
    #
    # When variance_known = TRUE, Wald-type Z-statistics are compared to
    # the standard normal distribution.
    #
    # When variance_known = FALSE, the sample variance-covariance matrix
    # is sampled from a Wishart distribution with nu = N - 2 degrees of
    # freedom and scale matrix Sigma_use / (N - 2). This avoids the need
    # to generate individual-level data. The diagonal elements supply
    # endpoint-specific sample variances used in the t-statistic
    # denominator, and p-values are obtained from the t distribution with
    # nu degrees of freedom.
    # ------------------------------------------------------------------
    p_mat <- matrix(NA_real_, nrow = nsim, ncol = K)
    if (variance_known) {
      for (j in seq_len(K)) {
        Z_j        <- hat_delta_overall[, j] / se_overall[j]
        p_mat[, j] <- 1 - stats::pnorm(Z_j)
      }
    } else {
      nu          <- N - 2L
      if (nu < 1L) {
        stop("variance_known = FALSE requires N >= 3 (degrees of freedom >= 1)")
      }
      # rWishart returns a K-by-K-by-nsim array of Wishart(nu, Sigma_use/nu)
      # draws. The diagonal element [j, j, i] is the sample variance for
      # endpoint j in simulation i.
      W_arr  <- stats::rWishart(nsim, df = nu, Sigma = Sigma_use / nu)
      hat_sd <- matrix(NA_real_, nrow = nsim, ncol = K)
      for (j in seq_len(K)) {
        hat_sd[, j] <- sqrt(W_arr[j, j, ])
      }
      for (j in seq_len(K)) {
        T_j        <- hat_delta_overall[, j] / (hat_sd[, j] / sqrt(N))
        p_mat[, j] <- 1 - stats::pt(T_j, df = nu)
      }
    }

    # Consistency flags are independent of mcp_method: compute once.
    # Method 1 and Method 2 criteria use only point estimates of treatment
    # effects, so they do not depend on whether the variance is known.
    consist_m1 <- matrix(FALSE, nrow = nsim, ncol = K)
    consist_m2 <- matrix(FALSE, nrow = nsim, ncol = K)
    for (j in seq_len(K)) {
      consist_m1[, j] <- hat_delta_s_arr[, j, 1] >
        gamma_M1 * hat_delta_overall[, j]
      consist_m2[, j] <- rowSums(hat_delta_s_arr[, j, ] > gamma_M2) == S
    }

    # ------------------------------------------------------------------
    # Apply all four MCP methods to the same p_mat.
    # For each method, compute Power, RCP_M1, RCP_M2.
    # ------------------------------------------------------------------
    all_methods <- c("bonferroni", "holm", "hochberg", "hommel")

    sim_results <- data.frame(
      mcp_method = all_methods,
      Power      = NA_real_,
      RCP_M1     = NA_real_,
      RCP_M2     = NA_real_,
      stringsAsFactors = FALSE
    )

    for (mth in all_methods) {
      signif_mth <- matrix(FALSE, nrow = nsim, ncol = K)
      for (i in seq_len(nsim)) {
        p_adj           <- stats::p.adjust(p_mat[i, ], method = mth)
        signif_mth[i, ] <- p_adj < alpha
      }

      any_signif   <- rowSums(signif_mth) >= 1
      any_joint_m1 <- rowSums(consist_m1 & signif_mth) >= 1
      any_joint_m2 <- rowSums(consist_m2 & signif_mth) >= 1

      pw  <- mean(any_signif)
      rm1 <- if (pw > 0) mean(any_joint_m1) / pw else NA_real_
      rm2 <- if (pw > 0) mean(any_joint_m2) / pw else NA_real_

      sim_results[sim_results$mcp_method == mth,
                  c("Power", "RCP_M1", "RCP_M2")] <- list(pw, rm1, rm2)
    }

    formula_result <- NULL
  }

  # ========== Output ==========
  out <- list(
    approach       = approach,
    nsim           = if (approach == "simulation") nsim else NULL,
    delta          = delta,
    Sigma          = Sigma_use,
    N              = N,
    fs             = fs,
    K              = K,
    gamma_M1       = gamma_M1,
    gamma_M2       = gamma_M2,
    alpha          = alpha,
    alpha_adj      = alpha_adj,
    variance_known = variance_known,
    formula_result = formula_result,
    sim_results    = sim_results
  )
  class(out) <- "rcpmcp_single"
  return(out)
}


#' Print Method for rcpmcp_single Objects
#'
#' @description
#' Print a summary of an \code{"rcpmcp_single"} object.
#'
#' For \code{approach = "formula"}, a single column of results is shown
#' (Bonferroni closed-form).
#'
#' For \code{approach = "simulation"}, results for all four multiplicity
#' adjustment procedures are stored in \code{$sim_results} and displayed
#' in a side-by-side table. The \code{mcp_method} argument selects which
#' procedure(s) to include and in what order.
#'
#' @param x An object of class \code{"rcpmcp_single"}.
#' @param digits Integer. Number of decimal places for probability values.
#'   Default is \code{4}.
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
#' @rdname rcpmcp_single
#' @export
print.rcpmcp_single <- function(x,
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

  if (is.matrix(x$delta)) {
    cat("   Treatment effect  : delta    = (K-by-S matrix)\n")
    rownames_show <- paste0("                                E", seq_len(nrow(x$delta)))
    colnames_show <- paste0("R", seq_len(ncol(x$delta)))
    delta_fmt <- formatC(x$delta, format = "f", digits = digits)
    cat(strrep(" ", 31), paste(formatC(colnames_show, width = max(nchar(delta_fmt)),
                                       flag = " "), collapse = "  "), "\n", sep = "")
    for (k in seq_len(nrow(x$delta))) {
      cat(rownames_show[k], "  ",
          paste(formatC(delta_fmt[k, ], width = max(nchar(delta_fmt)),
                        flag = " "), collapse = "  "), "\n", sep = "")
    }
  } else {
    cat(sprintf("   Treatment effect  : delta    = (%s)\n",
                paste(formatC(x$delta, format = "f", digits = digits),
                      collapse = ", ")))
  }
  cat(sprintf("   Std. deviation    : sd       = (%s)\n",
              paste(formatC(sqrt(diag(x$Sigma)), format = "f", digits = digits),
                    collapse = ", ")))
  cat(sprintf("   Total sample size : N        = %d\n",   x$N))
  cat(sprintf("   Regional props.   : fs       = (%s)\n",
              paste(x$fs, collapse = ", ")))
  cat(sprintf("   No. of endpoints  : K        = %d\n",   x$K))
  cat(sprintf("   Threshold (M1)    : gamma_M1 = %.4f\n", x$gamma_M1))
  cat(sprintf("   Threshold (M2)    : gamma_M2 = %.4f\n", x$gamma_M2))
  cat(sprintf("   Significance lvl  : alpha    = %.4f\n", x$alpha))
  cat(sprintf("   Adjusted level    : alpha/K  = %.4f\n", x$alpha_adj))
  variance_known <- if (is.null(x$variance_known)) TRUE else x$variance_known
  if (x$approach == "simulation") {
    if (variance_known) {
      cat("   Variance          : Known (Z-statistic, normal critical value)\n")
    } else {
      cat(sprintf(
        "   Variance          : Unknown (T-statistic, df = %d, Wishart-sampled)\n",
        x$N - 2L
      ))
    }
  }

  corr_use <- cov2cor(x$Sigma)
  is_indep <- all(abs(corr_use - diag(x$K)) < 1e-8)
  if (is_indep) {
    cat("   Endpoint corr.    : Independent (identity matrix)\n")
  } else {
    cat("   Endpoint corr.    : User-specified (see $Sigma)\n")
  }

  cat(strrep("-", 70), "\n")
  cat("Results:\n\n")

  fmt <- function(val) formatC(val, format = "f", digits = digits)

  # ---- Formula approach: single column ----
  if (x$approach == "formula") {
    cat(sprintf("   Power          : %s\n", fmt(x$formula_result$Power)))
    cat(sprintf("   RCP (Method 1) : %s\n", fmt(x$formula_result$RCP_M1)))
    cat(sprintf("   RCP (Method 2) : %s\n", fmt(x$formula_result$RCP_M2)))

    # ---- Simulation approach: side-by-side table ----
  } else {
    valid_methods <- c("bonferroni", "holm", "hochberg", "hommel")
    mcp_method    <- mcp_method[mcp_method %in% valid_methods]
    if (length(mcp_method) == 0L) {
      warning("No valid mcp_method specified; displaying all four methods.")
      mcp_method <- valid_methods
    }

    # Subset sim_results to requested methods, preserving supplied order
    sr <- x$sim_results[match(mcp_method, x$sim_results$mcp_method), ]

    # Column headers (first letter capitalised)
    cap     <- function(s) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
    headers <- vapply(sr$mcp_method, cap, character(1L))

    # Formatted value rows
    pw_row  <- vapply(sr$Power,  fmt, character(1L))
    rm1_row <- vapply(sr$RCP_M1, fmt, character(1L))
    rm2_row <- vapply(sr$RCP_M2, fmt, character(1L))

    # Column widths
    row_labels <- c("Power", "RCP (Method 1)", "RCP (Method 2)")
    label_w    <- max(nchar(row_labels))
    col_w      <- max(nchar(headers), nchar(pw_row),
                      nchar(rm1_row), nchar(rm2_row), 10L)

    # Header row
    cat(
      strrep(" ", label_w + 5L),
      paste(formatC(headers, width = col_w, flag = " "), collapse = "  "),
      "\n", sep = ""
    )

    # Result rows
    .print_row <- function(label, vals) {
      cat(sprintf("   %-*s  %s\n",
                  label_w, label,
                  paste(formatC(vals, width = col_w, flag = " "),
                        collapse = "  ")))
    }
    .print_row("Power",          pw_row)
    .print_row("RCP (Method 1)", rm1_row)
    .print_row("RCP (Method 2)", rm2_row)
  }

  cat("\n")
  invisible(x)
}
