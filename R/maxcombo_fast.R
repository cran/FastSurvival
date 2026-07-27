#' Fast Max-Combo Weighted Log-Rank Test for Two-Group Survival Data
#'
#' @description
#' Computes the max-combo test, the maximum over a set of Fleming-Harrington
#' weighted log-rank statistics, for comparing survival between two groups under
#' non-proportional hazards. The C++ backend evaluates every weighted numerator
#' and the full between-scheme covariance matrix in a single scan over the
#' pooled sorted data, and the p-value is obtained from the multivariate normal
#' distribution implied by the correlation of the component statistics. The test
#' is robust to the shape of the hazard difference because the most extreme of
#' several complementary weights is taken, with the multiplicity accounted for
#' through the joint distribution.
#'
#' @details
#' Each component is a Fleming-Harrington G(rho, gamma) weighted log-rank
#' statistic Z_k = U_k / sqrt(V_kk), where U_k = sum w_k (O_1 - E_1) and the
#' weight is w_k = S(t-)^rho_k (1 - S(t-))^gamma_k evaluated from the
#' left-continuous pooled Kaplan-Meier estimate S(t-). The between-scheme
#' covariance is V_ab = sum w_a w_b v, with v the hypergeometric variance
#' increment shared by all schemes, so the diagonal of V reproduces the
#' single-scheme weighted variances and the off-diagonal entries give the
#' correlation matrix R of the component Z-scores. The sign convention matches
#' \code{\link{survdiff_fast}}: a component Z is negative when the treatment
#' group is favored.
#'
#' The max-combo statistic and its p-value depend on \code{side}. When
#' \code{side = 1}, the statistic is the most negative component, min_k Z_k, so
#' that a negative value favors the treatment group in the same way as
#' \code{\link{survdiff_fast}} with \code{side = 1}. The one-sided p-value is
#' 1 - P(G_1 >= m, ..., G_K >= m) for G distributed as multivariate normal with
#' mean zero and correlation R, where m = min_k Z_k. When \code{side = 2}, the
#' statistic is max_k abs(Z_k) and the p-value is
#' 1 - P(-m <= G_1 <= m, ..., -m <= G_K <= m).
#'
#' The joint normal probability is evaluated by dimension. With a single weight
#' the univariate normal is used. With two or three weights and a one-sided
#' test, where the integration region is a half-space, the deterministic
#' \code{\link[mvtnorm]{TVPACK}} algorithm is used. For the two-sided test, whose
#' region is a bounded rectangle, and for four or more weights, the
#' quasi-Monte-Carlo \code{\link[mvtnorm]{GenzBretz}} algorithm is used, whose
#' precision is governed by \code{abseps} and \code{maxpts}. In a simulation
#' study the Monte Carlo error of the estimated rejection rate is driven by the
#' number of simulated trials rather than by the precision of each individual
#' p-value, so \code{abseps} can be loosened to speed up the four-weight case
#' with negligible effect on the operating characteristics.
#'
#' When \code{presorted = TRUE}, the inputs are assumed to be sorted in ascending
#' order of \code{time} and the internal \code{order()} call is skipped, which is
#' useful inside simulation loops where the data are generated in sorted order.
#'
#' @param time A numeric vector of follow-up times for all subjects.
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group A vector of group labels aligned with \code{time}.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group.
#' @param side An integer, either 1 or 2. If \code{side = 1}, the
#'   one-sided max-combo test for treatment benefit is computed. If
#'   \code{side = 2} (default), the two-sided test based on the maximum absolute
#'   component is computed.
#' @param rho A numeric vector of Fleming-Harrington first parameters, one per
#'   component weight. Defaults to \code{c(0, 0, 1, 1)}.
#' @param gamma A numeric vector of Fleming-Harrington second parameters, one
#'   per component weight, aligned with \code{rho}. Defaults to
#'   \code{c(0, 1, 0, 1)}. The default pairs are the standard four-weight
#'   max-combo: G(0,0) for proportional hazards, G(0,1) for late differences,
#'   G(1,0) for early differences, and G(1,1) for middle differences.
#' @param presorted A logical value. If \code{TRUE}, \code{time}, \code{event},
#'   and \code{group} are assumed to be sorted in ascending order of
#'   \code{time}, and the internal \code{order()} call is skipped. If
#'   \code{FALSE} (default), sorting is handled internally.
#' @param abseps A single positive numeric value, the absolute error tolerance
#'   passed to the multivariate normal integration. Defaults to 1e-5. Larger
#'   values speed up the four-or-more-weight case at the cost of p-value
#'   precision.
#' @param maxpts A single positive integer, the maximum number of function
#'   evaluations for the quasi-Monte-Carlo integration used when four or more
#'   weights are supplied. Defaults to 25000.
#'
#' @return An object of class \code{"maxcombo_fast"}, a named numeric vector of
#'   length two with elements \code{statistic} (the max-combo statistic;
#'   min_k Z_k when \code{side = 1}, so a negative value favors treatment, and
#'   max_k abs(Z_k) when \code{side = 2}) and \code{p.value}. The component
#'   Z-scores are stored in the attribute \code{z}, their correlation matrix in
#'   \code{corr}, the Fleming-Harrington parameters in \code{rho} and
#'   \code{gamma}, the requested \code{side}, and the total sample size in
#'   \code{n}. Returns \code{NA} values (still with class
#'   \code{"maxcombo_fast"}) when any component variance is zero or not finite.
#'
#' @examples
#' library(survival)
#'
#' # Standard four-weight max-combo, one-sided
#' fit <- maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
#' fit["statistic"]
#' fit["p.value"]
#'
#' # Two-sided test
#' maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 1, side = 2)
#'
#' # Custom weight set: proportional plus late-difference only
#' maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 1,
#'               rho = c(0, 0), gamma = c(0, 1))
#'
#' # presorted = TRUE: sort once outside, reuse inside a loop
#' ord <- order(ovarian$futime)
#' maxcombo_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
#'               control = 1, presorted = TRUE)
#'
#' \donttest{
#' # Cross-check against simtrial::maxcombo (one-sided)
#' if (requireNamespace("simtrial", quietly = TRUE)) {
#'   df <- data.frame(
#'     stratum   = "All",
#'     treatment = ifelse(ovarian$rx == 2, "experimental", "control"),
#'     tte       = ovarian$futime,
#'     event     = ovarian$fustat
#'   )
#'   simtrial::maxcombo(df, rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1),
#'                      return_corr = TRUE)$p_value
#' }
#' }
#'
#' @references
#' Karrison, T. G. (2016). Versatile tests for comparing survival curves based
#' on weighted log-rank statistics. \emph{The Stata Journal}, \emph{16}(3),
#' 678-690.
#'
#' Lin, R. S., Lin, J., Roychoudhury, S., et al. (2020). Alternative analysis
#' methods for time to event endpoints under nonproportional hazards: a
#' comparative analysis. \emph{Statistics in Biopharmaceutical Research},
#' \emph{12}(2), 187-198.
#'
#' @seealso
#' \code{\link{survdiff_fast}} for the single-scheme weighted log-rank test.
#'
#' @importFrom stats pnorm cov2cor
#' @importFrom mvtnorm pmvnorm GenzBretz TVPACK
#' @export
maxcombo_fast <- function(time, event, group, control, side = 2,
                          rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1),
                          presorted = FALSE, abseps = 1e-5, maxpts = 25000) {
  # Input validation
  n_total <- length(time)
  if (length(event) != n_total || length(group) != n_total) {
    stop("'time', 'event', and 'group' must have the same length")
  }
  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  if (length(rho) != length(gamma)) {
    stop("'rho' and 'gamma' must have the same length")
  }
  if (length(rho) < 1L) {
    stop("at least one Fleming-Harrington weight must be supplied")
  }
  if (sum(event) == 0L) {
    stop("No events observed in the data")
  }

  nw   <- length(rho)
  side <- as.integer(side)

  # Treatment indicator: 1 = treatment, 0 = control
  if (is.factor(group)) group <- as.character(group)
  j <- as.integer(group != control)

  # Sort pooled data by time when not presorted
  if (!presorted) {
    ord   <- order(time)
    time  <- time[ord]
    event <- as.integer(event[ord])
    j     <- j[ord]
  } else {
    event <- as.integer(event)
  }

  # C++ core: single scan -> list(O1, U (length nw), V (nw x nw))
  res  <- combo_logrank_core(time, event, j, as.numeric(rho), as.numeric(gamma))
  Umat <- res$U
  Vmat <- res$V

  # Component labels, e.g. "FH(0,1)"
  lab <- sprintf("FH(%g,%g)", rho, gamma)

  dV <- diag(Vmat)
  if (any(!is.finite(dV)) || any(dV <= 0)) {
    z_na <- rep(NA_real_, nw)
    names(z_na) <- lab
    na_out <- c(statistic = NA_real_, p.value = NA_real_)
    return(structure(na_out, z = z_na, corr = NULL,
                     rho = rho, gamma = gamma, side = side, n = n_total,
                     control = control,
                     class = "maxcombo_fast"))
  }

  z_vec <- Umat / sqrt(dV)
  names(z_vec) <- lab
  corr_mat <- cov2cor(Vmat)
  dimnames(corr_mat) <- list(lab, lab)

  # Max-combo statistic and joint normal tail probability by dimension
  if (side == 1L) {
    # One-sided: treatment benefit makes a component Z negative, so the
    # strongest evidence is the most negative component. To match the sign
    # convention of survdiff_fast (a negative statistic favors treatment,
    # HR < 1), the reported statistic is min(Z). The p-value uses the
    # lower-tail form of simtrial::pvalue_maxcombo, lower = min(Z) with
    # upper = +Inf.
    m_obs <- min(z_vec)
    lower <- rep(min(z_vec), nw)
    upper <- rep(Inf, nw)
  } else {
    # Two-sided: the most extreme component in absolute value, with the
    # rectangle [-m, m] in every coordinate.
    m_obs <- max(abs(z_vec))
    lower <- rep(-m_obs, nw)
    upper <- rep(m_obs, nw)
  }

  tvpack_ok <- all(lower == -Inf) || all(upper == Inf)
  if (nw == 1L) {
    # Univariate: integral of the standard normal over [lower, upper]
    joint <- pnorm(upper) - pnorm(lower)
  } else if (nw <= 3L && tvpack_ok) {
    joint <- pmvnorm(lower = lower, upper = upper, corr = corr_mat,
                     algorithm = TVPACK(abseps = abseps))[1L]
  } else {
    joint <- pmvnorm(lower = lower, upper = upper, corr = corr_mat,
                     algorithm = GenzBretz(maxpts = maxpts, abseps = abseps,
                                           releps = 0))[1L]
  }

  p_val <- 1 - as.numeric(joint)
  out   <- c(statistic = as.numeric(m_obs), p.value = p_val)

  structure(out, z = z_vec, corr = corr_mat,
            rho = rho, gamma = gamma, side = side, n = n_total,
            control = control,
            class = "maxcombo_fast")
}
