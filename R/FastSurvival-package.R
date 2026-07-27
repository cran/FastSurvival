#' @description
#' FastSurvival provides fast alternatives to the standard survival analysis
#' functions in the \pkg{survival} package, together with a simulation layer
#' for designing and evaluating time-to-event trials. Every function is
#' designed for repeated evaluation inside large simulation loops, such as
#' adaptive sample-size re-estimation, probability-of-success calculations, and
#' regional consistency evaluation in multi-regional clinical trials. Core
#' computations are implemented in C++ via \pkg{Rcpp}.
#'
#' @details
#' The estimation and testing functions return S3-class objects that are
#' internally named numeric vectors, so a result can be used directly in
#' arithmetic, subsetting, and aggregation after stripping the class with
#' \code{\link[base]{unclass}}. Each class has a \code{print()} method that
#' formats the result similarly to the corresponding \pkg{survival} package
#' output.
#'
#' Estimation and testing:
#' \describe{
#'   \item{\code{\link{survfit_fast}}}{Kaplan-Meier survival probability,
#'     standard error, and confidence interval at a single specified time
#'     point. The C++ backend locates the evaluation cutoff by binary search
#'     and accumulates the Kaplan-Meier product and the Greenwood sum in a
#'     single scan over event positions only.}
#'   \item{\code{\link{survdiff_fast}}}{Log-rank test for two-group survival
#'     data, including Fleming-Harrington, modestly-weighted, Gehan-Breslow,
#'     and Tarone-Ware weights, and stratified and stratified-weighted
#'     variants. The C++ backend uses a two-pointer merge scan over pooled
#'     sorted vectors.}
#'   \item{\code{\link{coxph_fast}}}{Closed-form hazard ratio estimator via
#'     the Pike-Halley Estimator method, with Wald confidence interval. The
#'     estimator anchors at the Pike estimate and applies a single analytic
#'     Halley correction to the Cox partial likelihood score.}
#'   \item{\code{\link{rmst_fast}}}{Restricted mean survival time for a single
#'     group or a two-group comparison, with difference and ratio contrasts.}
#'   \item{\code{\link{wmst_fast}}}{Window mean survival time over an
#'     interval, generalizing \code{\link{rmst_fast}} to a nonzero lower
#'     window limit.}
#'   \item{\code{\link{milestone_fast}}}{Two-group comparison of Kaplan-Meier
#'     survival at a milestone timepoint, with Wald, log-log, and MOVER
#'     inference methods.}
#'   \item{\code{\link{medsurv_fast}}}{Median survival time for a single group
#'     or a two-group difference, with a kernel-hazard variance method and a
#'     local-constant-hazard method.}
#'   \item{\code{\link{maxcombo_fast}}}{Max-combo test over a set of
#'     Fleming-Harrington weighted log-rank statistics, with the joint p-value
#'     obtained from the implied multivariate normal distribution.}
#'   \item{\code{\link{rmw_fast}}}{Robust modestly-weighted log-rank test, the
#'     maximum of the standard log-rank and a modestly-weighted log-rank
#'     statistic.}
#'   \item{\code{\link{wkm_fast}}}{Weighted Kaplan-Meier (Pepe-Fleming) test,
#'     the weighted integrated difference between two Kaplan-Meier curves.}
#'   \item{\code{\link{ahsw_fast}}}{Average hazard with survival weight,
#'     reporting the ratio and difference contrasts.}
#'   \item{\code{\link{ahr_fast}}}{Kalbfleisch-Prentice average hazard ratio
#'     over a restricted interval, with a test on the group-share scale and an
#'     equivalent test and confidence interval on the log scale.}
#' }
#'
#' Simulation and sequential analysis:
#' \describe{
#'   \item{\code{\link{simdata_fast}}}{Individual patient data simulator for
#'     one-group, two-group, and K-arm time-to-event trials. Supports
#'     piecewise uniform accrual, piecewise exponential survival and dropout
#'     times, subgroups, and correlated two-endpoint illness-death
#'     generation. Random number generation uses \pkg{dqrng}.}
#'   \item{\code{\link{analysis_fast}}}{Interim and final analyses of
#'     simulated trials at a set of information times or calendar cutoffs,
#'     returning one row of test statistics per replicate and look.}
#'   \item{\code{\link{pairwise_fast}}}{Pairwise comparisons of each
#'     experimental arm against a shared control in a multi-arm trial.}
#'   \item{\code{\link{simsummary_fast}}}{Aggregation of replicate-level
#'     results into operating characteristics, given efficacy and futility
#'     boundaries computed elsewhere.}
#' }
#'
#' Visualization:
#' \describe{
#'   \item{\code{\link{gen_scenario_fast}}}{Design-stage assembly of survival,
#'     dropout, and accrual scenarios, with \code{plot()} and \code{print()}
#'     methods.}
#'   \item{\code{\link{kmcurve_fast}}}{Analysis-stage Kaplan-Meier curves with
#'     a risk table, with \code{plot()} and \code{print()} methods.}
#' }
#'
#' @references
#' Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
#' estimator. Manuscript under review.
#'
#' Collett, D. (2014). \emph{Modelling Survival Data in Medical Research}
#' (3rd ed.). Chapman and Hall/CRC.
#'
#' @useDynLib FastSurvival, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"
