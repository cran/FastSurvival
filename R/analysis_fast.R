#' Fast Sequential Analysis of Simulated Trial Data
#'
#' @description
#' Performs interim or sequential analyses of simulated two-group
#' time-to-event data at one or more analysis times ("looks"). Each look is
#' defined either by a target cumulative number of events
#' (information-based timing) or by a calendar time (calendar-based timing).
#' At every look the data are administratively censored at the corresponding
#' calendar cutoff, and the requested statistics are computed for each
#' simulated trial by reusing \code{\link{survdiff_fast}},
#' \code{\link{coxph_fast}}, \code{\link{rmst_fast}}, \code{\link{survfit_fast}},
#' \code{\link{maxcombo_fast}}, \code{\link{ahsw_fast}},
#' \code{\link{milestone_fast}}, \code{\link{rmw_fast}},
#' \code{\link{ahr_fast}}, \code{\link{medsurv_fast}},
#' \code{\link{wkm_fast}}, and \code{\link{wmst_fast}}. Optionally the same
#' statistics are also reported within each subgroup. The censoring and time
#' sorting are handled by a C++ backend, and the analysis cores are called with
#' \code{presorted = TRUE}, so each look avoids a redundant sort.
#'
#' @details
#' The input \code{data} is the data frame returned by
#' \code{\link{simdata_fast}} for a two-group trial. The columns \code{sim},
#' \code{group}, \code{accrual_time}, \code{tte}, and \code{event} are
#' required.
#'
#' For a look at calendar time \code{cutoff}, each subject with accrual time
#' \code{a} contributes only if enrolled by then (\code{a <= cutoff}). The
#' observed time at the look is \code{min(tte, cutoff - a)} and the observed
#' event indicator is the original \code{event} when the natural event or
#' dropout occurred on or before \code{cutoff} (\code{a + tte <= cutoff}), and
#' zero otherwise (administrative censoring at the look).
#'
#' When \code{event.looks} is supplied, the calendar cutoff for a target of
#' \code{d} events is the calendar time of the \code{d}-th event in that
#' simulated trial, counted over the whole trial population. If a simulation
#' contains fewer than \code{d} events, the target is never reached: the full
#' data are used, \code{reached} is \code{FALSE}, and \code{cutoff} is
#' \code{NA}. When \code{time.looks} is supplied, the cutoff is the specified
#' calendar time and \code{reached} is always \code{TRUE}. In both cases the
#' cutoff is determined once on the whole population and then used for the
#' overall analysis and for every subgroup analysis at that look.
#'
#' Exactly one of \code{event.looks} and \code{time.looks} must be supplied.
#'
#' The statistics are selected with \code{stat}, which may name one or more of
#' \code{"logrank"}, \code{"coxph"}, \code{"rmst"}, \code{"km"},
#' \code{"maxcombo"}, \code{"ahsw"}, \code{"milestone"}, \code{"rmw"},
#' \code{"ahr"}, \code{"medsurv"}, \code{"wkm"}, and \code{"wmst"}.
#'
#' The \code{"logrank"} statistic is configurable. By default it is the
#' ordinary unweighted, unstratified two-group log-rank test and reproduces the
#' behavior of earlier versions of this function exactly. A non-default
#' \code{weight} selects a weighted log-rank test (Fleming-Harrington,
#' modestly-weighted, Gehan-Breslow, or Tarone-Ware) for non-proportional
#' hazards, and a non-\code{NULL} \code{strata} selects the stratified test,
#' summing the per-stratum contributions. The two options combine to give the
#' stratified weighted log-rank test. Whatever configuration is chosen, the
#' result is written to the same \code{logrank.z}, \code{logrank.chisq}, and
#' \code{logrank.p} columns; the unweighted unstratified case is the \code{K = 1}
#' single-stratum degenerate form of the general statistic. The columns
#' \code{rho}, \code{gamma}, and \code{t_star} parametrize \code{weight}, and
#' \code{strata} names one or more subgroup columns of \code{data} used as the
#' stratification variable. The stratification is determined on the whole cut
#' data, independently of the \code{population} marginalization, so a stratified
#' overall analysis is the canonical primary test; in a single-subgroup
#' population the stratum is constant and the stratified test degenerates to the
#' ordinary one within that subset.
#'
#' The \code{"maxcombo"} statistic is the max-combo test of \code{mc.rho} and
#' \code{mc.gamma} Fleming-Harrington weights. Its \code{maxcombo.stat} is the
#' most extreme component (\code{min} of the component Z-scores when
#' \code{side = 1}, so that a negative value favors treatment, and the maximum
#' absolute component when \code{side = 2}), and \code{maxcombo.p} is the joint
#' multivariate-normal p-value, which already follows \code{side}.
#'
#' The \code{"ahsw"} statistic is the average hazard with survival weight of Uno
#' and Horiguchi on the window from 0 to \code{tau}. It reports the per-group
#' average hazards, the ratio (RAH) and difference (DAH) contrasts with their
#' confidence intervals, and two-sided p-values for both contrasts. The AHSW
#' p-values are always two-sided and do not depend on \code{side}, matching
#' \code{\link{ahsw_fast}}.
#'
#' The \code{"milestone"} statistic compares the Kaplan-Meier survival
#' probabilities of the two groups at the milestone timepoint \code{tau}. It
#' reports the per-group survival, the difference (treatment minus control) with
#' its confidence interval, the test statistic, and the p-value. The inference
#' method is selected with \code{ms.method} (\code{"wald"}, \code{"loglog"}, or
#' \code{"mover"}), matching \code{\link{milestone_fast}}; the benefit direction
#' is a positive difference (higher treatment survival), so a positive Z favors
#' treatment.
#'
#' The \code{"rmw"} statistic is the robust modestly-weighted log-rank test of
#' Magirr and Ohrn, the maximum of the standard log-rank component and a single
#' modestly-weighted component with survival-threshold \code{s_star}. Its
#' \code{rmw.stat} is the most extreme standardized component (the minimum when
#' \code{side = 1}, so a negative value favors treatment, and the maximum
#' absolute component when \code{side = 2}), and \code{rmw.p} is the joint
#' bivariate-normal p-value, which already follows \code{side}, matching
#' \code{\link{rmw_fast}}.
#'
#' The \code{"ahr"} statistic is the Kalbfleisch-Prentice average hazard ratio
#' over the window from 0 to \code{tau}. It reports the average hazard ratio
#' (treatment relative to control), the two group shares of the total hazard,
#' the test statistic on the group-share scale, and the p-value. The benefit
#' direction is an average hazard ratio below 1 (a negative Z), as in
#' \code{\link{ahr_fast}}.
#'
#' The \code{"medsurv"} statistic compares the non-parametric Kaplan-Meier
#' median survival times of the two groups. It reports the per-group medians,
#' their difference (treatment minus control) with a plain-scale confidence
#' interval, the Wald statistic, and the p-value. The variance method is
#' selected with \code{medsurv.method} (\code{"km"} or \code{"nph"}), matching
#' \code{\link{medsurv_fast}}; \code{medsurv.bw} optionally overrides the kernel
#' bandwidth used by \code{medsurv.method = "km"}. The benefit direction is a
#' positive difference (a longer median under treatment), so a positive Z favors
#' treatment.
#'
#' The \code{"wkm"} statistic is the weighted Kaplan-Meier (Pepe-Fleming) test.
#' It reports the weighted integrated survival difference (treatment minus
#' control) with its confidence interval, the Wald statistic, and the p-value.
#' The weight is selected with \code{wkm.weight} (\code{"PF"}, \code{"sqrtPF"},
#' or \code{"constant"}), matching \code{\link{wkm_fast}}. The benefit direction
#' is a positive weighted difference (a positive Z favors treatment).
#'
#' The \code{"wmst"} statistic compares the window mean survival time, the area
#' under the Kaplan-Meier curve between \code{wmst.tau1} and \code{wmst.tau2}.
#' It reports the per-group WMST values, their difference (treatment minus
#' control) with its confidence interval, the Wald statistic, and the p-value,
#' matching \code{\link{wmst_fast}}. With \code{wmst.tau1 = 0} the window mean
#' survival time is the restricted mean survival time at \code{wmst.tau2}. The
#' benefit direction is a positive difference (a positive Z favors treatment).
#'
#' When \code{by.subgroup = TRUE}, the output is given in long form with a
#' \code{population} column. Each \code{(sim, look)} produces one row for the
#' whole trial (\code{population = "overall"}) plus one row per subgroup level
#' of each subgroup factor in \code{data}. Subgroup factors are the columns
#' named \code{subgroup} or \code{subgroup1}, \code{subgroup2}, and so on. Each
#' factor is marginalized separately, so a factor with levels 1 and 2 yields
#' the populations \code{subgroup_1} and \code{subgroup_2} (or
#' \code{subgroup1_1}, \code{subgroup1_2}, and so on for numbered factors).
#' The look cutoff is always determined on the whole population, so subgroup
#' rows at a given look share the same \code{cutoff}, \code{reached}, and
#' \code{look.value}; their \code{n.enrolled} and \code{n.event} are the counts
#' within that subgroup. When \code{by.subgroup = FALSE} (default) the output
#' has no \code{population} column and one row per \code{(sim, look)}, matching
#' the whole-population analysis.
#'
#' @param data A data frame from \code{\link{simdata_fast}} for a two-group
#'   trial, containing at least \code{sim}, \code{group}, \code{accrual_time},
#'   \code{tte}, and \code{event}.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group.
#' @param event.looks A numeric vector of target cumulative event counts, one
#'   per look. Mutually exclusive with \code{time.looks}.
#' @param time.looks A numeric vector of calendar times, one per look.
#'   Mutually exclusive with \code{event.looks}.
#' @param stat A character vector naming the statistics to compute. Any subset
#'   of \code{"logrank"}, \code{"coxph"}, \code{"rmst"}, \code{"km"},
#'   \code{"maxcombo"}, \code{"ahsw"}, \code{"milestone"}, \code{"rmw"},
#'   \code{"ahr"}, \code{"medsurv"}, \code{"wkm"}, and \code{"wmst"}. Defaults to
#'   \code{"logrank"}.
#' @param tau A single positive numeric value, the restriction horizon for
#'   \code{"rmst"}, the truncation time for \code{"ahsw"} and \code{"ahr"}, and
#'   the milestone timepoint for \code{"milestone"}. Required only when
#'   \code{"rmst"}, \code{"ahsw"}, \code{"milestone"}, or \code{"ahr"} is
#'   requested.
#' @param t.eval A single positive numeric value, the landmark time for
#'   \code{"km"}. Required only when \code{"km"} is requested.
#' @param conf.level A single numeric value in (0, 1), the confidence level for
#'   \code{"coxph"}, \code{"rmst"}, and \code{"ahsw"}. Defaults to 0.95.
#' @param side An integer, either 1 or 2. When \code{side = 2} (default),
#'   two-sided p-values \code{2 pnorm(-|z|)} are reported for log-rank, Cox, and
#'   RMST, and the two-sided max-combo test is used. When \code{side = 1}, the
#'   one-sided p-value in the direction of treatment benefit is reported for
#'   each of those statistics, and the one-sided max-combo test is used. The
#'   test statistics themselves are always reported with their natural sign, so
#'   the choice of \code{side} affects only the p-value columns. For log-rank
#'   and Cox the benefit direction is a negative Z (the one-sided p-value is the
#'   lower tail \code{pnorm(z)}); for RMST it is a positive Z (the upper tail
#'   \code{pnorm(-z)}). The AHSW p-values are always two-sided. For
#'   group-sequential boundary comparisons (for example with gsDesign or
#'   rpact), align the sign of the reported Z with the boundary convention
#'   before comparing.
#' @param by.subgroup A logical value. When \code{TRUE}, the analysis is also
#'   reported within each subgroup level, and the output gains a
#'   \code{population} column (long form). When \code{FALSE} (default), only
#'   the whole-population analysis is returned and no \code{population} column
#'   is added. Requires at least one subgroup column in \code{data}.
#' @param weight A character string naming the weight scheme for the
#'   \code{"logrank"} statistic. \code{"logrank"} (default) is the ordinary
#'   unweighted test. \code{"fh"} is the Fleming-Harrington G(rho, gamma) test,
#'   \code{"mwlrt"} is the modestly-weighted log-rank test (requires
#'   \code{t_star}), \code{"gehan"} is Gehan-Breslow, and \code{"tarone-ware"}
#'   is Tarone-Ware. See \code{\link{survdiff_fast}}.
#' @param rho A numeric Fleming-Harrington first parameter, used only when
#'   \code{weight = "fh"}. Defaults to 0.
#' @param gamma A numeric Fleming-Harrington second parameter, used only when
#'   \code{weight = "fh"}. Defaults to 0.
#' @param t_star A single non-negative numeric value, the timepoint of the
#'   modestly-weighted log-rank test. Required only when \code{weight =
#'   "mwlrt"}.
#' @param strata An optional character vector naming one or more subgroup
#'   columns of \code{data} to use as the stratification variable for the
#'   \code{"logrank"} statistic. \code{NULL} (default) gives the unstratified
#'   test. When several columns are named their interaction defines the strata.
#'   Stratification applies only to the \code{"logrank"} statistic; the other
#'   statistics ignore it.
#' @param ms.method A character string naming the inference method for the
#'   \code{"milestone"} statistic, one of \code{"wald"} (default),
#'   \code{"loglog"}, or \code{"mover"}. See \code{\link{milestone_fast}}.
#' @param s_star A single numeric value in (0, 1], the survival-probability
#'   threshold of the modestly-weighted component of the \code{"rmw"} statistic.
#'   The weight is capped at \code{1 / s_star}. Defaults to 0.5. See
#'   \code{\link{rmw_fast}}.
#' @param mc.rho A numeric vector of Fleming-Harrington first parameters for the
#'   \code{"maxcombo"} statistic, one per component weight. Defaults to
#'   \code{c(0, 0, 1, 1)}.
#' @param mc.gamma A numeric vector of Fleming-Harrington second parameters for
#'   the \code{"maxcombo"} statistic, aligned with \code{mc.rho}. Defaults to
#'   \code{c(0, 1, 0, 1)}.
#' @param abseps A single positive numeric value, the absolute error tolerance
#'   passed to the multivariate normal integration of the \code{"maxcombo"}
#'   p-value. Defaults to 1e-5.
#' @param maxpts A single positive integer, the maximum number of function
#'   evaluations for the quasi-Monte-Carlo integration used by the
#'   \code{"maxcombo"} p-value when four or more weights are supplied. Defaults
#'   to 25000.
#' @param medsurv.method A character string naming the variance method for the
#'   \code{"medsurv"} statistic, one of \code{"km"} (default) or \code{"nph"}.
#'   See \code{\link{medsurv_fast}}.
#' @param medsurv.bw An optional single positive numeric value, the kernel
#'   bandwidth for the median hazard used by \code{medsurv.method = "km"},
#'   applied to both groups. \code{NULL} (default) uses a Silverman type default
#'   computed per group from the cut data, matching \code{\link{medsurv_fast}}.
#' @param wkm.weight A character string naming the weight for the \code{"wkm"}
#'   statistic, one of \code{"PF"} (default, Pepe-Fleming combined censoring
#'   weight), \code{"sqrtPF"}, or \code{"constant"}. See \code{\link{wkm_fast}}.
#' @param wmst.tau1 A single non-negative numeric value, the lower window limit
#'   for the \code{"wmst"} statistic. Defaults to 0.
#' @param wmst.tau2 A single positive numeric value, the upper window limit for
#'   the \code{"wmst"} statistic, which must exceed \code{wmst.tau1}. \code{NULL}
#'   (default) falls back to \code{tau}. Required (through either argument) only
#'   when \code{"wmst"} is requested.
#'
#' @return A data frame. When \code{by.subgroup = FALSE}, it has
#'   \code{nsim * length(looks)} rows. When \code{by.subgroup = TRUE}, it has
#'   \code{nsim * length(looks) * (1 + total subgroup levels)} rows and an
#'   extra \code{population} column placed after \code{look.value}. The common
#'   columns are \code{sim}, \code{look} (1-based look index), \code{look.type}
#'   (\code{"event"} or \code{"time"}), \code{look.value} (the requested event
#'   count or calendar time), optionally \code{population}, \code{cutoff} (the
#'   calendar time used, \code{NA} when an event target was not reached),
#'   \code{reached}, \code{n.enrolled}, \code{n.event}, \code{n.dropout} (the
#'   number of enrolled subjects whose dropout occurred on or before the cutoff)
#'   and \code{n.pipeline} (\code{n.enrolled - n.event - n.dropout}, the
#'   subjects still in follow-up at the cutoff), followed by the
#'   columns of the requested statistics. A statistic that cannot be computed
#'   for a row (no events, or an empty group) is \code{NA}. The statistic
#'   columns are \code{logrank.z}, \code{logrank.chisq}, and \code{logrank.p}
#'   for \code{"logrank"}; \code{cox.coef}, \code{cox.hr}, \code{cox.se},
#'   \code{cox.z}, \code{cox.p}, \code{cox.lower}, and \code{cox.upper} for
#'   \code{"coxph"}; \code{rmst.ctrl}, \code{rmst.trt}, \code{rmst.diff},
#'   \code{rmst.diff.lower}, \code{rmst.diff.upper}, \code{rmst.z}, and
#'   \code{rmst.p} for \code{"rmst"}; \code{km.surv.ctrl} and
#'   \code{km.surv.trt} for \code{"km"}; \code{maxcombo.stat} and
#'   \code{maxcombo.p} for \code{"maxcombo"}; and \code{ahsw.ah.ctrl},
#'   \code{ahsw.ah.trt}, \code{ahsw.rah}, \code{ahsw.rah.lower},
#'   \code{ahsw.rah.upper}, \code{ahsw.p.rah}, \code{ahsw.dah},
#'   \code{ahsw.dah.lower}, \code{ahsw.dah.upper}, and \code{ahsw.p.dah} for
#'   \code{"ahsw"}; \code{milestone.surv.ctrl}, \code{milestone.surv.trt},
#'   \code{milestone.diff}, \code{milestone.diff.lower},
#'   \code{milestone.diff.upper}, \code{milestone.z}, and \code{milestone.p}
#'   for \code{"milestone"}; \code{rmw.stat} and \code{rmw.p} for \code{"rmw"};
#'   and \code{ahr.ahr}, \code{ahr.theta.ctrl}, \code{ahr.theta.trt},
#'   \code{ahr.z}, and \code{ahr.p} for \code{"ahr"};
#'   \code{medsurv.ctrl}, \code{medsurv.trt}, \code{medsurv.diff},
#'   \code{medsurv.diff.lower}, \code{medsurv.diff.upper}, \code{medsurv.z}, and
#'   \code{medsurv.p} for \code{"medsurv"}; \code{wkm.wdiff}, \code{wkm.lower},
#'   \code{wkm.upper}, \code{wkm.z}, and \code{wkm.p} for \code{"wkm"}; and
#'   \code{wmst.ctrl}, \code{wmst.trt}, \code{wmst.diff},
#'   \code{wmst.diff.lower}, \code{wmst.diff.upper}, \code{wmst.z}, and
#'   \code{wmst.p} for \code{"wmst"}. The Z columns
#'   \code{logrank.z}, \code{cox.z}, and
#'   \code{rmst.z} carry the natural sign of each test, and the p-value columns
#'   follow \code{side} except for the AHSW p-values, which are two-sided.
#'
#' @examples
#' df <- simdata_fast(
#'   nsim       = 50,
#'   n          = c(150, 150),
#'   a.time     = c(0, 12),
#'   a.rate     = 300 / 12,
#'   e.hazard   = list(list(0.10, 0.07), 0.05),
#'   prevalence = c(0.5, 0.5),
#'   seed       = 1
#' )
#'
#' # Whole-population analysis at two event-based looks
#' res1 <- analysis_fast(df, control = 1, event.looks = c(80, 140))
#' head(res1)
#'
#' # Stratified log-rank on the subgroup factor
#' res2 <- analysis_fast(df, control = 1, time.looks = 24,
#'                       stat = "logrank", strata = "subgroup")
#' head(res2)
#'
#' # Fleming-Harrington G(0, 1) weighted log-rank for delayed effects
#' res3 <- analysis_fast(df, control = 1, time.looks = 24,
#'                       stat = "logrank", weight = "fh", rho = 0, gamma = 1)
#' head(res3)
#'
#' # Max-combo and AHSW together
#' res4 <- analysis_fast(df, control = 1, time.looks = 24,
#'                       stat = c("maxcombo", "ahsw"), tau = 18)
#' head(res4)
#'
#' # Milestone survival, robust modestly-weighted, and average hazard ratio
#' res5 <- analysis_fast(df, control = 1, time.looks = 24,
#'                       stat = c("milestone", "rmw", "ahr"), tau = 18,
#'                       ms.method = "loglog", s_star = 0.5, side = 1)
#' head(res5)
#'
#' @seealso \code{\link{simdata_fast}}, \code{\link{survdiff_fast}},
#'   \code{\link{coxph_fast}}, \code{\link{rmst_fast}},
#'   \code{\link{survfit_fast}}, \code{\link{maxcombo_fast}},
#'   \code{\link{ahsw_fast}}, \code{\link{milestone_fast}},
#'   \code{\link{rmw_fast}}, \code{\link{ahr_fast}},
#'   \code{\link{medsurv_fast}}, \code{\link{wkm_fast}},
#'   \code{\link{wmst_fast}}.
#'
#' @importFrom stats pnorm qnorm cov2cor
#' @export
analysis_fast <- function(data, control,
                          event.looks = NULL, time.looks = NULL,
                          stat = "logrank",
                          tau = NULL, t.eval = NULL,
                          conf.level = 0.95,
                          side = 2,
                          by.subgroup = FALSE,
                          weight = c("logrank", "fh", "mwlrt", "gehan",
                                     "tarone-ware"),
                          rho = 0, gamma = 0, t_star = NULL,
                          strata = NULL,
                          ms.method = c("wald", "loglog", "mover"),
                          s_star = 0.5,
                          mc.rho = c(0, 0, 1, 1), mc.gamma = c(0, 1, 0, 1),
                          abseps = 1e-5, maxpts = 25000,
                          medsurv.method = c("km", "nph"), medsurv.bw = NULL,
                          wkm.weight = c("PF", "sqrtPF", "constant"),
                          wmst.tau1 = 0, wmst.tau2 = NULL) {

  weight <- match.arg(weight)
  ms.method <- match.arg(ms.method)
  medsurv.method <- match.arg(medsurv.method)
  wkm.weight <- match.arg(wkm.weight)

  # ---- Input validation --------------------------------------------------
  req_cols <- c("sim", "group", "accrual_time", "tte", "event")
  if (!is.data.frame(data) || !all(req_cols %in% names(data))) {
    stop("'data' must be a data frame with columns: ",
         paste(req_cols, collapse = ", "))
  }
  has_event <- !is.null(event.looks)
  has_time  <- !is.null(time.looks)
  if (has_event == has_time) {
    stop("supply exactly one of 'event.looks' or 'time.looks'")
  }
  looks     <- if (has_event) event.looks else time.looks
  look_type <- if (has_event) 0L else 1L
  if (length(looks) < 1L || any(!is.finite(looks)) || any(looks <= 0)) {
    stop("'looks' must be positive and finite")
  }

  allowed_stat <- c("logrank", "coxph", "rmst", "km", "maxcombo", "ahsw",
                    "milestone", "rmw", "ahr", "medsurv", "wkm", "wmst")
  if (!all(stat %in% allowed_stat)) {
    stop("'stat' must be a subset of: ", paste(allowed_stat, collapse = ", "))
  }
  stat <- unique(stat)

  if (conf.level <= 0 || conf.level >= 1) stop("'conf.level' must be in (0, 1)")
  if (length(side) != 1L || !side %in% c(1, 2)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  side <- as.integer(side)
  if (("rmst" %in% stat || "ahsw" %in% stat || "milestone" %in% stat ||
       "ahr" %in% stat) &&
      (is.null(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0)) {
    stop("'tau' must be a single positive value when 'rmst', 'ahsw', ",
         "'milestone', or 'ahr' is requested")
  }
  if ("km" %in% stat &&
      (is.null(t.eval) || length(t.eval) != 1L || !is.finite(t.eval) ||
       t.eval <= 0)) {
    stop("'t.eval' must be a single positive value when 'km' is requested")
  }
  if (weight == "mwlrt" &&
      (is.null(t_star) || length(t_star) != 1L || !is.finite(t_star) ||
       t_star < 0)) {
    stop("'t_star' must be a single non-negative value when weight = 'mwlrt'")
  }
  if ("maxcombo" %in% stat && length(mc.rho) != length(mc.gamma)) {
    stop("'mc.rho' and 'mc.gamma' must have the same length")
  }
  if ("rmw" %in% stat &&
      (length(s_star) != 1L || !is.finite(s_star) || s_star <= 0 ||
       s_star > 1)) {
    stop("'s_star' must be a single value in (0, 1] when 'rmw' is requested")
  }
  if (length(wmst.tau1) != 1L || !is.finite(wmst.tau1) || wmst.tau1 < 0) {
    stop("'wmst.tau1' must be a single non-negative value")
  }
  wmst.tau2_use <- if (is.null(wmst.tau2)) tau else wmst.tau2
  if ("wmst" %in% stat) {
    if (is.null(wmst.tau2_use) || length(wmst.tau2_use) != 1L ||
        !is.finite(wmst.tau2_use) || wmst.tau2_use <= wmst.tau1) {
      stop("'wmst.tau2' (or 'tau') must be a single value greater than ",
           "'wmst.tau1' when 'wmst' is requested")
    }
  }
  if (!is.null(medsurv.bw) &&
      (length(medsurv.bw) != 1L || !is.finite(medsurv.bw) ||
       medsurv.bw <= 0)) {
    stop("'medsurv.bw' must be a single positive value or NULL")
  }

  # ---- Subgroup columns and population definitions -----------------------
  sub_cols <- grep("^subgroup([0-9]+)?$", names(data), value = TRUE)
  if (length(sub_cols) > 1L) {
    sfx      <- suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", sub_cols)))
    sub_cols <- sub_cols[order(sfx, sub_cols)]
  }
  if (by.subgroup && length(sub_cols) == 0L) {
    stop("'by.subgroup = TRUE' but 'data' has no subgroup columns")
  }

  use_strata <- !is.null(strata)
  if (use_strata) {
    if (!is.character(strata) || length(strata) < 1L ||
        !all(strata %in% names(data))) {
      stop("'strata' must name one or more columns of 'data'")
    }
  }

  # Population list: overall first, then one entry per (subgroup column, level).
  pop_defs <- list(list(col = NA_integer_, level = NA_integer_,
                        label = "overall"))
  sub_levels <- list()
  if (by.subgroup) {
    for (ci in seq_along(sub_cols)) {
      cn   <- sub_cols[ci]
      levs <- sort(unique(data[[cn]]))
      sub_levels[[cn]] <- levs
      for (lv in levs) {
        pop_defs[[length(pop_defs) + 1L]] <-
          list(col = ci - 1L, level = as.integer(lv),
               label = paste0(cn, "_", lv))
      }
    }
  }
  n_pop  <- length(pop_defs)
  pop_col   <- vapply(pop_defs, function(p) {
    if (is.na(p$col)) -1L else as.integer(p$col)
  }, integer(1))
  pop_level <- vapply(pop_defs, function(p) {
    if (is.na(p$level)) NA_integer_ else as.integer(p$level)
  }, integer(1))
  pop_label <- vapply(pop_defs, function(p) p$label, character(1))

  # ---- Order data by simulation; build sim row offsets -------------------
  # simdata_fast() emits rows grouped by simulation, so the sim column is
  # already non-decreasing; in that case skip order() and the per-column gather
  # and pass the columns through unchanged (the kernel only needs rows
  # contiguous by simulation). Unsorted external data falls back to ordering.
  sim_vec <- data$sim
  if (is.unsorted(sim_vec)) {
    ord   <- order(sim_vec)
    reidx <- function(x) x[ord]
  } else {
    reidx <- function(x) x
  }
  sim_s   <- reidx(sim_vec)
  sim_ids <- unique(sim_s)
  nsim    <- length(sim_ids)
  # 0-based offsets: rows of simulation k are [sim_ptr[k], sim_ptr[k + 1]).
  counts  <- tabulate(match(sim_s, sim_ids), nbins = nsim)
  sim_ptr <- as.integer(c(0L, cumsum(counts)))

  accrual <- reidx(as.numeric(data$accrual_time))
  tte     <- reidx(as.numeric(data$tte))
  event   <- reidx(as.integer(data$event))

  grp <- data$group
  if (is.factor(grp)) grp <- as.character(grp)
  j_all <- reidx(as.integer(grp != control))

  # Subgroup matrix (N x n_subcols), integer labels aligned to ordered data.
  if (length(sub_cols) > 0L) {
    sub_mat <- matrix(0L, nrow = nrow(data), ncol = length(sub_cols))
    for (ci in seq_along(sub_cols)) {
      sub_mat[, ci] <- reidx(as.integer(data[[sub_cols[ci]]]))
    }
  } else {
    sub_mat <- matrix(0L, nrow = nrow(data), ncol = 1L)
  }

  # Stratum identifiers (integer codes), aligned to ordered data.
  if (use_strata) {
    st_lab <- if (length(strata) == 1L) {
      as.character(data[[strata]])
    } else {
      do.call(paste, c(lapply(strata, function(cn) data[[cn]]), sep = "."))
    }
    strata_int <- reidx(as.integer(factor(st_lab)))
  } else {
    strata_int <- integer(0)
  }

  # ---- Statistic flags and weight scheme code ----------------------------
  do_logrank  <- "logrank"  %in% stat
  do_coxph    <- "coxph"    %in% stat
  do_rmst     <- "rmst"     %in% stat
  do_km       <- "km"       %in% stat
  do_maxcombo <- "maxcombo" %in% stat
  do_ahsw     <- "ahsw"     %in% stat
  do_milestone <- "milestone" %in% stat
  do_rmw       <- "rmw"       %in% stat
  do_ahr       <- "ahr"       %in% stat
  do_medsurv   <- "medsurv"   %in% stat
  do_wkm       <- "wkm"       %in% stat
  do_wmst      <- "wmst"      %in% stat

  weight_scheme <- switch(weight,
                          logrank       = -1L,
                          fh            = 0L,
                          mwlrt         = 1L,
                          gehan         = 2L,
                          `tarone-ware` = 3L)
  t_star_v <- if (weight == "mwlrt") t_star else 0
  tau_v    <- if (is.null(tau)) 0 else tau
  teval_v  <- if (is.null(t.eval)) 0 else t.eval
  s_star_v <- if (is.null(s_star)) 0.5 else as.numeric(s_star)
  wkm_weight_v <- switch(wkm.weight, PF = 0L, sqrtPF = 1L, constant = 2L)
  medsurv_bw_v <- if (is.null(medsurv.bw)) -1 else as.numeric(medsurv.bw)
  wmst_tau1_v  <- as.numeric(wmst.tau1)
  wmst_tau2_v  <- if (is.null(wmst.tau2_use)) 0 else as.numeric(wmst.tau2_use)

  # ---- Call the fused kernel ---------------------------------------------
  core <- analysis_loop_core(
    sim_ptr, accrual, tte, event, j_all,
    look_type, as.numeric(looks),
    pop_col, pop_level, sub_mat,
    if (use_strata) strata_int else integer(0), use_strata,
    do_logrank, do_coxph, do_rmst, do_km, do_maxcombo, do_ahsw,
    do_milestone, do_rmw, do_ahr,
    do_medsurv, do_wkm, do_wmst,
    weight_scheme, rho, gamma, t_star_v,
    as.numeric(mc.rho), as.numeric(mc.gamma),
    tau_v, teval_v, s_star_v,
    wmst_tau1_v, wmst_tau2_v, wkm_weight_v, medsurv_bw_v
  )

  total <- nsim * length(looks) * n_pop
  z_mult <- qnorm(1 - (1 - conf.level) / 2)

  # p-value helper matching analysis_fast's p_from_z.
  p_from_z <- function(z, eff_dir) {
    out <- rep(NA_real_, length(z))
    ok  <- !is.na(z)
    if (side == 2L) {
      out[ok] <- 2 * pnorm(-abs(z[ok]))
    } else if (eff_dir < 0) {
      out[ok] <- pnorm(z[ok])
    } else {
      out[ok] <- pnorm(-z[ok])
    }
    out
  }

  # ---- Assemble the output frame in analysis_fast column order -----------
  out <- list()
  # Kernel cell order is (sim outer, then look, then population), so the
  # identifier columns are built to match that nesting exactly.
  out$sim <- rep(sim_ids, each = length(looks) * n_pop)
  out$look <- rep(rep(seq_along(looks), each = n_pop), times = nsim)
  out$look.value <- rep(rep(looks, each = n_pop), times = nsim)
  if (by.subgroup) {
    out$population <- rep(pop_label, times = nsim * length(looks))
  }
  out$cutoff     <- core$cutoff
  out$reached    <- as.logical(core$reached)
  out$n.enrolled <- core$n_enrolled
  out$n.event    <- core$n_event
  out$n.dropout  <- core$n_dropout
  out$n.pipeline <- out$n.enrolled - out$n.event - out$n.dropout

  if (do_logrank) {
    num <- core$logrank[, 1]
    var <- core$logrank[, 2]
    z   <- ifelse(is.finite(var) & var > 0, num / sqrt(var), NA_real_)
    out$logrank.z     <- z
    out$logrank.chisq <- z * z
    out$logrank.p     <- p_from_z(z, eff_dir = -1)
  }
  if (do_coxph) {
    theta0 <- core$coxph[, 1]; U0 <- core$coxph[, 2]
    I0     <- core$coxph[, 3]; J0 <- core$coxph[, 4]
    delta  <- U0 / I0 - (J0 * U0 * U0) / (2 * I0 * I0 * I0)
    theta  <- theta0 * exp(delta)
    se     <- 1 / sqrt(I0)
    coef   <- log(theta)
    z_v    <- ifelse(is.finite(se) & se > 0, coef / se, NA_real_)
    out$cox.coef  <- coef
    out$cox.hr    <- theta
    out$cox.se    <- se
    out$cox.z     <- z_v
    out$cox.p     <- p_from_z(z_v, eff_dir = -1)
    out$cox.lower <- exp(coef - z_mult * se)
    out$cox.upper <- exp(coef + z_mult * se)
  }
  if (do_rmst) {
    r0 <- core$rmst[, 1]; v0 <- core$rmst[, 2]
    r1 <- core$rmst[, 3]; v1 <- core$rmst[, 4]
    est_diff <- r1 - r0
    se_diff  <- sqrt(v1 + v0)
    z_diff   <- ifelse(is.finite(se_diff) & se_diff > 0,
                       est_diff / se_diff, NA_real_)
    out$rmst.ctrl       <- r0
    out$rmst.trt        <- r1
    out$rmst.diff       <- est_diff
    out$rmst.diff.lower <- est_diff - z_mult * se_diff
    out$rmst.diff.upper <- est_diff + z_mult * se_diff
    out$rmst.z          <- z_diff
    out$rmst.p          <- p_from_z(z_diff, eff_dir = 1)
  }
  if (do_km) {
    out$km.surv.ctrl <- core$km[, 1]
    out$km.surv.trt  <- core$km[, 2]
  }
  if (do_maxcombo) {
    nw  <- length(mc.rho)
    lab <- sprintf("FH(%g,%g)", mc.rho, mc.gamma)
    stat_vec <- rep(NA_real_, total)
    p_vec    <- rep(NA_real_, total)
    Umat_all <- core$mc_U
    Vmat_all <- core$mc_V
    for (r in seq_len(total)) {
      U <- Umat_all[r, ]
      if (anyNA(U)) next
      V  <- matrix(Vmat_all[r, ], nrow = nw, ncol = nw)
      dV <- diag(V)
      if (any(!is.finite(dV)) || any(dV <= 0)) next
      zc   <- U / sqrt(dV)
      corr <- cov2cor(V)
      if (side == 1L) {
        m_obs <- min(zc)
        lower <- rep(m_obs, nw); upper <- rep(Inf, nw)
      } else {
        m_obs <- max(abs(zc))
        lower <- rep(-m_obs, nw); upper <- rep(m_obs, nw)
      }
      if (nw == 1L) {
        joint <- pnorm(upper) - pnorm(lower)
      } else if (nw <= 3L) {
        joint <- mvtnorm::pmvnorm(lower = lower, upper = upper, corr = corr,
                                  algorithm = mvtnorm::TVPACK(abseps = abseps))[1L]
      } else {
        joint <- mvtnorm::pmvnorm(lower = lower, upper = upper, corr = corr,
                                  algorithm = mvtnorm::GenzBretz(
                                    maxpts = maxpts, abseps = abseps,
                                    releps = 0))[1L]
      }
      stat_vec[r] <- as.numeric(m_obs)
      p_vec[r]    <- 1 - as.numeric(joint)
    }
    out$maxcombo.stat <- stat_vec
    out$maxcombo.p    <- p_vec
  }
  if (do_ahsw) {
    a0 <- core$ahsw[, 1]; vQ0 <- core$ahsw[, 2]; vU0 <- core$ahsw[, 3]; n0 <- core$ahsw[, 4]
    a1 <- core$ahsw[, 5]; vQ1 <- core$ahsw[, 6]; vU1 <- core$ahsw[, 7]; n1 <- core$ahsw[, 8]
    ok <- is.finite(a0) & is.finite(a1) & a0 > 0 & a1 > 0 &
      is.finite(vQ0) & is.finite(vQ1) & is.finite(vU0) & is.finite(vU1)
    log_rah <- ifelse(ok, log(a1 / a0), NA_real_)
    se_rah  <- ifelse(ok, sqrt(vQ1 / n1 + vQ0 / n0), NA_real_)
    dah     <- ifelse(ok, a1 - a0, NA_real_)
    se_dah  <- ifelse(ok, sqrt(vU1 / n1 + vU0 / n0), NA_real_)
    out$ahsw.ah.ctrl   <- ifelse(ok, a0, NA_real_)
    out$ahsw.ah.trt    <- ifelse(ok, a1, NA_real_)
    out$ahsw.rah       <- ifelse(ok, exp(log_rah), NA_real_)
    out$ahsw.rah.lower <- ifelse(ok, exp(log_rah - z_mult * se_rah), NA_real_)
    out$ahsw.rah.upper <- ifelse(ok, exp(log_rah + z_mult * se_rah), NA_real_)
    out$ahsw.p.rah     <- ifelse(ok & se_rah > 0,
                                 2 * pnorm(-abs(log_rah) / se_rah), NA_real_)
    out$ahsw.dah       <- dah
    out$ahsw.dah.lower <- ifelse(ok, dah - z_mult * se_dah, NA_real_)
    out$ahsw.dah.upper <- ifelse(ok, dah + z_mult * se_dah, NA_real_)
    out$ahsw.p.dah     <- ifelse(ok & se_dah > 0,
                                 2 * pnorm(-abs(dah) / se_dah), NA_real_)
  }

  if (do_milestone) {
    s0 <- core$milestone[, 1]; v0 <- core$milestone[, 2]
    s1 <- core$milestone[, 3]; v1 <- core$milestone[, 4]
    se0 <- sqrt(v0); se1 <- sqrt(v1)
    diff_est <- s1 - s0

    # One-sample confidence limits per group under the chosen transform,
    # vectorized over rows (mirrors milestone_fast's one_sample_ci).
    one_ci <- function(surv, se) {
      lo <- rep(NA_real_, length(surv)); up <- lo
      if (ms.method == "wald") {
        lo <- surv - z_mult * se
        up <- surv + z_mult * se
      } else {
        ok_s <- is.finite(surv) & is.finite(se) & surv > 0 & surv < 1
        if (ms.method == "mover") {
          g_se <- se / surv
          lo[ok_s] <- (surv * exp(-z_mult * g_se))[ok_s]
          up[ok_s] <- (surv * exp(z_mult * g_se))[ok_s]
        } else {
          g_se <- se / (surv * abs(log(surv)))
          lo[ok_s] <- (surv^exp(z_mult * g_se))[ok_s]
          up[ok_s] <- (surv^exp(-z_mult * g_se))[ok_s]
        }
      }
      list(lo = lo, up = up)
    }
    ci0 <- one_ci(s0, se0); ci1 <- one_ci(s1, se1)

    if (ms.method == "wald") {
      se_d  <- sqrt(v0 + v1)
      lower <- diff_est - z_mult * se_d
      upper <- diff_est + z_mult * se_d
      stat  <- ifelse(is.finite(se_d) & se_d > 0, diff_est / se_d, NA_real_)
    } else {
      sigma_l <- (ci1$lo - s1)^2 + (ci0$up - s0)^2
      sigma_u <- (ci1$up - s1)^2 + (ci0$lo - s0)^2
      lower   <- diff_est - sqrt(sigma_l)
      upper   <- diff_est + sqrt(sigma_u)
      if (ms.method == "loglog") {
        ok_m <- is.finite(s0) & is.finite(s1) & s0 > 0 & s0 < 1 &
          s1 > 0 & s1 < 1
        g0  <- log(-log(s0)); g1 <- log(-log(s1))
        vg0 <- v0 / (s0 * log(s0))^2
        vg1 <- v1 / (s1 * log(s1))^2
        stat <- ifelse(ok_m, (g1 - g0) / sqrt(vg0 + vg1), NA_real_)
      } else {
        bind <- ifelse(diff_est >= 0, sigma_l, sigma_u)
        stat <- ifelse(is.finite(bind) & bind > 0,
                       z_mult * diff_est / sqrt(bind), NA_real_)
      }
    }

    out$milestone.surv.ctrl  <- s0
    out$milestone.surv.trt   <- s1
    out$milestone.diff       <- diff_est
    out$milestone.diff.lower <- lower
    out$milestone.diff.upper <- upper
    out$milestone.z          <- stat
    out$milestone.p          <- p_from_z(stat, eff_dir = 1)
  }

  if (do_rmw) {
    U_lr <- core$rmw[, 1]; V_lr <- core$rmw[, 2]
    U_mw <- core$rmw[, 3]; V_mw <- core$rmw[, 4]; Cuv <- core$rmw[, 5]
    stat_vec <- rep(NA_real_, total)
    p_vec    <- rep(NA_real_, total)
    for (r in seq_len(total)) {
      vl <- V_lr[r]; vm <- V_mw[r]
      if (!is.finite(vl) || vl <= 0 || !is.finite(vm) || vm <= 0) next
      z_lr <- U_lr[r] / sqrt(vl)
      z_mw <- U_mw[r] / sqrt(vm)
      rho_r   <- Cuv[r] / sqrt(vl * vm)
      rho_use <- max(min(rho_r, 1 - 1e-10), -1 + 1e-10)
      corr    <- matrix(c(1, rho_use, rho_use, 1), 2L, 2L)
      if (side == 1L) {
        m_obs <- min(z_lr, z_mw)
        lower <- c(m_obs, m_obs); upper <- c(Inf, Inf)
      } else {
        m_obs <- max(abs(z_lr), abs(z_mw))
        lower <- c(-m_obs, -m_obs); upper <- c(m_obs, m_obs)
      }
      tvpack_ok <- all(lower == -Inf) || all(upper == Inf)
      alg <- if (tvpack_ok) mvtnorm::TVPACK() else mvtnorm::Miwa()
      inside <- as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper,
                                            corr = corr, algorithm = alg))
      stat_vec[r] <- as.numeric(m_obs)
      p_vec[r]    <- min(max(1 - inside, 0), 1)
    }
    out$rmw.stat <- stat_vec
    out$rmw.p    <- p_vec
  }

  if (do_ahr) {
    th1   <- core$ahr[, 1]; th2  <- core$ahr[, 2]
    ahr_e <- core$ahr[, 3]; vth2 <- core$ahr[, 5]
    ok <- is.finite(ahr_e) & ahr_e > 0 & is.finite(vth2) & vth2 > 0
    se_th <- ifelse(ok, sqrt(vth2), NA_real_)
    z_ahr <- ifelse(ok, (th2 - 0.5) / se_th, NA_real_)
    out$ahr.ahr        <- ifelse(ok, ahr_e, NA_real_)
    out$ahr.theta.ctrl <- ifelse(ok, th1, NA_real_)
    out$ahr.theta.trt  <- ifelse(ok, th2, NA_real_)
    out$ahr.z          <- z_ahr
    out$ahr.p          <- p_from_z(z_ahr, eff_dir = -1)
  }

  if (do_medsurv) {
    med0 <- core$medsurv[, 1]; gw0 <- core$medsurv[, 2]; kh0 <- core$medsurv[, 3]
    nv0  <- core$medsurv[, 4]; lh0 <- core$medsurv[, 5]
    med1 <- core$medsurv[, 6]; gw1 <- core$medsurv[, 7]; kh1 <- core$medsurv[, 8]
    nv1  <- core$medsurv[, 9]; lh1 <- core$medsurv[, 10]
    if (medsurv.method == "nph") {
      var0 <- ifelse(is.finite(lh0) & lh0 > 0, nv0 / (lh0 * lh0), NA_real_)
      var1 <- ifelse(is.finite(lh1) & lh1 > 0, nv1 / (lh1 * lh1), NA_real_)
    } else {
      var0 <- ifelse(is.finite(kh0) & kh0 > 0, gw0 / (kh0 * kh0), NA_real_)
      var1 <- ifelse(is.finite(kh1) & kh1 > 0, gw1 / (kh1 * kh1), NA_real_)
    }
    diff_v  <- med1 - med0
    se_diff <- sqrt(var0 + var1)
    z_med   <- ifelse(is.finite(se_diff) & se_diff > 0, diff_v / se_diff,
                      NA_real_)
    ci_ok   <- is.finite(diff_v) & is.finite(se_diff)
    out$medsurv.ctrl       <- med0
    out$medsurv.trt        <- med1
    out$medsurv.diff       <- diff_v
    out$medsurv.diff.lower <- ifelse(ci_ok, diff_v - z_mult * se_diff, NA_real_)
    out$medsurv.diff.upper <- ifelse(ci_ok, diff_v + z_mult * se_diff, NA_real_)
    out$medsurv.z          <- z_med
    out$medsurv.p          <- p_from_z(z_med, eff_dir = 1)
  }

  if (do_wkm) {
    num_raw  <- core$wkm[, 1]; variance <- core$wkm[, 2]
    n1_wkm   <- core$wkm[, 3]; n2_wkm   <- core$wkm[, 4]
    n_tot    <- n1_wkm + n2_wkm
    scale_w  <- sqrt(n1_wkm * n2_wkm / n_tot)
    ok       <- is.finite(variance) & variance > 0 &
                is.finite(scale_w) & scale_w > 0
    se_w     <- ifelse(ok, sqrt(variance) / scale_w, NA_real_)
    z_wkm    <- ifelse(ok, num_raw / se_w, NA_real_)
    ci_ok    <- is.finite(num_raw) & is.finite(se_w)
    out$wkm.wdiff <- num_raw
    out$wkm.lower <- ifelse(ci_ok, num_raw - z_mult * se_w, NA_real_)
    out$wkm.upper <- ifelse(ci_ok, num_raw + z_mult * se_w, NA_real_)
    out$wkm.z     <- z_wkm
    out$wkm.p     <- p_from_z(z_wkm, eff_dir = 1)
  }

  if (do_wmst) {
    w0 <- core$wmst[, 1]; vw0 <- core$wmst[, 2]
    w1 <- core$wmst[, 3]; vw1 <- core$wmst[, 4]
    diff_v  <- w1 - w0
    se_diff <- sqrt(vw0 + vw1)
    z_wmst  <- ifelse(is.finite(se_diff) & se_diff > 0, diff_v / se_diff,
                      NA_real_)
    out$wmst.ctrl       <- w0
    out$wmst.trt        <- w1
    out$wmst.diff       <- diff_v
    out$wmst.diff.lower <- diff_v - z_mult * se_diff
    out$wmst.diff.upper <- diff_v + z_mult * se_diff
    out$wmst.z          <- z_wmst
    out$wmst.p          <- p_from_z(z_wmst, eff_dir = 1)
  }

  df <- list2DF(out)
  df
}
