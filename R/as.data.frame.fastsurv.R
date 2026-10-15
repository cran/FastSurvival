#' Data Frame, Estimates, and Confidence Limits of a Result
#'
#' @description
#' The results of the estimation and testing functions have, besides their
#' own class, the class \code{"fastsurv"}. Its \code{as.data.frame()} method
#' returns every estimate and test of a result as a row of a data frame with
#' the same columns for all the functions, so that results can be stacked with
#' \code{rbind()}. \code{coef()} returns the estimates and \code{confint()}
#' their confidence limits.
#'
#' @details
#' The data frame has the columns \code{term}, \code{estimate},
#' \code{std.error}, \code{conf.low}, \code{conf.high}, \code{statistic}, and
#' \code{p.value}, with one row for each quantity in the table below. The
#' terms are the names of the elements of the result. The estimate, its
#' standard error, and its confidence limits are on the same scale. The
#' statistic and the p-value are those of the test that the function reports
#' for that quantity, following its \code{side}, as printed by the
#' \code{print()} method. Entries that the function does not compute are
#' \code{NA}.
#'
#' \tabular{ll}{
#'   \code{\link{survfit_fast}} \tab \code{surv} \cr
#'   \code{\link{survdiff_fast}} \tab \code{logrank}, or the name of the
#'     weight (test only) \cr
#'   \code{\link{coxph_fast}} \tab \code{coef} (the log hazard ratio) \cr
#'   \code{\link{rmst_fast}} \tab \code{rmst}; or \code{rmst.ctrl},
#'     \code{rmst.trt}, \code{diff}, \code{log.ratio} \cr
#'   \code{\link{wmst_fast}} \tab \code{wmst}; or \code{wmst.control},
#'     \code{wmst.treatment}, \code{diff} \cr
#'   \code{\link{medsurv_fast}} \tab \code{median}; or \code{median.control},
#'     \code{median.treatment}, \code{diff} \cr
#'   \code{\link{milestone_fast}} \tab \code{surv.control},
#'     \code{surv.treatment}, \code{diff} \cr
#'   \code{\link{maxcombo_fast}} \tab \code{maxcombo} and the component
#'     statistics, such as \code{FH(0,1)} (test only) \cr
#'   \code{\link{rmw_fast}} \tab \code{rmw}, \code{logrank}, \code{mwlrt}
#'     (test only) \cr
#'   \code{\link{wkm_fast}} \tab \code{wdiff} \cr
#'   \code{\link{ahsw_fast}} \tab \code{ah.ctrl}, \code{ah.trt},
#'     \code{log.rah}, \code{dah} \cr
#'   \code{\link{ahr_fast}} \tab \code{log.ahr}, \code{theta}
#' }
#'
#' Ratios are reported on the log scale, on which their confidence intervals
#' are symmetric and their tests are computed: \code{coef} of
#' \code{coxph_fast()} is the log hazard ratio (as from \code{coef()} of
#' \code{survival::coxph()} with the group variable as a factor), and
#' \code{log.ratio}, \code{log.rah}, and \code{log.ahr} are the logarithms of
#' the ratio of restricted mean survival times, the ratio of average hazards,
#' and the average hazard ratio. Take
#' \code{exp()} of the estimate and of the confidence limits for the ratio.
#' \code{rmst_fast()} and \code{ahsw_fast()} do not store the standard error of
#' a log ratio or of the difference of average hazards; it is recovered from
#' the confidence limits, which are symmetric on these scales, as
#' \code{(conf.high - conf.low) / (2 q)}, where \code{q} is the normal quantile
#' of the confidence level, and the statistic of \code{log.rah} and \code{dah}
#' is the estimate divided by this standard error.
#'
#' The rows of a test only carry the statistic: the value of
#' \code{survdiff_fast()} (the Z-score when \code{side = 1}, the chi-square
#' when \code{side = 2}), the max-combo and rMW statistics with their
#' components, which have no p-value of their own. For \code{milestone_fast()},
#' the statistic of \code{diff} is that of \code{method}, which for
#' \code{"loglog"} is on the complementary log-log scale; a standard error is
#' given only for \code{method = "wald"}, whose interval is symmetric. For
#' \code{ahr_fast()}, \code{theta} is the treatment group's share of the total
#' hazard with its standard error and the primary test of
#' \code{ahr_fast()}, and \code{log.ahr} carries the equivalent test on the log
#' scale.
#'
#' \code{coef()} returns the \code{estimate} column named by \code{term}, and
#' \code{confint()} the confidence limits as a matrix, for the rows that are
#' estimates. Both give an error for a test only. The limits are those computed
#' by the function at its \code{conf.level}; they cannot be recomputed at
#' another level, so \code{level} must be \code{NULL} or equal to that
#' \code{conf.level}.
#'
#' @param x,object A result of an estimation or testing function, of class
#'   \code{"fastsurv"}.
#' @param row.names \code{NULL} or a character vector of row names.
#' @param optional Ignored; the column names are fixed.
#' @param parm The estimates for which limits are returned, as terms or as
#'   positions among the estimates. All estimates by default.
#' @param level \code{NULL} (default) or the confidence level at which the
#'   result was computed.
#' @param ... Further arguments, currently ignored.
#'
#' @return \code{as.data.frame()} returns a data frame with the columns
#'   described in Details. \code{coef()} returns a named numeric vector, and
#'   \code{confint()} a matrix with one row per estimate and two columns
#'   labeled by the percentages of the limits, as \code{stats::confint()}.
#'
#' @seealso \code{\link{formula_interface}}.
#'
#' @examples
#' library(survival)
#'
#' fit <- coxph_fast(Surv(futime, fustat) ~ rx, data = ovarian)
#' as.data.frame(fit)
#' coef(fit)
#' exp(confint(fit))
#'
#' as.data.frame(rmst_fast(Surv(futime, fustat) ~ rx, data = ovarian,
#'                         tau = 600))
#'
#' # Results of a loop (with the vector interface) stacked into a data frame
#' res <- lapply(1:3, function(s) {
#'   set.seed(s)
#'   i <- sample(nrow(ovarian), replace = TRUE)
#'   fit_s <- coxph_fast(ovarian$futime[i], ovarian$fustat[i], ovarian$rx[i],
#'                       control = 1)
#'   cbind(sim = s, as.data.frame(fit_s))
#' })
#' do.call(rbind, res)
#'
#' @method as.data.frame fastsurv
#' @export
as.data.frame.fastsurv <- function(x, row.names = NULL, optional = FALSE,
                                   ...) {
  out <- fastsurv_table(x)
  attr(out, "estimates")  <- NULL
  attr(out, "conf.level") <- NULL
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}
