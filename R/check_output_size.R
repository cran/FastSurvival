#' Check that a simulated data set fits in an R vector
#'
#' Internal helper that stops with an informative error when the number of
#' rows of a simulated data set, \code{nsim} times the total sample size,
#' exceeds the largest integer, which is the limit of the C++ kernels.
#'
#' @param nsim The number of simulated trials.
#' @param n_grp An integer vector of per-group sample sizes.
#'
#' @return \code{NULL}, invisibly. An error is raised when the limit is
#'   exceeded.
#'
#' @keywords internal
#' @noRd
check_output_size <- function(nsim, n_grp) {
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1) {
    stop("'nsim' must be a positive whole number")
  }
  if (as.numeric(nsim) * sum(as.numeric(n_grp)) > .Machine$integer.max) {
    stop("nsim * sum(n) exceeds ", .Machine$integer.max, " rows; ",
         "reduce 'nsim' and simulate in batches")
  }
  invisible(NULL)
}
