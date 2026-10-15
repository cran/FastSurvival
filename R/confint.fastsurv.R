#' @rdname as.data.frame.fastsurv
#' @importFrom stats confint
#' @method confint fastsurv
#' @export
confint.fastsurv <- function(object, parm, level = NULL, ...) {
  tab <- fastsurv_table(object)
  est <- attr(tab, "estimates")
  if (!any(est)) {
    stop("an object of class '", class(object)[1L], "' is a test without an ",
         "estimate; use as.data.frame() for its statistic and p-value")
  }
  cl <- attr(tab, "conf.level")

  # The limits were computed by the function at its conf.level and cannot be
  # recomputed here for another level
  if (!is.null(level)) {
    if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
        abs(level - cl) > sqrt(.Machine$double.eps)) {
      stop("the confidence limits were computed at conf.level = ", format(cl),
           "; call ", class(object)[1L], "() again with conf.level = ",
           format(level))
    }
  }

  tab <- tab[est, , drop = FALSE]
  if (!missing(parm)) {
    if (is.character(parm)) {
      idx <- match(parm, tab$term)
      if (anyNA(idx)) {
        stop("unknown 'parm': ", paste(parm[is.na(idx)], collapse = ", "))
      }
    } else if (is.numeric(parm)) {
      idx <- as.integer(parm)
      if (anyNA(idx) || any(idx != parm) || any(idx < 1L | idx > nrow(tab))) {
        stop("'parm' must be whole numbers between 1 and ", nrow(tab))
      }
    } else {
      stop("'parm' must be a character or numeric vector")
    }
    tab <- tab[idx, , drop = FALSE]
  }

  a   <- (1 - cl) / 2
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE,
                      digits = 3), "%")
  matrix(c(tab$conf.low, tab$conf.high), ncol = 2L,
         dimnames = list(tab$term, pct))
}
