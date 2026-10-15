#' @rdname as.data.frame.fastsurv
#' @importFrom stats coef setNames
#' @method coef fastsurv
#' @export
coef.fastsurv <- function(object, ...) {
  tab <- fastsurv_table(object)
  est <- attr(tab, "estimates")
  if (!any(est)) {
    stop("an object of class '", class(object)[1L], "' is a test without an ",
         "estimate; use as.data.frame() for its statistic and p-value")
  }
  setNames(tab$estimate[est], tab$term[est])
}
