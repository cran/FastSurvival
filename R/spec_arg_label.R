#' Name of a hazard specification in error messages
#'
#' Internal helper for the checks of the illness-death model in
#' \code{\link{simdata_fast}}. The median form of a hazard (for example
#' \code{h01.median}) is converted to the hazard before these checks, so the
#' message names both forms of a hazard argument.
#'
#' @param nm A single character, the internal name of the argument, such as
#'   \code{"h01.hazard"} or \code{"switch.prop"}.
#'
#' @return A single character, the quoted name, followed for a hazard by the
#'   quoted name of its median form in parentheses.
#'
#' @keywords internal
#' @noRd
spec_arg_label <- function(nm) {
  if (grepl("\\.hazard$", nm)) {
    paste0("'", nm, "' (or '", sub("\\.hazard$", ".median", nm), "')")
  } else {
    paste0("'", nm, "'")
  }
}
