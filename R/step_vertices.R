#' Vertices of a step function
#'
#' Internal helper that returns the vertices of a right-continuous step
#' function, whose value \code{ys[i]} holds from \code{xs[i]} to
#' \code{xs[i + 1]}, ending at the point \code{(xs[k], ys[k])}. Used to draw a
#' shaded band between two step functions with \code{polygon()}.
#'
#' @param xs A numeric vector of abscissae in ascending order.
#' @param ys A numeric vector of values, of the same length as \code{xs}.
#'
#' @return A list with the vectors \code{x} and \code{y} of the vertices,
#'   each of length \code{2 * length(xs) - 1}.
#'
#' @keywords internal
#' @noRd
step_vertices <- function(xs, ys) {
  k <- length(xs)
  list(x = c(xs[1L], rep(xs[-1L], each = 2L)),
       y = c(rep(ys[-k], each = 2L), ys[k]))
}
