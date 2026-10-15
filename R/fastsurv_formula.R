#' Evaluate the formula interface of the analysis functions
#'
#' Internal helper for the formula interface of the estimation and testing
#' functions and of \code{kmcurve_fast()}. It evaluates a formula of the form
#' \code{Surv(time, event) ~ group + strata(s)} in \code{data} and returns the
#' vectors that the vector interface of the calling function takes. The
#' response is read from the call to \code{Surv()} without calling it, so the
#' survival package is not needed.
#'
#' @param formula A formula with the response \code{Surv(time, event)}.
#' @param data A data frame or list in which the variables of the formula are
#'   evaluated, or \code{NULL} for the environment of the formula.
#' @param event The value given as the second argument of the calling
#'   function, or \code{NULL}. A data frame given there is used as
#'   \code{data}.
#' @param group The value given as \code{group}, or \code{NULL}.
#' @param strata The value given as \code{strata}, or \code{NULL}.
#' @param group_mode One of \code{"none"}, \code{"optional"}, and
#'   \code{"required"}: whether the right-hand side must not, may, or must
#'   name a group variable.
#' @param strata_ok Whether \code{strata()} terms are allowed.
#' @param fun The name of the calling function, used in error messages.
#'
#' @return A list with the elements \code{time}, \code{event}, \code{group}
#'   (\code{NULL} without a group variable), \code{strata} (\code{NULL}
#'   without \code{strata()}), and \code{control} (the first level of the
#'   group variable, or \code{NULL}).
#'
#' @keywords internal
#' @noRd
fastsurv_formula <- function(formula, data = NULL, event = NULL, group = NULL,
                             strata = NULL,
                             group_mode = c("none", "optional", "required"),
                             strata_ok = FALSE, fun = "this function") {
  group_mode <- match.arg(group_mode)

  # A data frame given as the second argument is the data, as in
  # coxph(formula, data); any other value there is a mistake
  if (!is.null(event)) {
    if (is.null(data) && is.list(event)) {
      data <- event
    } else {
      stop("with a formula, the event indicator is given in Surv(); the ",
           "second argument must be the data frame or omitted")
    }
  }
  if (!is.null(group)) {
    stop("with a formula, the group variable is given on the right-hand side ",
         "of the formula, not in 'group'")
  }
  if (!is.null(strata)) {
    stop("with a formula, strata are given by strata() in the formula, not ",
         "in 'strata'")
  }
  if (!is.null(data) && !is.list(data)) {
    stop("'data' must be a data frame or a list")
  }
  if (length(formula) != 3L) {
    stop("the formula must have the response Surv(time, event)")
  }

  # Variables not in 'data' are looked up from the environment of the formula
  # (or, for a formula without one, from the caller of the analysis function)
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame(2L)
  eval_var <- function(e) {
    if (is.null(data)) eval(e, env) else eval(e, data, env)
  }

  # Name of the function called by an expression, also for pkg::fun
  head_name <- function(e) {
    h <- e[[1L]]
    if (is.name(h)) return(as.character(h))
    if (is.call(h) && length(h) == 3L && is.name(h[[1L]]) &&
        as.character(h[[1L]]) %in% c("::", ":::")) {
      return(as.character(h[[3L]]))
    }
    ""
  }

  # Response: Surv(time, event), with positional or named arguments
  lhs <- formula[[2L]]
  if (!is.call(lhs) || head_name(lhs) != "Surv") {
    stop("the response of the formula must be Surv(time, event)")
  }
  s_args <- as.list(lhs)[-1L]
  s_nm   <- names(s_args)
  if (is.null(s_nm)) s_nm <- rep("", length(s_args))
  if (length(s_args) != 2L || !all(s_nm %in% c("", "time", "event")) ||
      anyDuplicated(s_nm[s_nm != ""]) > 0L) {
    stop("the response must be Surv(time, event) for right-censored data; ",
         "other forms of Surv() are not supported")
  }
  i_time  <- match("time", s_nm)
  i_event <- match("event", s_nm)
  free    <- which(s_nm == "")
  if (is.na(i_time)) {
    i_time <- free[1L]
    free   <- free[-1L]
  }
  if (is.na(i_event)) i_event <- free[1L]

  # Right-hand side: terms separated by '+', each a group variable, a
  # strata() term, or the intercept 1
  split_plus <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
      c(split_plus(e[[2L]]), split_plus(e[[3L]]))
    } else {
      list(e)
    }
  }
  bad_ops <- c("*", ":", "^", "-", "/", "|", "%in%", "+", "(")
  rhs_msg <- paste0("the right-hand side of the formula must be 1 or a ",
                    "single group variable",
                    if (strata_ok) ", optionally with strata()" else "")
  g_expr  <- NULL
  s_exprs <- list()
  for (e in split_plus(formula[[3L]])) {
    if (is.numeric(e) && length(e) == 1L && e == 1) next
    if (is.name(e) && identical(as.character(e), ".")) {
      stop("'.' is not supported in the formula; name the group variable")
    }
    if (is.call(e) && head_name(e) == "strata") {
      if (!strata_ok) {
        stop("strata() is supported only by survdiff_fast() and coxph_fast()")
      }
      st_args <- as.list(e)[-1L]
      st_nm   <- names(st_args)
      if (length(st_args) == 0L || (!is.null(st_nm) && any(st_nm != ""))) {
        stop("strata() must contain one or more variables and no named ",
             "arguments")
      }
      s_exprs <- c(s_exprs, st_args)
      next
    }
    if ((is.call(e) && head_name(e) %in% bad_ops) || !(is.name(e) ||
                                                      is.call(e))) {
      stop(rhs_msg, "; interactions, other formula operators, and ",
           "constants are not supported")
    }
    if (!is.null(g_expr)) {
      stop(rhs_msg, "; it names more than one variable")
    }
    g_expr <- e
  }
  if (group_mode == "none" && !is.null(g_expr)) {
    stop(fun, " is for a single group; use Surv(time, event) ~ 1")
  }
  if (group_mode == "required" && is.null(g_expr)) {
    stop("the right-hand side of the formula must name the group variable")
  }

  # Evaluate the variables in 'data' (or the environment of the formula)
  time   <- eval_var(s_args[[i_time]])
  event  <- eval_var(s_args[[i_event]])
  grp    <- if (is.null(g_expr)) NULL else eval_var(g_expr)
  s_vals <- lapply(s_exprs, eval_var)

  n    <- length(time)
  lens <- c(length(event), if (!is.null(grp)) length(grp),
            vapply(s_vals, length, integer(1L)))
  if (any(lens != n)) {
    stop("the variables of the formula must have the same length")
  }
  all_vars <- c(list(time, event), if (!is.null(grp)) list(grp), s_vals)
  if (any(vapply(all_vars, anyNA, logical(1L)))) {
    stop("the variables of the formula must not contain missing values; ",
         "remove the incomplete rows first, for example with na.omit()")
  }
  # A status coded 1 (censored) and 2 (event), as in data sets of the
  # survival package, is not an event indicator; Surv(time, status == 2) is
  if (is.numeric(event) && length(event) > 0L && all(event %in% c(1, 2)) &&
      any(event == 2)) {
    stop("the event indicator in Surv() must be 0 (censored) or 1 (event); ",
         "for a status coded 1 and 2, as in the survival package, use ",
         "Surv(time, status == 2)")
  }

  # Several stratification variables are combined through their integer
  # codes, so that labels containing the separator cannot merge two strata
  st <- NULL
  if (length(s_vals) == 1L) {
    st <- s_vals[[1L]]
  } else if (length(s_vals) > 1L) {
    codes <- lapply(s_vals, function(v) match(v, sort(unique(v))))
    st    <- do.call(paste, c(codes, sep = "."))
  }

  # The control group defaults to the first level of the group variable, the
  # reference level of a model formula
  control <- NULL
  if (!is.null(grp)) {
    control <- if (is.factor(grp)) {
      levels(droplevels(grp))[1L]
    } else {
      sort(unique(grp))[1L]
    }
  }

  list(time = time, event = event, group = grp, strata = st,
       control = control)
}
