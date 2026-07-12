#' @keywords internal
#' @noRd
.spbp_interp_surv <- function(t_old, s_old, t_new) {
  t_old <- as.numeric(t_old)
  s_old <- as.numeric(s_old)
  ok <- is.finite(t_old) & is.finite(s_old)
  t_old <- t_old[ok]
  s_old <- s_old[ok]
  if (!length(t_old)) {
    return(rep(NA_real_, length(t_new)))
  }
  if (length(unique(t_old)) < 2L) {
    return(rep(s_old[1L], length(t_new)))
  }
  stats::approx(t_old, s_old, xout = t_new, rule = 2, ties = mean)$y
}

#' @keywords internal
#' @noRd
.spbp_censored_pred_types <- function() {
  c("survival", "time", "linear_pred")
}

#' @keywords internal
#' @noRd
.spbp_curve_conf_types <- function() {
  c("log", "log-log", "plain")
}

#' @keywords internal
#' @noRd
.spbp_resolve_predict_mode <- function(type) {
  if (is.null(type)) {
    return(list(mode = "curve", conf_type = "log"))
  }
  type <- as.character(type)[1L]
  if (type %in% .spbp_censored_pred_types()) {
    return(list(mode = type, conf_type = "log"))
  }
  if (type %in% .spbp_curve_conf_types()) {
    return(list(mode = "curve", conf_type = type))
  }
  if (identical(type, "curve")) {
    return(list(mode = "curve", conf_type = "log"))
  }
  stop(
    "'type' must be one of ",
    paste(c("curve", .spbp_censored_pred_types(), .spbp_curve_conf_types()), collapse = ", "),
    ".",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.spbp_newdata_matrix <- function(object, newdata) {
  p <- length(object$coefficients)
  if (p == 0L) {
    return(matrix(numeric(0), nrow = nrow(newdata), ncol = 0L))
  }
  stats::model.matrix(
    object = object$formula[-2L],
    xlev = object$xlevels,
    data = newdata
  )[, -1L, drop = FALSE]
}

#' @keywords internal
#' @noRd
.spbp_predict_linear_pred <- function(object, newdata) {
  if (is.null(newdata)) {
    newdata <- .spbp_training_data(object)
  }
  X <- .spbp_newdata_matrix(object, newdata)
  p <- length(object$coefficients)
  if (p == 0L) {
    return(rep(0, nrow(newdata)))
  }
  as.vector(X %*% object$coefficients)
}

#' @keywords internal
#' @noRd
.spbp_survival_at_times <- function(object, newdata, eval_time) {
  eval_time <- sort(unique(as.numeric(eval_time)))
  eval_time <- eval_time[is.finite(eval_time) & eval_time >= 0]
  if (!length(eval_time)) {
    stop("'eval_time' must contain at least one finite non-negative value.", call. = FALSE)
  }
  if (is.null(newdata)) {
    pr <- survfit(
      object,
      times = eval_time,
      tidy = TRUE
    )
    surv <- pr$surv
    if (length(surv) != length(eval_time)) {
      surv <- .spbp_interp_surv(pr$time, pr$surv, eval_time)
    }
    return(surv)
  }

  pr <- survfit(
    object,
    newdata = newdata,
    times = eval_time,
    tidy = TRUE
  )
  if (!"id" %in% names(pr)) {
    stop("Unable to obtain survival predictions for new data.", call. = FALSE)
  }
  ids <- sort(unique(pr$id))
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(eval_time))
  for (j in seq_along(ids)) {
    rows <- pr$id == ids[j]
    t_j <- pr$time[rows]
    s_j <- pr$surv[rows]
    out[j, ] <- .spbp_interp_surv(t_j, s_j, eval_time)
  }
  out
}

#' @keywords internal
#' @noRd
.spbp_predict_median_time <- function(object, newdata) {
  if (is.null(newdata)) {
    times <- .spbp_default_survfit_times(object)
    pr <- survfit(object, times = times, tidy = TRUE)
    idx <- which(pr$surv <= 0.5)[1L]
    if (is.na(idx)) {
      return(NA_real_)
    }
    return(pr$time[idx])
  }

  n <- nrow(newdata)
  out <- rep(NA_real_, n)
  times <- .spbp_default_survfit_times(object)
  pr <- survfit(
    object,
    newdata = newdata,
    times = times,
    tidy = TRUE
  )
  ids <- sort(unique(pr$id))
  for (j in seq_along(ids)) {
    rows <- pr$id == ids[j]
    t_j <- pr$time[rows]
    s_j <- pr$surv[rows]
    idx <- which(s_j <= 0.5)[1L]
    out[j] <- if (is.na(idx)) NA_real_ else t_j[idx]
  }
  out
}

#' @keywords internal
#' @noRd
.spbp_predict_survival_censored <- function(object, newdata, eval_time) {
  eval_time <- sort(unique(as.numeric(eval_time)))
  eval_time <- eval_time[is.finite(eval_time) & eval_time >= 0]
  if (!length(eval_time)) {
    stop("'eval_time' must be specified for type = 'survival'.", call. = FALSE)
  }
  if (is.null(newdata)) {
    newdata <- .spbp_training_data(object)
  }
  surv_mat <- .spbp_survival_at_times(object, newdata, eval_time)
  if (is.vector(surv_mat)) {
    surv_mat <- matrix(surv_mat, nrow = 1L)
  }
  pred_list <- lapply(seq_len(nrow(surv_mat)), function(i) {
    data.frame(
      .eval_time = eval_time,
      .pred_survival = surv_mat[i, ],
      stringsAsFactors = FALSE
    )
  })
  out <- data.frame(.pred = I(pred_list), stringsAsFactors = FALSE)
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}

#' @keywords internal
#' @noRd
.spbp_predict_time_censored <- function(object, newdata) {
  if (is.null(newdata)) {
    newdata <- .spbp_training_data(object)
  }
  out <- data.frame(
    .pred_time = .spbp_predict_median_time(object, newdata),
    stringsAsFactors = FALSE
  )
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}

#' @keywords internal
#' @noRd
.spbp_predict_linear_pred_censored <- function(object, newdata) {
  if (is.null(newdata)) {
    newdata <- .spbp_training_data(object)
  }
  out <- data.frame(
    .pred_linear_pred = .spbp_predict_linear_pred(object, newdata),
    stringsAsFactors = FALSE
  )
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
