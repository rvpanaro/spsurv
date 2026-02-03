#'
#' @keywords internal
.handler1 <- function(Call) {
  aux <- match(c("formula", "data"), names(Call), nomatch = 0)
  return(aux)
}

#'
#' @keywords internal
.handler2 <- function(temp, formula) {
  id <- NULL
  if (!is.null(attr(temp$formula, "specials")$frailty)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty]
    rand <- 1 ## gamma
  } else if (!is.null(attr(temp$formula, "specials")$frailty.gamma)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty.gamma
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.gamma]
    rand <- 1 ## gamma
  } else if (!is.null(attr(temp$formula, "specials")$frailty.gauss)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty.gauss
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.gauss]
    rand <- 2 ## gauss
  } else if (!is.null(attr(temp$formula, "specials")$frailty.t)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty.t
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.t]
    rand <- 3 ## t-student
  } else {
    rand <- 0
    frailty_idx <- NULL
  }
  list(rand = rand, id = id, frailty_idx = frailty_idx)
}

#'
#' @keywords internal
.handler3 <- function(priors) {
  if (length(priors$beta) > 0) {
    betap <- lapply(priors$beta, read_prior)
  } else {
    betap <- list(c("normal", "0", "5"))
  }
  priordist_beta <- sapply(betap, `[[`, 1)
  location_beta <- sapply(betap, `[[`, 2)
  scale_beta <- sapply(betap, `[[`, 3)

  if (length(priors$gamma) > 0) {
    gammap <- lapply(priors$gamma, read_prior)
  } else {
    gammap <- list(c("lognormal", "0", "5"))
  }
  priordist_gamma <- sapply(gammap, `[[`, 1)
  location_gamma <- sapply(gammap, `[[`, 2)
  scale_gamma <- sapply(gammap, `[[`, 3)

  if (length(priors$frailty) > 0) {
    frailtyp <- lapply(priors$frailty, read_prior)
  } else {
    frailtyp <- list(c("gamma", "1", "1"))
  }
  priordist_frailty <- sapply(frailtyp, `[[`, 1)
  par1_frailty <- sapply(frailtyp, `[[`, 2)
  par2_frailty <- sapply(frailtyp, `[[`, 3)

  list(
    priordist_beta = priordist_beta,
    location_beta = location_beta,
    scale_beta = scale_beta,
    priordist_gamma = priordist_gamma,
    location_gamma = location_gamma,
    scale_gamma = scale_gamma,
    priordist_frailty = priordist_frailty,
    par1_frailty = par1_frailty,
    par2_frailty = par2_frailty
  )
}

#'
#' @keywords internal
.handler4 <- function(mf, Y, type, Terms, formula) {
  if (nrow(mf) == 0) stop("Only missing observations")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  if (type != "right" && type != "counting") {
    stop(paste("spsurv doesn't support \"", type,
      "\" survival data",
      sep = ""
    ))
  }

  if (length(attr(Terms, "variables")) > 2) { # a ~1 formula has length 2
    ytemp <- terms.inner(formula)[1:2]
    xtemp <- terms.inner(formula)[-c(1, 2)]
    if (any(!is.na(match(xtemp, ytemp)))) {
      stop("a variable appears on both the left and right sides of
                  the formula")
    }
  }
}
