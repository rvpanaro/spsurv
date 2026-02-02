## --------------- Degree error handling ---------------
handler1 <- function() {
  e <- parent.frame()

  # variable names in parent frame

  vnames <- objects(, envir = e)[!objects(, envir = e) %in% "degree"]

  # "sourcing" the parent.frame
  for (n in vnames) assign(n, get(n, e))

  aux <- match(c("formula", "data"),
    names(Call),
    nomatch = 0
  )

  e$aux <- aux
}

## --------------- Frailty handling ---------------
handler2 <- function() {
  e <- parent.frame()
  # variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for (n in vnames) assign(n, get(n, e))

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
  e$rand <- rand
  e$id <- id
  e$frailty_idx <- frailty_idx
}

## --------------- Priors handling ---------------
handler3 <- function() {
  e <- parent.frame()
  vnames <- objects(, envir = e)

  for (n in vnames) assign(n, get(n, e))

  if (length(priors$beta) > 0) {
    betap <- lapply(priors$beta, read_prior)
  } else {
    betap <- list(c("normal", "0", "5"))
  }
  e$priordist_beta <- sapply(betap, `[[`, 1)
  e$location_beta <- sapply(betap, `[[`, 2)
  e$scale_beta <- sapply(betap, `[[`, 3)

  if (length(priors$gamma) > 0) {
    gammap <- lapply(priors$gamma, read_prior)
  } else {
    gammap <- list(c("lognormal", "0", "5"))
  }
  e$priordist_gamma <- sapply(gammap, `[[`, 1)
  e$location_gamma <- sapply(gammap, `[[`, 2)
  e$scale_gamma <- sapply(gammap, `[[`, 3)

  if (length(priors$frailty) > 0) {
    frailtyp <- lapply(priors$frailty, read_prior)
  } else {
    frailtyp <- list(c("gamma", "1", "1"))
  }
  e$priordist_frailty <- sapply(frailtyp, `[[`, 1)
  e$par1_frailty <- sapply(frailtyp, `[[`, 2)
  e$par2_frailty <- sapply(frailtyp, `[[`, 3)
}

## --------------- Model Frame error handling ---------------
handler4 <- function() {
  e <- parent.frame()
  # variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for (n in vnames) assign(n, get(n, e))

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
