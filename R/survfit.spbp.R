survfit.spbp <- function(spbp, tau = max(spbp$y[,1]), conf.int = .95,
                         newdata = NULL){

  if(is.null(newdata)) newdata <- spbp$means
  else   newdata <- model.matrix(spbp, data = newdata)

  if (is.vector(newdata, "numeric")) {
     if (is.null(names(newdata))) {
      stop("Newdata argument must be a data frame")
    }
    newdata <- data.frame(as.list(newdata))
  }

  Call <- match.call() ## call
  n <- length(spbp$y[, 1]) ## dataset dimension= n x cols
  tab <- table(spbp$y[, 1], spbp$y[, 2])
  col_sum <- tab[, 1] + tab[, 2]

  ### time, n.risk, n.event and n.censor
  time <- as.numeric(rownames(tab))
  n.risk <- as.numeric(n  - c(0 , cumsum(col_sum)[-length(col_sum)]))
  n.event <- as.numeric(tab[, 2])
  n.censor <- as.numeric(tab[, 1])
  ##

  ### surv, cumhaz
  aux <- paste('~', paste(attr(spbp$terms, "term.labels"),
                          collapse = '+'))
  X <- as.matrix(newdata, ncol = spbp$q)
  coef <- spbp$coefficients
  beta <- matrix(coef[1:spbp$q], ncol = 1)
  gamma <- matrix(coef[(spbp$q + 1):length(coef)], ncol = 1)

  degree <- length(gamma)
  model <- spbp$model
  var  <- spbp$var

  ## cumhaz.f

  cumhaz <- cumhaz(t = time, coef = coef, degree = degree,
                            tau = tau, model = model, features = X)
  if(nrow(cumhaz) == 1) cumhaz <- as.numeric(cumhaz)

  surv <- exp(-cumhaz)
  ##

  ### std.err, lower, upper
  grad <- grad(t = time, coef = coef, degree = degree,
                       tau = tau, model = model, features = X)
 print(grad)
  std.err <- NULL
  for(i in 1:length(time)){
    std.err[i] <- sqrt(t(grad[i,]) %*% (var/n) %*% grad[i, ])
    print(std.err[i])
  }

  ## significance  level
  alpha <- 1 - conf.int
  upper <- surv + qnorm(1-(alpha/2)) * std.err
  lower <- surv - qnorm(1-(alpha/2)) * std.err
  ##

  output <- list(n = n,
                 time = time,
                 n.risk = n.risk,
                 n.event = n.event,
                 n.censor = n.censor,
                 surv = surv,
                 type = 'right',
                 cumhaz = cumhaz,
                 std.err= std.err,
                 lower = lower,
                 upper = upper,
                 conf.int = conf.int,
                 call = Call)
  class(output) <- c("survfit.spbp", "survfit")
  return(output)
}

