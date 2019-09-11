survfit.spbp <- function(spbp, tau = max(spbp$y[,1]), conf.int = .95,
                         newdata = rep(0, ncol(model.matrix.spbp(spbp)))){
  Call <- match.call()
  tab <- as.data.frame(table(spbp$y[,1]))
  n <- length(spbp$y[,1])
  time <- as.numeric(levels(tab[, 1]))
  n.risk <- n - cumsum( tab[, 2]) + tab[1,2]
  n.event <- tab[, 2]

  alpha <- (1 - conf_level); ## signif level
  q <- spbp$q
  beta <- matrix(spbp$coef[1:q], ncol = 1) ## coef. vector
  gamma <- matrix(spbp$coef[(q + 1):length(spbp$coef)], ncol = 1)
  degree <- length(gamma)

  X <- newdata

  time_scaled <- spbp$y[,1] / tau ## saled time
  grad <- matrix(NA, nrow = n, ncol = length(spbp$coefficients))

  print(newdata)
  output <- list(n = n,
                 time = time,
                 n.risk,
                 n.event,
                 n.censor,
                 surv,
                 type = 'right',
                 cumhaz,
                 std.err,
                 lower = lower,
                 upper = upper,
                 conf.int = conf.int,
                 call = Call
                 )
  return()
}

