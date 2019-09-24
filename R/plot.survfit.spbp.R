plot.survfit.spbp <-
  function(){

 ####
  if(spbp$model == "ph"){

    b <-  bp(time = time, m = degree, tau = tau)$b
    B <-  bp(time = time, m = degree, tau = tau)$B
    surv <-  exp(-(B %*% gamma) * exp(pred))

    for(j in (1:q)){ ## coef indexes
      grad[, j] = (B %*% gamma) * ( as.matrix(X[, -j]) %*% as.vector(beta[-j]) + X[,j]) *
        exp(pred) * surv
    }
    for(k in (1:(length(spbp$coef)-q))){ ## polynomial indexes
        grad[, (k + q)] = (-B[, -k] %*% gamma[-k]  + B[, k]) *
          exp(pred) * surv
    }
  }
  else if(spbp$model == "po"){

    cat('void')

  }
  else{
    #
    #     tau = max(time / exp(pred))
    #     time_scaled = spbp$y[,1] / tau
    #
    #     b = bp(time_scaled, m =  degree, tau = tau)$b
    #     B = bp(time_scaled, m =  degree, tau = tau)$B
    #     grad = NULL ## gradient
    #
    #     grad[j] = -(gamma %*% b) * (exp(spbp$linar_predictors) ^ 2) / time
    #     grad[k] = t(gamma[-k]) %*% B[-k, -k] + B[k, k]
    cat('void')
  }
  for(i in 1:n){
    # surv_var[i, ] = t(grad[i, ]) %*% spbp$var %*% grad[i, ] ## delta method
    surv_var[i] = t(grad[i, ]) %*% spbp$var %*% grad[i, ]
  }

  lower = surv - qnorm(1-(alpha/2)) * sqrt(surv_var/n)
  upper = surv + qnorm(1-(alpha/2)) * sqrt(surv_var/n)
  print(lower)
  ord = order(time)
  plot(time[ord], surv[ord])
  H <- function(t){
    exp(-(bp(time = t, m = degree)$b %*% gamma) * exp(pred))
  }
  # curve(, add= T)
  lines(time[ord], upper[ord], lty = 2)
  lines(time[ord], lower[ord], lty = 2)
}
