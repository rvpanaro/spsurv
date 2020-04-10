#' Inverse Transform Sampling To Generate Time-to-event Data From Parametric Models
#'
#' @export
#' @description Random survival times generation for the weibull or
#'  log-logistic distributions with parameters `scale` and `shape`.
#' @param n integer; sample size
#' @param beta vector of regression coefficients
#' @param lambda  event scale parameter
#' @param x matrix of independent variables (columns)
#' @param k event shape parameter
#' @param censor censoring percentage; approx. mean(status)
#' @param dist "weibull" or "log-logistic"
#'
#' @details rsurv returns weibull or log-logistic(log-logistic) randomly
#' generated survival times. According to Collett (2003), the
#' accelerated failure time model encompasses a wide variety of parametric
#' models, including weibull and log-logistic models.
#'
#' @return data.frame of `ncol(x) +2` columns in which the
#'  survival times are the response variable denoted by `y`,
#'   `status` indicates failure (0 = failure) and dependent variables
#'   are appended to the next columns.
#'
#'@examples rows <- 200
#'
#' categorical <- rbinom(rows, size = 3, prob = .5)
#' x <- data.frame(numerical = rnorm(rows),
#'            cat0 = as.numeric(categorical == 0),
#'            cat1 = as.numeric(categorical == 1),
#'            cat2 = as.numeric(categorical == 2),
#'            cat3 = as.numeric(categorical == 3))
#'
#' newdata <- rsurv(n = rows, beta = c(1, -2, .5, .1, 1),
#'   features = x, model = 'ph', dist = 'weibull')
#'
#' @seealso \code{\link[spsurv]{spbp}}

rsurv <- function(n,
  beta = c(2, -1),
  lambda = .5,
  k = 2,
  x = data.frame(rnorm(n), rnorm(n)),
  dist = c("weibull", "log-logistic"),
  censor = .25){

  delta <- sample(c(rep(T, ceiling(n * censor)),
                    rep(F, floor(n * (1-censor)))))

  if(n %% 1 != 0) stop("n must be integer")
  if(dim(x)[1] != n) stop("Lengths differ: dim(x)[1] = ",
                          dim(x)[1]," but n = ", n)
  if(class(x) != "data.frame") stop("x must be data.frame")

  x <-  model.matrix(~., data = x)[, -1]

  eta <- x %*% beta
  dist <- match.arg(dist)

  if(dist == "weibull"){
      t <- rweibull(n, scale = lambda * exp(eta),
                    shape = k)
    }
  else{
      t <- exp(rlogis(n, location = -(log(lambda) + eta)/k,
                      scale =  exp(location)))
  }
  c <- runif(n = n, min = 0, max =  t)
  y <- t # failure time
  y[delta] <- c[delta] # censoring time
  status <- as.numeric(delta) # event indicator

  if(mean(status) == 0){
    stop("No censoring was generated with the chosen arguments")
  }
  else{
    message("The censoring percentage is ",
            mean(status) * 100, "%, generated from ",
            dist, " distribution.")
  }
  db <- data.frame(time = y, status = status, x = x)
  attr(db, "censoring") <- mean(status)
  attr(db, "dist") <- dist
  attr(db, "scale") <- scale
  attr(db, "shape") <- k
  return(db)
}
