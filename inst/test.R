rm(list = ls(all.names = TRUE))
## Weibull simulated data

### Inverse transform method
n <- 100
beta1 <- 2; beta2 = -1
lambdaT <- .002 # baseline hazard
lambdaC <- .004  # hazard of censoring

x1 <- rnorm(n,0)
x2 <- rnorm(n,0)
# true event time
t <- rweibull(n, shape = 1, scale = lambdaT*exp(-beta1*x1-beta2*x2))
c <- rweibull(n, shape = 1, scale = lambdaC)   #censoring time
time <- pmin(t,c)  #observed time is min of censored and true
event <- time == t   # set to 1 if event is observed

dat <- data.frame(time = time, status = as.numeric(event), x1, x2)
head(dat)

# Use `survreg` in order to check whether the code above has generated a weibull dataset correctly.

### survival
library(survival) #standard model
fit_survival <- survreg(Surv(time, status) ~ x1 + x2, data = dat)
summary(fit_survival)

### spsurv - Accelerated Failure Time (AFT) model
# source('guidelines_load.R')


library(spsurv) #bernstein polynomial based regression
fit <- spbp(Surv(time, event) ~ x1 + x2, data = dat,
            model = 'po', approach = 'mle')
summary(fit)
fit$var


# Larynx dataset

```{r, eval = F}
library(spsurv) ##
library(survival)

data("larynx")
str(larynx)

fit <- spbp(Surv(time, delta) ~  factor(stage) + age,
            approach = 'mle', model = 'ph', data = larynx)
p<- survfit(fit)
summary(p)
summary(fit)

fit_survival <- coxph(Surv(time, delta) ~ factor(stage) + age, data = larynx)
fit_survival$coefficients
summary(fit_survival)
summary(p)
new <- data.frame(age = c(77,88), stage = c(1,2))
p <- survfit(fit_survival, newdata = new)
p2 <- survfit(fit, newdata = new)
names(fit)

names(fit_survival)
summary(p)
View(p)

```


