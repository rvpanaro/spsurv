rm( list = ls(all = T))

# install.packages('bpcp')
require(bpcp)
require(spsurv)

data(leuk2)
table(leuk2$treatment)

# Table 1 gives the ordered
# times for two samples of individuals; censored values are denoted with asterisks.

sort(Surv(leuk2[leuk2$treatment == "6-MP", ]$time, leuk2[leuk2$treatment == "6-MP", ]$status))
sort(Surv(leuk2[leuk2$treatment == "placebo", ]$time, leuk2[leuk2$treatment == "placebo", ]$status))

spsurv::spbp(Surv(time, status) ~ treatment, data = leuk2, approach = 'mle')

t_10 <- leuk2$time -10

spsurv::spbp(Surv(time, status) ~ treatment + t_10, data = leuk2, approach = 'mle')
