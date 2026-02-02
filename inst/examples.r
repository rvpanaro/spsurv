library("KMsurv")

data("larynx", package = "KMsurv")
print(larynx)

library("spsurv")

mod_bp <- spbp(
  formula = Surv(time, delta) ~ age + factor(stage),
  model = "ph", data = larynx, init = 0
)
summary(mod_bp)

library("survminer")
newdata <- data.frame(age = 70, stage = unique(larynx$stage))
fit_bp_new <- survfit(mod_bp, newdata)
summary(fit_bp_new)
ggsurvplot(fit_bp_new, data = larynx, conf.int = F, surv.median.line = "hv")

mod_bp2 <- bpph(
  formula = Surv(time, delta) ~ age + factor(stage),
  approach = "bayes", data = larynx, init = "0"
)
summary(mod_bp2)

library("bayesplot")
mcmc_pairs(mod_bp2$posterior$beta, off_diag_fun = "hex")

library("survival")
mod <- coxph(Surv(time, delta) ~ age + factor(stage), data = larynx)
res_bp <- residuals(mod_bp)
res_bp2 <- residuals(mod_bp2)
plot(data.frame(res_bp, res_bp2, residuals(mod)), pch = 19)

library("survival")
data("cancer", package = "survival")
veteran2 <- veteran[veteran$prior == 0, ]
veteran2$celltype <- factor(
  x = veteran2$celltype,
  levels = c("large", "adeno", "smallcell", "squamous")
)
print(veteran2)

mod_bp3 <- bppo(
  formula = Surv(time, status) ~ karno + factor(celltype),
  approach = "bayes", data = veteran2
)
mod_bp4 <- bpaft(
  formula = Surv(time, status) ~ karno + factor(celltype),
  approach = "bayes", data = veteran2
)

newdata2 <- data.frame(karno = c(30, 70), celltype = "squamous")
fit_bp3 <- survfit(mod_bp3, newdata = newdata2)
fit_bp4 <- survfit(mod_bp4, newdata = newdata2)

library("survminer")
ggsurvplot(list(BPAFT = fit_bp4, BPPO = fit_bp3), combine = TRUE, conf.int = T)

library("ggplot2")
library("patchwork")

km_bp3 <- survfit(Surv(resid(mod_bp3, "cox-snell"), veteran2$status) ~ 1)
km_bp4 <- survfit(Surv(resid(mod_bp4, "cox-snell"), veteran2$status) ~ 1)

p1 <- ggplot() +
  geom_step(aes(x = km_bp3$time, y = -log(km_bp3$surv))) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x = "Residuals BPPO", y = "Cumulative hazard") +
  geom_abline(col = c("blue"))

p2 <- ggplot() +
  geom_step(aes(x = km_bp4$time, y = -log(km_bp4$surv))) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x = "Residuals BPAFT", y = "Cumulative hazard") +
  geom_abline(col = c("blue"))

p1 + p2

mod_bp3_null <- bppo(formula = Surv(time, status) ~ 1, approach = "mle", data = veteran2)
mod_bp4_null <- bpaft(formula = Surv(time, status) ~ 1, approach = "bayes", data = veteran2)

p3 <- ggplot() +
  geom_point(aes(
    x = veteran2$karno, residuals(mod_bp3_null),
    color = factor(veteran2$status)
  )) +
  ylim(-4, 2) +
  geom_smooth(aes(x = veteran2$karno, residuals(mod_bp3_null)),
              method = "loess"
  ) +
  labs(y = "Martingale residuals BPPO", x = "PS")

p4 <- ggplot() +
  geom_point(aes(
    x = veteran2$karno, residuals(mod_bp4_null),
    color = factor(veteran2$status)
  )) +
  ylim(-4, 2) +
  geom_smooth(aes(x = veteran2$karno, residuals(mod_bp4_null)),
              method = "loess"
  ) +
  ylab("Martingale residuals BPAFT") +
  xlab("PS")

p3 + p4

p5 <- ggplot() +
  geom_point(aes(
    x = veteran2$karno, resid(mod_bp3, type = "deviance"),
    color = factor(veteran2$status)
  )) +
  labs(y = "Deviance residuals BPPO", x = "PS") +
  theme_survminer()

p6 <- ggplot() +
  geom_point(aes(
    x = veteran2$karno, resid(mod_bp4, type = "deviance"),
    color = factor(veteran2$status)
  )) +
  labs(y = "Deviance residuals BPAFT", x = "PS") +
  theme_survminer()

p5 + p6
