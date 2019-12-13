devtools::install()
devtools::load_all(".")

data("veteran") ## imports from survival package


fit <- spbp(Surv(time, status) ~ karno + factor(celltype),
            data = veteran, approach = "bayes", model = "po", chains = 1, iter = 1000)
print(fit)
summary(fit)
