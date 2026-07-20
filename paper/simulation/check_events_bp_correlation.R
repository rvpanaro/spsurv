#!/usr/bin/env Rscript
# Check correlation between events/BP-coefficient ratio and inflated MLE SEs.

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
source(file.path(paths$sim_dir, "R", "simulation_io.R"), local = TRUE)
source(file.path(paths$render_dir, "mc-paper-summary.R"), local = TRUE)

rep <- read_results_table(find_latest_output("results", paths$sim_output_dir))
cen <- read_censoring_table(find_latest_output("censoring", paths$sim_output_dir))

rep1 <- merge(
  rep,
  cen,
  by = c("nsize", "gdist", "approach", "model", "rep"),
  all.x = TRUE
)
rep1$n_events <- rep1$nsize * rep1$event_proportion
rep1$bp_degree <- as.integer(ceiling(rep1$nsize^0.5))
rep1$n_bp_coef <- rep1$bp_degree
rep1$events_per_bp <- rep1$n_events / rep1$n_bp_coef
rep1$log10_se <- log10(pmax(rep1$se, .Machine$double.eps))
rep1$huge_se <- rep1$se > 1e10

mle <- rep1[rep1$approach == "mle", , drop = FALSE]
age <- mle[mle$par == "age", , drop = FALSE]

gen_cell_label <- function(gdist, model) {
  fam <- ifelse(gdist == "weibull", "W", "LL")
  cls <- ifelse(model == "ph", "PH", ifelse(model == "po", "PO", "AFT"))
  paste0(fam, cls)
}
mle$generator_cell <- gen_cell_label(mle$gdist, mle$model)
mle$par <- ifelse(mle$par == "sexm", "binary", ifelse(mle$par == "age", "continuous", mle$par))
age <- mle[mle$par == "continuous", , drop = FALSE]

cat("=== Replicate level (MLE, all parameters) ===\n")
cat("N rows:", nrow(mle), "\n")
cat("Cor(events_per_bp, log10(se)):", cor(mle$events_per_bp, mle$log10_se, use = "complete.obs"), "\n")
cat("Cor(n_events, log10(se)):", cor(mle$n_events, mle$log10_se, use = "complete.obs"), "\n")
cat("Cor(n_bp_coef, log10(se)):", cor(mle$n_bp_coef, mle$log10_se, use = "complete.obs"), "\n")
cat("Cor(events_per_bp, huge_se):", cor(mle$events_per_bp, as.integer(mle$huge_se), use = "complete.obs"), "\n")

cat("\nMean events_per_bp by huge_se flag:\n")
print(aggregate(events_per_bp ~ huge_se, mle, function(x) {
  c(mean = mean(x), median = median(x), n = length(x))
}))

cat("\n=== Replicate level (MLE, age only) ===\n")
cat("Cor(events_per_bp, log10(se)):", cor(age$events_per_bp, age$log10_se, use = "complete.obs"), "\n")
cat("Cor(events_per_bp, huge_se):", cor(age$events_per_bp, as.integer(age$huge_se), use = "complete.obs"), "\n")
cat("Huge SE rate by nsize:\n")
print(aggregate(huge_se ~ nsize, age, function(x) mean(x)))

cat("\nHuge SE rate by events_per_bp quartile (age):\n")
age$epb_q <- cut(
  age$events_per_bp,
  breaks = quantile(age$events_per_bp, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE
)
print(aggregate(huge_se ~ epb_q, age, function(x) c(rate = mean(x), n = length(x))))

summ <- summarize_bp_mcsim_replicates(rep)

cell_epb <- aggregate(
  events_per_bp ~ nsize + generator_cell + approach + model + par,
  mle,
  mean,
  na.rm = TRUE
)
names(cell_epb)[names(cell_epb) == "events_per_bp"] <- "mean_events_per_bp"
summ2 <- merge(
  summ,
  cell_epb,
  by = c("nsize", "generator_cell", "approach", "model", "par"),
  all.x = TRUE
)
summ2$nsize <- as.integer(as.character(summ2$nsize))
summ2$bp_degree <- as.integer(ceiling(summ2$nsize^0.5))
summ2$log10_calib <- log10(pmax(summ2$calib, .Machine$double.eps))

cat("\n=== Summary cell level (MLE) ===\n")
mle_sum <- summ2[summ2$approach == "MLE", , drop = FALSE]
cat("Cor(mean_events_per_bp, log10(calib)):", cor(mle_sum$mean_events_per_bp, mle_sum$log10_calib, use = "complete.obs"), "\n")
cat("Cor(mean_events_per_bp, calib>1e10):", cor(mle_sum$mean_events_per_bp, as.integer(mle_sum$calib > 1e10), use = "complete.obs"), "\n")
print(
  mle_sum[order(mle_sum$mean_events_per_bp), c(
    "nsize", "generator_cell", "model", "par", "mean_events_per_bp", "bp_degree", "calib"
  )],
  row.names = FALSE
)

fit <- stats::glm(
  huge_se ~ events_per_bp + factor(nsize) + factor(model) + factor(gdist),
  data = age,
  family = stats::binomial
)
cat("\nLogistic (age MLE): huge_se ~ events_per_bp + nsize + model + gdist\n")
print(summary(fit)$coefficients["events_per_bp", , drop = FALSE])

cat("\nWithin nsize=50 (age MLE): cor(events_per_bp, huge_se) =",
    cor(age$events_per_bp[age$nsize == 50], as.integer(age$huge_se[age$nsize == 50]), use = "complete.obs"), "\n")
cat("Within nsize=100 (age MLE): cor(events_per_bp, huge_se) =",
    cor(age$events_per_bp[age$nsize == 100], as.integer(age$huge_se[age$nsize == 100]), use = "complete.obs"), "\n")

cat("\nSpearman (age MLE): events_per_bp vs log10(se):",
    cor(age$events_per_bp, age$log10_se, method = "spearman", use = "complete.obs"), "\n")

cat("\nHuge SE rate by model (age MLE):\n")
print(aggregate(huge_se ~ model, age, mean))

cat("\nBayes: cor(events_per_bp, log10(se)) (age):",
    cor(rep1$events_per_bp[rep1$approach == "bayes" & rep1$par == "age"],
        rep1$log10_se[rep1$approach == "bayes" & rep1$par == "age"], use = "complete.obs"), "\n")
