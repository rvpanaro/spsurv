# Regenerate paper/spsurv.TeX illustration figures (001--006) and verbatim outputs.
# Figures 007--008 need simulation artifacts in paper/; use inst/render-paper-figures.R
# for the full 001--008 pipeline. PDFs are written to figures/ for paper/spsurv.TeX.
# Usage (from package root):
#   Rscript inst/render-tex-assets.R

pkg_root <- Sys.getenv("SPSURV_ROOT", unset = normalizePath("..", winslash = "/"))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  pkg_root <- normalizePath(".", winslash = "/")
}
setwd(pkg_root)

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(KMsurv)
  library(survival)
  library(ggsurvfit)
  library(ggplot2)
  library(patchwork)
})

fig_dir <- file.path(pkg_root, "figures")
paper_dir <- file.path(pkg_root, "paper")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paper_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(pkg_root, "inst", "tex-verbatim-wrap.R"))

.bp_tidy_survfit_long <- function(sf, curve_labels) {
  td <- ggsurvfit::tidy_survfit(sf)
  k <- length(curve_labels)
  parts <- vector("list", k)
  for (j in seq_len(k)) {
    parts[[j]] <- data.frame(
      time = td$time,
      estimate = td[[paste0("estimate.", j)]],
      conf.low = td[[paste0("conf.low.", j)]],
      conf.high = td[[paste0("conf.high.", j)]],
      series = curve_labels[j],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, parts)
}

tex_lines <- function(x, width = 65L) {
  x <- capture.output(x)
  if (!length(x)) {
    return("")
  }
  wrapped <- tex_verbatim_wrap_lines(paste0("#> ", x), width = width)
  paste(wrapped, collapse = "\n")
}

save_fig <- function(plot, filename, width = 7, height = 5) {
  path <- file.path(fig_dir, filename)
  grDevices::pdf(path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  message("Wrote ", path)
  invisible(path)
}

# --- Larynx MLE -----------------------------------------------------------
data("larynx", package = "KMsurv")
larynx$stage <- factor(larynx$stage)

bpph_fit <- bpph(
  Surv(time, delta) ~ age + stage,
  data = larynx,
  approach = "mle"
)

fragments <- list(
  larynx_setup = paste(
    c(
      'data("larynx", package = "KMsurv")',
      "larynx$stage <- factor(larynx$stage)",
      "",
      "bpph_fit <- bpph(",
      "  Surv(time, delta) ~ age + stage,",
      "  data = larynx,",
      '  approach = "mle"',
      ")"
    ),
    collapse = "\n"
  ),
  larynx_tidy_code = 'print(bpph_fit, what = "tidy")',
  larynx_tidy_out = tex_lines({
    print(bpph_fit, what = "tidy")
  }),
  larynx_glance_code = 'print(bpph_fit, what = "glance")',
  larynx_glance_out = tex_lines({
    print(bpph_fit, what = "glance")
  }),
  larynx_print_code = "print(bpph_fit)",
  larynx_print_out = tex_lines({
    print(bpph_fit)
  })
)

nd <- data.frame(
  age = 65,
  stage = factor(levels(larynx$stage))
)

fragments$larynx_survfit_code <- paste(
  c(
    "nd <- data.frame(",
    "  age   = 65,",
    "  stage = factor(levels(larynx$stage))",
    ")",
    "",
    "survfit(bpph_fit, newdata = nd)"
  ),
  collapse = "\n"
)
fragments$larynx_survfit_out <- tex_lines({
  survfit(bpph_fit, newdata = nd)
})

plot_times_larynx <- seq(0, max(larynx$time), length.out = 121)
fit_bp_ml <- survfit(bpph_fit, newdata = nd, times = plot_times_larynx)
curve_labs <- paste("Stage", levels(nd$stage))
d_l <- .bp_tidy_survfit_long(fit_bp_ml, curve_labs)
p_larynx <- ggplot(d_l, aes(x = time, y = estimate, color = series)) +
  geom_line(linewidth = 0.45) +
  labs(x = "Time (years)", y = "Survival probability", color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
save_fig(p_larynx, "001_larynx_survfit.pdf", width = 7, height = 4.5)

# --- Larynx Bayes ---------------------------------------------------------
set.seed(1)
bpph_fit_bayes <- bpph(
  Surv(time, delta) ~ age + stage,
  data = larynx,
  approach = "bayes"
)

fragments$larynx_bayes_code <- paste(
  c(
    "set.seed(1)",
    "bpph_fit_bayes <- bpph(",
    "  Surv(time, delta) ~ age + stage,",
    "  data = larynx,",
    '  approach = "bayes"',
    ")",
    "",
    "survfit(bpph_fit_bayes, newdata = nd)",
    "",
    "plot_times <- seq(0, max(larynx$time), length.out = 121)",
    "sf_bayes <- survfit(",
    "  bpph_fit_bayes,",
    "  newdata = nd,",
    "  times = plot_times,",
    '  type = "log",',
    '  interval.type = "hpd"',
    ")"
  ),
  collapse = "\n"
)
fragments$larynx_bayes_survfit_out <- tex_lines({
  survfit(bpph_fit_bayes, newdata = nd)
})
fragments$larynx_bayes_print_code <- "print(bpph_fit_bayes)"
fragments$larynx_bayes_print_out <- tex_lines({
  print(bpph_fit_bayes)
})
fragments$larynx_bayes_tidy_code <- 'print(bpph_fit_bayes, what = "tidy")'
fragments$larynx_bayes_tidy_out <- tex_lines({
  print(bpph_fit_bayes, what = "tidy")
})
fragments$larynx_bayes_glance_code <- 'print(bpph_fit_bayes, what = "glance")'
fragments$larynx_bayes_glance_out <- tex_lines({
  print(bpph_fit_bayes, what = "glance")
})

# --- Larynx residuals (fig 002) -------------------------------------------
mod_cox <- coxph(
  Surv(time, delta) ~ age + stage,
  data = larynx
)

fragments$larynx_resid_code <- paste(
  c(
    "mod_cox <- coxph(",
    "  Surv(time, delta) ~ age + stage,",
    "  data = larynx",
    ")",
    "",
    "resid_df <- data.frame(",
    '  ml_bpph = residuals(bpph_fit, type = "martingale"),',
    "  bayes_bpph = residuals(bpph_fit_bayes,",
    '    type = "martingale"),',
    '  cox = residuals(mod_cox, type = "martingale")',
    ")"
  ),
  collapse = "\n"
)

resid_df <- data.frame(
  ml_bpph = residuals(bpph_fit, type = "martingale"),
  bayes_bpph = residuals(bpph_fit_bayes, type = "martingale"),
  cox = residuals(mod_cox, type = "martingale")
)

p_res1 <- ggplot(resid_df, aes(x = cox, y = ml_bpph)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Cox martingale", y = "MLE BPPH martingale", title = "(i) MLE vs Cox") +
  theme_bw()
p_res2 <- ggplot(resid_df, aes(x = cox, y = bayes_bpph)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Cox martingale", y = "Bayes BPPH martingale", title = "(ii) Bayes vs Cox") +
  theme_bw()
p_res3 <- ggplot(resid_df, aes(x = ml_bpph, y = bayes_bpph)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "MLE BPPH martingale", y = "Bayes BPPH martingale", title = "(iii) MLE vs Bayes") +
  theme_bw()
save_fig(p_res1 + p_res2 + p_res3, "002_larynx_bpxcox_residuals.pdf", width = 10, height = 3.5)

# --- Veteran ----------------------------------------------------------------
data("veteran", package = "survival")
veteran2 <- veteran[veteran$prior == 0, ]
veteran2$celltype <- factor(
  veteran2$celltype,
  levels = c("large", "adeno", "smallcell", "squamous")
)

f <- Surv(time, status) ~ karno + celltype
fit_po <- bppo(f, data = veteran2, approach = "bayes")
fit_aft <- bpaft(f, data = veteran2, approach = "bayes")

fragments$veteran_setup_code <- paste(
  c(
    'data("veteran", package = "survival")',
    "",
    "# Pettitt (1984) analytic subsample: no prior therapy",
    "veteran2 <- veteran[veteran$prior == 0, ]",
    "veteran2$celltype <- factor(",
    "  veteran2$celltype,",
    '  levels = c("large", "adeno", "smallcell", "squamous")',
    ")",
    "",
    "f <- Surv(time, status) ~ karno + celltype",
    "",
    'fit_po  <- bppo(f, data = veteran2, approach = "bayes")',
    'fit_aft <- bpaft(f, data = veteran2, approach = "bayes")'
  ),
  collapse = "\n"
)

newdata2 <- data.frame(karno = c(30, 70), celltype = "squamous")
vet_plot_times <- seq(0, max(veteran2$time), length.out = 241)
fit_bp3 <- survfit(fit_po, newdata = newdata2, times = vet_plot_times)
fit_bp4 <- survfit(fit_aft, newdata = newdata2, times = vet_plot_times)

curve_cmp <- c(
  paste0("BPAFT (Karnofsky ", c(30, 70), ")"),
  paste0("BPPO (Karnofsky ", c(30, 70), ")")
)
d_cmp <- rbind(
  .bp_tidy_survfit_long(fit_bp4, curve_cmp[1:2]),
  .bp_tidy_survfit_long(fit_bp3, curve_cmp[3:4])
)
d_cmp$series <- factor(d_cmp$series, levels = curve_cmp)
p_vet_surv <- ggplot(d_cmp, aes(x = time, y = estimate, color = series)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = series), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.45) +
  labs(x = "Time (days)", y = "Survival probability", color = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
save_fig(p_vet_surv, "003_veteran_bpaft_bppo_survival.pdf", width = 7, height = 4.5)

mod_cox_vet <- coxph(
  Surv(time, status) ~ karno + celltype,
  data = veteran2
)
sf_cox <- survival::survfit(mod_cox_vet, newdata = newdata2, times = vet_plot_times)
d_cox <- .bp_tidy_survfit_long(
  sf_cox,
  paste0("Cox PH (Karnofsky ", c(30, 70), ")")
)
d_all <- rbind(d_cmp, d_cox)
d_all$series <- factor(
  d_all$series,
  levels = c(curve_cmp, paste0("Cox PH (Karnofsky ", c(30, 70), ")"))
)
p_vet_cox <- ggplot(d_all, aes(x = time, y = estimate, color = series)) +
  geom_ribbon(
    data = subset(d_all, !grepl("^Cox", as.character(series))),
    aes(ymin = conf.low, ymax = conf.high, fill = series),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.45) +
  labs(x = "Time (days)", y = "Survival probability", color = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
save_fig(p_vet_cox, "004_veteran_bpaft_bppo_cox.pdf", width = 7.5, height = 5)

fragments$veteran_diag_code <- paste(
  c(
    "mod_cox_vet <- coxph(",
    "  Surv(time, status) ~ karno + celltype,",
    "  data = veteran2",
    ")",
    "",
    'fit_po_mle  <- bppo(f, data = veteran2, approach = "mle")',
    'fit_aft_mle <- bpaft(f, data = veteran2, approach = "mle")',
    "",
    'mr_bppo <- residuals(fit_po_mle, type = "martingale")',
    'mr_bpaft <- residuals(fit_aft_mle, type = "martingale")',
    'mr_cox <- residuals(mod_cox_vet, type = "martingale")',
    "",
    "bppo_null_fit <- bppo(",
    "  Surv(time, status) ~ 1,",
    "  data = veteran2,",
    '  approach = "mle"',
    ")",
    'mr_bppo_null <- residuals(bppo_null_fit, type = "martingale")'
  ),
  collapse = "\n"
)

fit_po_mle <- bppo(f, data = veteran2, approach = "mle")
fit_aft_mle <- bpaft(f, data = veteran2, approach = "mle")

km_bp3 <- survfit(Surv(resid(fit_po_mle, "cox-snell"), veteran2$status) ~ 1)
km_bp4 <- survfit(Surv(resid(fit_aft_mle, "cox-snell"), veteran2$status) ~ 1)

p1 <- ggplot() +
  geom_step(aes(x = km_bp3$time, y = -log(km_bp3$surv))) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x = "Residuals BPPO", y = "Cumulative hazard", title = "Cox-Snell BPPO") +
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_bw()

p2 <- ggplot() +
  geom_step(aes(x = km_bp4$time, y = -log(km_bp4$surv))) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(x = "Residuals BPAFT", y = "Cumulative hazard", title = "Cox-Snell BPAFT") +
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_bw()

mr_bppo <- residuals(fit_po_mle, type = "martingale")
mr_bpaft <- residuals(fit_aft_mle, type = "martingale")
mr_cox_vet <- residuals(mod_cox_vet, type = "martingale")

p_fm_bp3_k <- ggplot() +
  geom_point(aes(x = veteran2$karno, y = mr_bppo, color = factor(veteran2$status))) +
  geom_smooth(aes(x = veteran2$karno, y = mr_bppo), method = "loess", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-4, 2)) +
  labs(title = "BPPO: Karnofsky", x = "Karnofsky score", y = "Martingale residuals") +
  theme_bw() +
  theme(legend.position = "none")

p_fm_bp4_k <- ggplot() +
  geom_point(aes(x = veteran2$karno, y = mr_bpaft, color = factor(veteran2$status))) +
  geom_smooth(aes(x = veteran2$karno, y = mr_bpaft), method = "loess", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-4, 2)) +
  labs(title = "BPAFT: Karnofsky", x = "Karnofsky score", y = "Martingale residuals") +
  theme_bw() +
  theme(legend.position = "none")

p_fm_bp3_ct <- ggplot(data.frame(celltype = veteran2$celltype, mr = mr_bppo)) +
  geom_boxplot(aes(x = celltype, y = mr), outlier.alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "BPPO: cell type", x = "Cell type", y = "Martingale residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fm_bp4_ct <- ggplot(data.frame(celltype = veteran2$celltype, mr = mr_bpaft)) +
  geom_boxplot(aes(x = celltype, y = mr), outlier.alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "BPAFT: cell type", x = "Cell type", y = "Martingale residuals") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bppo_null_fit <- bppo(Surv(time, status) ~ 1, data = veteran2, approach = "mle")
mr_bppo_null <- residuals(bppo_null_fit, type = "martingale")

p3 <- ggplot() +
  geom_point(aes(
    x = veteran2$karno, y = mr_bppo_null,
    color = factor(veteran2$status)
  )) +
  coord_cartesian(ylim = c(-4, 2)) +
  geom_smooth(aes(x = veteran2$karno, y = mr_bppo_null),
    method = "loess", color = "black", se = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "BPPO null vs Karnofsky",
    y = "Martingale residuals", x = "Karnofsky score"
  ) +
  theme_bw() +
  theme(legend.position = "none")

mod_bp4_null <- bpaft(Surv(time, status) ~ 1, data = veteran2, approach = "mle")
mr_bpaft_null <- residuals(mod_bp4_null, type = "martingale")
p4 <- ggplot() +
  geom_point(aes(
    x = veteran2$karno, y = mr_bpaft_null,
    color = factor(veteran2$status)
  )) +
  coord_cartesian(ylim = c(-4, 2)) +
  geom_smooth(aes(x = veteran2$karno, y = mr_bpaft_null),
    method = "loess", color = "black", se = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "BPAFT null vs Karnofsky",
    y = "Martingale residuals", x = "Karnofsky score"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_mart <- (p1 | p2) /
  (p_fm_bp3_k | p_fm_bp4_k) /
  (p_fm_bp3_ct | p_fm_bp4_ct) /
  (p3 | p4)
save_fig(p_mart, "005_veteran_martingale.pdf", width = 10, height = 12)

source(file.path(pkg_root, "inst", "render-larynx-degree-comparison.R"))
render_larynx_degree_comparison(
  output = file.path(fig_dir, "006_larynx_degree_comparison.pdf")
)

mcsim_rds <- file.path(paper_dir, "bp-mcsim-results.rds")
if (file.exists(mcsim_rds)) {
  source(file.path(pkg_root, "inst", "render-bp-mcsim-abc.R"))
  render_bp_mcsim_abc(
    results_rds = mcsim_rds,
    output = file.path(fig_dir, "007_bp_mcsim_abc.pdf")
  )
} else {
  message(
    "Skipping 007_bp_mcsim_abc.pdf (missing paper/bp-mcsim-results.rds)"
  )
}

degree_csv <- file.path(paper_dir, "degree_llph_bpph.csv")
if (file.exists(degree_csv)) {
  source(file.path(pkg_root, "inst", "render-degree-llph-bpph.R"))
  render_degree_llph_bpph(
    csv_path = degree_csv,
    output = file.path(fig_dir, "008_degree_llph_bpph.pdf")
  )
} else {
  message(
    "Skipping 008_degree_llph_bpph.pdf (missing paper/degree_llph_bpph.csv)"
  )
}

mc_table_script <- file.path(pkg_root, "inst", "rebuild-mc-appendix-tables.R")
if (file.exists(mcsim_rds) || file.exists(degree_csv)) {
  status <- system2("Rscript", mc_table_script)
  if (!identical(status, 0L)) {
    stop("rebuild-mc-appendix-tables.R failed", call. = FALSE)
  }
}

# --- Persist fragments for TeX patching -----------------------------------
frag_path <- file.path(pkg_root, "paper", "tex-fragments.rds")
saveRDS(fragments, frag_path)
message("Saved verbatim fragments to ", frag_path)

rebuild_script <- file.path(pkg_root, "inst", "rebuild-tex-illustrations.R")
if (file.exists(rebuild_script) && identical(Sys.getenv("SPSURV_REBUILD_TEX"), "1")) {
  system2("Rscript", rebuild_script)
}
