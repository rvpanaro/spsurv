# Regenerate figures/005_veteran_martingale.pdf.
# Usage (from package root):
#   Rscript inst/render-veteran-martingale.R

pkg_root <- Sys.getenv("SPSURV_ROOT", unset = normalizePath("..", winslash = "/"))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  pkg_root <- normalizePath(".", winslash = "/")
}
setwd(pkg_root)

suppressPackageStartupMessages({
  devtools::load_all(pkg_root, quiet = TRUE)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

fig_dir <- file.path(pkg_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

data("veteran", package = "survival")
veteran2 <- veteran[veteran$prior == 0, ]
veteran2$celltype <- factor(
  veteran2$celltype,
  levels = c("large", "adeno", "smallcell", "squamous")
)
f <- Surv(time, status) ~ karno + celltype
fit_po_mle <- bppo(f, data = veteran2, approach = "mle")
fit_aft_mle <- bpaft(f, data = veteran2, approach = "mle")

km_bp3 <- survfit(Surv(resid(fit_po_mle, "cox-snell"), veteran2$status) ~ 1)
km_bp4 <- survfit(Surv(resid(fit_aft_mle, "cox-snell"), veteran2$status) ~ 1)

p1 <- ggplot() +
  geom_step(aes(x = km_bp3$time, y = -log(km_bp3$surv)), linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(title = "Cox-Snell BPPO", x = "Residuals BPPO", y = "Cumulative hazard") +
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_bw()

p2 <- ggplot() +
  geom_step(aes(x = km_bp4$time, y = -log(km_bp4$surv)), linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(title = "Cox-Snell BPAFT", x = "Residuals BPAFT", y = "Cumulative hazard") +
  geom_abline(slope = 1, intercept = 0, col = "blue") +
  theme_bw()

mr_bppo <- residuals(fit_po_mle, type = "martingale")
mr_bpaft <- residuals(fit_aft_mle, type = "martingale")

p_fm_bp3_k <- ggplot() +
  geom_point(aes(x = veteran2$karno, y = mr_bppo, color = factor(veteran2$status)), size = 1.8, alpha = 0.8) +
  geom_smooth(aes(x = veteran2$karno, y = mr_bppo), method = "loess", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-4, 2)) +
  labs(title = "BPPO: Karnofsky", x = "Karnofsky score", y = "Martingale residuals") +
  theme_bw() +
  theme(legend.position = "none")

p_fm_bp4_k <- ggplot() +
  geom_point(aes(x = veteran2$karno, y = mr_bpaft, color = factor(veteran2$status)), size = 1.8, alpha = 0.8) +
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
mod_bp4_null <- bpaft(Surv(time, status) ~ 1, data = veteran2, approach = "mle")
mr_bppo_null <- residuals(bppo_null_fit, type = "martingale")
mr_bpaft_null <- residuals(mod_bp4_null, type = "martingale")

p3 <- ggplot() +
  geom_point(aes(x = veteran2$karno, y = mr_bppo_null, color = factor(veteran2$status)), size = 1.8, alpha = 0.8) +
  coord_cartesian(ylim = c(-4, 2)) +
  geom_smooth(aes(x = veteran2$karno, y = mr_bppo_null), method = "loess", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "BPPO null vs Karnofsky", y = "Martingale residuals", x = "Karnofsky score") +
  theme_bw() +
  theme(legend.position = "none")

p4 <- ggplot() +
  geom_point(aes(x = veteran2$karno, y = mr_bpaft_null, color = factor(veteran2$status)), size = 1.8, alpha = 0.8) +
  coord_cartesian(ylim = c(-4, 2)) +
  geom_smooth(aes(x = veteran2$karno, y = mr_bpaft_null), method = "loess", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "BPAFT null vs Karnofsky", y = "Martingale residuals", x = "Karnofsky score") +
  theme_bw() +
  theme(legend.position = "none")

p_mart <- (p1 | p2) / (p_fm_bp3_k | p_fm_bp4_k) / (p_fm_bp3_ct | p_fm_bp4_ct) / (p3 | p4)

out_fig <- file.path(fig_dir, "005_veteran_martingale.pdf")

grDevices::pdf(out_fig, width = 10, height = 12, onefile = TRUE)
print(p_mart)
grDevices::dev.off()

message("Wrote ", out_fig)
