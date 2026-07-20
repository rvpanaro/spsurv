# Rebuild paper/spsurv.TeX illustration verbatim blocks from paper/tex-fragments.rds
for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
pkg_root <- paths$pkg_root

fragments <- readRDS(paths$tex_fragments)
tex_path <- paths$tex_path
tex <- readLines(tex_path, warn = FALSE)

source(file.path(paths$render_dir, "tex-verbatim-wrap.R"))

wrap_code <- function(code, width = 65L) {
  lines <- tex_verbatim_wrap_lines(strsplit(code, "\n", fixed = TRUE)[[1L]], width = width)
  c("{\\scriptsize\\begin{verbatim}", lines, "\\end{verbatim}}")
}

wrap_out <- function(out, width = 65L) {
  lines <- tex_verbatim_wrap_lines(strsplit(out, "\n", fixed = TRUE)[[1L]], width = width)
  c("{\\scriptsize\\begin{verbatim}", lines, "\\end{verbatim}}")
}

illustration_start <- grep("Model fit and coefficients", tex, fixed = TRUE)[1L]
residual_start <- grep("\\\\subsubsection\\{Residual diagnostics\\}", tex, fixed = FALSE)[1L]
if (is.na(illustration_start) || is.na(residual_start)) {
  stop("Could not locate illustration section markers.")
}

before <- tex[seq_len(illustration_start - 1L)]
after <- tex[seq.int(residual_start, length(tex))]

middle <- c(
  "\\subsubsection{Model fit and coefficients}\\label{model-fit-and-coefficients}",
  "",
  "The BPPH model described in Section~\\ref{sec:models} was fitted to the",
  "\\texttt{larynx} data using maximum likelihood estimation. Categorical",
  "covariates such as \\texttt{stage} are converted to factors in the data",
  "before model fitting rather than through \\texttt{factor()} in the formula,",
  "so coefficient labels appear as \\texttt{stage2}, \\texttt{stage3}, and",
  "\\texttt{stage4}. Numeric output in the examples below is rounded to two",
  "decimal places for display. The corresponding hazard-ratio estimates are",
  "reported below.",
  "",
  wrap_code(fragments$larynx_setup),
  "",
  wrap_code(fragments$larynx_tidy_code),
  "",
  wrap_out(fragments$larynx_tidy_out),
  "",
  wrap_code(fragments$larynx_glance_code),
  "",
  wrap_out(fragments$larynx_glance_out),
  "",
  "\\needspace{6\\baselineskip}",
  "",
  "The estimated hazard ratio for Stage IV versus Stage I at the same age",
  "is 6.05, with 95\\% limits 2.61",
  "and 14.02 from \\texttt{tidy(bpph\\_fit, conf.int = TRUE, exponentiate = TRUE)},",
  "indicating higher mortality for Stage IV than for Stage",
  "I. The \\texttt{tidy()} and \\texttt{glance()} methods from \\pkg{generics} provide",
  "broom-style summaries: \\texttt{print(..., what = \"tidy\")}",
  "returns one row per coefficient with exponentiated intervals, and",
  "\\texttt{print(..., what = \"glance\")} returns a one-row table of model-level",
  "statistics including sample size, log-likelihood, likelihood-ratio test, and",
  "information criteria. For maximum-likelihood fits, \\texttt{glance()} and \\texttt{AIC()}",
  "count both regression and Bernstein (\\texttt{gamma}) parameters",
  "(\\texttt{.spbp\\_nparams()}); the likelihood-ratio \\texttt{df} in",
  "\\texttt{summary()} remains the number of regression terms only.",
  "For models fitted by maximum likelihood,",
  "\\texttt{print()} with the default \\texttt{what = \"summary\"} and \\texttt{summary()}",
  "report coefficient estimates, Wald 95\\% confidence limits on the log scale",
  "(\\(2.5\\%\\) and \\(97.5\\%\\)), standard errors, \\(z\\) statistics, and \\(p\\)-values in a",
  "compact layout, followed by exponentiated confidence intervals on the same",
  "percentiles and model-level log-likelihood and information criteria at the foot of the summary.",
  "",
  wrap_code(fragments$larynx_print_code),
  "",
  wrap_out(fragments$larynx_print_out),
  "",
  "\\subsubsection{Survival prediction}\\label{survival-prediction}",
  "",
  "Predicted survival curves for \\texttt{spbp} objects are obtained with",
  "\\texttt{survfit()} in \\pkg{spsurv} (\\texttt{survfit.spbp}). The method",
  "follows the familiar \\texttt{survfit} generic interface while replacing",
  "\\pkg{survival} \\citep{survival:2000} point estimates with Bernstein-polynomial",
  "predictions; a reference Cox object from \\pkg{survival} supplies only the",
  "time grid and \\texttt{survfit} object structure. Predicted curves for each",
  "stage at the reference age are obtained as follows. For maximum-likelihood",
  "fits on small samples such as \\texttt{larynx}, the default Bernstein degree",
  "$m=\\lceil\\sqrt{n}\\rceil$ can leave the $\\boldsymbol{\\gamma}$ information",
  "block ill-conditioned (Section~\\ref{sec:delta}); \\texttt{survfit.spbp} still",
  "computes delta-method bands but issues a single warning that they may be",
  "unreliable and advising a lower",
  "\\texttt{degree}. Fitting, \\texttt{summary()}, and \\texttt{vcov()} do not emit",
  "this warning. Passing",
  "\\texttt{tidy = TRUE} returns a long data frame suitable for \\pkg{ggplot2}.",
  "By default, \\texttt{survfit.spbp} evaluates curves on a dense time grid; users may pass",
  "an explicit \\texttt{times} vector for plotting.",
  "",
  wrap_code(fragments$larynx_survfit_code),
  "",
  wrap_out(fragments$larynx_survfit_out),
  "",
  "Figure~\\ref{fig:larynx-survfit-plot} shows predicted survival across",
  "stages at age 65 years (rounded from the sample mean age",
  "64.61 years), computed with \\texttt{survfit.spbp} in \\pkg{spsurv} and",
  "displayed with \\pkg{ggplot2}. Survival decreases with",
  "increasing stage. Where a curve crosses 0.5, that time is an estimated",
  "median survival on the time scale of the data; curves that remain above 0.5 on the plotted",
  "follow-up have no median crossing visible on the grid shown.",
  "Panel~(a) uses the default Bernstein degree",
  "$m=\\lceil\\sqrt{n}\\rceil=10$ on the 90-patient \\texttt{larynx} sample;",
  "the $\\boldsymbol{\\gamma}$ information block is then ill-conditioned",
  "(Section~\\ref{sec:delta}), so \\texttt{survfit.spbp} reports pointwise",
  "delta-method bands together with a warning that they may be unreliable",
  "(printed at this step only).",
  "Panel~(b) contrasts the same covariate profiles under three fits:",
  "$m=3$ maximum likelihood (stable $\\boldsymbol{\\gamma}$ information and",
  "95\\% delta-method bands without an instability warning), default",
  "$m=10$ maximum likelihood (ill-conditioned $\\boldsymbol{\\gamma}$ block and",
  "unreliable delta-method bands with a \\texttt{survfit()} warning), and",
  "default $m=10$ Bayesian (posterior HPD bands at the same degree, without",
  "the delta method).",
  "Lowering $m$ trades baseline flexibility for identifiability on small",
  "samples and is one adjustment suggested by the \\texttt{survfit()} warning;",
  "panel~(b) shows that refitting with \\texttt{approach = \"bayes\"} at the",
  "same degree is the other main route for curve-wise uncertainty. The next",
  "subsection gives posterior summaries and code for the Bayesian fit.",
  "",
  "\\begin{figure}[!ht]",
  "\\centering",
  "\\begin{minipage}{0.48\\linewidth}",
  "\\centering",
  paste0(paper_figure_tex_graphics(paths$figures$fig_001), "\\\\[0.4em]"),
  "{\\footnotesize (a) Default $m=10$ MLE with delta-method bands.}",
  "\\end{minipage}\\hfill",
  "\\begin{minipage}{0.48\\linewidth}",
  "\\centering",
  paste0(paper_figure_tex_graphics(paths$figures$fig_006), "\\\\[0.4em]"),
  "{\\footnotesize (b) $m=3$ MLE, $m=10$ MLE, and $m=10$ Bayes (top to bottom).}",
  "\\end{minipage}",
  "\\caption{Estimated survival curves from maximum-likelihood BPPH fits for a",
  "patient aged 65 years in the \\texttt{larynx} data (\\pkg{KMsurv}).",
  "Time is in years. Curves use \\texttt{survfit.spbp} on a dense time grid",
  "and \\pkg{ggplot2}. Panel~(a): default Bernstein degree with",
  "ill-conditioned $\\boldsymbol{\\gamma}$ information---delta-method",
  "intervals are reported with an instability warning from \\texttt{survfit()}.",
  "Panel~(b): three-way comparison at age 65 years---$m=3$ MLE with stable",
  "delta-method bands, default $m=10$ MLE with an instability warning on",
  "delta-method bands, and default $m=10$ Bayes with posterior HPD bands at",
  "the same polynomial degree. Lower curves indicate poorer prognosis; stage",
  "ordering is monotone from highest survival (Stage~I) to lowest (Stage~IV).}",
  "\\label{fig:larynx-survfit-plot}",
  "\\end{figure}",
  "",
  "\\subsubsection{Bayesian analysis}\\label{bayesian-analysis}",
  "",
  "As an alternative to lowering Bernstein \\texttt{degree} when maximum-likelihood",
  "delta-method bands are unreliable (Figure~\\ref{fig:larynx-survfit-plot},",
  "panel~(b)), the same model can be refitted with",
  "\\texttt{approach\\ =\\ \"bayes\"}; the posterior",
  "summaries agree closely with the maximum likelihood estimates, as the",
  "priors are weakly informative by default, using the prior specification",
  "of Section \\ref{sec:priors}. Unlike the maximum-likelihood fit above,",
  "the Bayesian \\texttt{survfit.spbp} summary returns finite credible limits",
  "for median survival when the posterior supports them. Curve-wise bands",
  "use posterior draws of \\((\\boldsymbol{\\beta},\\boldsymbol{\\gamma})\\) with",
  "\\texttt{interval.type = \"hpd\"} by default and optional monotone enforcement",
  "over time. For the Bayesian fits, we ran 4 chains (the \\texttt{spbp.default()}",
  "default) without overriding \\texttt{iter},",
  "\\texttt{warmup}, \\texttt{thin}, \\texttt{adapt\\_delta}, or",
  "\\texttt{max\\_treedepth}, and the stored fitted-model \\texttt{call} lists",
  "\\texttt{refresh = 0} (Section~\\ref{sec:mcmc}).",
  "",
  wrap_code(fragments$larynx_bayes_code),
  "",
  wrap_out(fragments$larynx_bayes_survfit_out),
  "",
  "\\needspace{6\\baselineskip}",
  "",
  "Posterior summaries provide highest posterior density (HPD) intervals",
  "and model comparison via deviance information criterion and widely",
  "applicable information criterion (WAIC) when available, with the same PH",
  "interpretation. For Bayesian fits (\\texttt{approach\\ =\\ \"bayes\"}), \\texttt{print()}",
  "lists posterior means, HPD credible intervals on the log and exponentiated",
  "scales (\\(2.5\\%\\) and \\(97.5\\%\\)), and posterior standard deviations for",
  "regression coefficients, and reports deviance information",
  "criterion \\citep{spiegelhalter2002bayesian}, WAIC \\citep{Vehtari:2017, watanabe2013widely}, and log pseudo-marginal likelihood",
  "\\citep{geisser1979predictive, Ibrahim:2014} when the posterior",
  "log-likelihood is available. \\texttt{summary()} gives the same HPD intervals",
  "as \\texttt{print()} for regression and exponentiated coefficients; the",
  "exponentiated intervals also appear in \\texttt{print(..., what = \"tidy\")} and",
  "\\texttt{tidy(..., conf.int = TRUE, exponentiate = TRUE)}, while",
  "\\texttt{print(..., what = \"glance\")} and \\texttt{glance()} collect WAIC, DIC,",
  "and LPML in one row.",
  "",
  wrap_code(fragments$larynx_bayes_print_code),
  "",
  wrap_out(fragments$larynx_bayes_print_out),
  "",
  wrap_code(fragments$larynx_bayes_tidy_code),
  "",
  wrap_out(fragments$larynx_bayes_tidy_out),
  "",
  wrap_code(fragments$larynx_bayes_glance_code),
  "",
  wrap_out(fragments$larynx_bayes_glance_out),
  "",
  "Extractors such as \\texttt{tidy()}, \\texttt{glance()}, \\texttt{coef()},",
  "\\texttt{model.matrix()}, and \\texttt{credint()} follow",
  "the conventions of standard fitted objects and integrate with the",
  "\\pkg{generics} / broom workflow. When comparing models,",
  "posterior means and HPD limits from \\texttt{tidy()} can be read alongside",
  "WAIC and log pseudo-marginal likelihood from \\texttt{glance()}, for example when",
  "choosing between PO and AFT formulations that both fit the data reasonably well.",
  "",
  wrap_code(fragments$larynx_resid_code)
)

# Veteran blocks live in the tail; rebuild them inline once.
veteran_alt_idx <- grep(
  "\\\\subsubsection\\{Alternatives to proportional hazards\\}",
  after,
  fixed = FALSE
)[1L]
veteran_setup_idx <- grep(
  'data("veteran", package = "survival")',
  after,
  fixed = TRUE
)[1L]
if (is.na(veteran_alt_idx) || is.na(veteran_setup_idx)) {
  stop("Could not locate veteran illustration markers.")
}
veteran_setup_start <- veteran_setup_idx
repeat {
  prev <- veteran_setup_start - 1L
  if (prev < 1L) {
    break
  }
  if (grepl("\\\\begin\\{verbatim\\}|\\\\footnotesize|\\\\scriptsize", after[prev])) {
    veteran_setup_start <- prev
  } else {
    break
  }
}
veteran_setup_end <- veteran_setup_idx +
  grep("\\end{verbatim}", after[veteran_setup_idx:length(after)], fixed = TRUE)[1L] - 1L

veteran_diag_idx <- grep("mr_bppo_null <- residuals(bppo_null_fit", after, fixed = TRUE)[1L]
veteran_diag_start <- veteran_diag_idx
while (veteran_diag_start > 1L && !grepl("begin\\{verbatim\\}", after[veteran_diag_start])) {
  veteran_diag_start <- veteran_diag_start - 1L
}
veteran_diag_end <- veteran_diag_start +
  grep("\\end{verbatim}", after[veteran_diag_start:length(after)], fixed = TRUE)[1L] - 1L

after <- c(
  after[seq_len(veteran_setup_start - 1L)],
  wrap_code(fragments$veteran_setup_code),
  after[seq.int(veteran_setup_end + 1L, veteran_diag_start - 1L)],
  wrap_code(fragments$veteran_diag_code),
  after[seq.int(veteran_diag_end + 1L, length(after))]
)

tex <- c(before, middle, after)
writeLines(tex, tex_path)
message("Rebuilt illustration section in ", tex_path)
