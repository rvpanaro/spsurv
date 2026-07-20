# Paper assets for `spsurv.TeX`

Everything for the manuscript and its simulations lives under `paper/`.

## Directory map

```
paper/
├── spsurv.TeX                 manuscript
├── paths.R                    paths + figure registry (single source of truth)
├── README.md
├── bp-mcsim-results.rds       aggregated MC summaries (Table tab:bp-mcsim-abc)
├── degree_llph_bpph.csv       degree-sensitivity summaries (fig 008, appendix table)
├── tex-fragments.rds          cached verbatim blocks (generated)
├── *.pdf                      figures 001–008 (generated; same names as in spsurv.TeX)
├── simulation/
│   ├── monte_carlo_regression.R
│   ├── build_bp_mcsim_rds.R
│   ├── build_degree_llph_bpph_csv.R
│   ├── build-all-summaries.R    aggregate both MC summary files
│   ├── R/simulation_io.R
│   └── output/                raw replicate tables (generated)
└── render/
    ├── render-paper-figures.R   main entry point
    ├── render-tex-assets.R      illustrations 001–006 + verbatim fragments
    ├── render-mc-figures.R        simulation summaries → 007–008 PDFs + appendix tables
    ├── sync-tex-figures.R         verify/patch \\includegraphics paths in spsurv.TeX
    ├── render-bp-mcsim-abc.R
    ├── render-degree-llph-bpph.R
    ├── render-larynx-degree-comparison.R
    ├── mc-paper-summary.R
    ├── rebuild-mc-appendix-tables.R
    ├── rebuild-tex-illustrations.R
    └── tex-verbatim-wrap.R
```

## Figure registry

All PDF filenames are defined once in `paper/paths.R` (`paper_figure_specs()`).
Render scripts write `paper/<pdf>`; `spsurv.TeX` references the same bare filenames.
`sync-tex-figures.R` checks that every registry PDF appears in `\includegraphics{...}`.

## Full production workflow

All artifacts stay under `paper/`. Temp files go to `paper/logs/tmp/`.

```bash
# Overnight production (recommended): preflight + parallel sim + figures; skips PDF
./paper/run-overnight.sh

# One-shot (includes PDF if spsurv.bbl exists)
./paper/run-full-workflow.sh

# Or stepwise:
Rscript paper/run-full-workflow.R --simulation-only   # degree + MC + aggregate
Rscript paper/run-full-workflow.R --skip-simulation   # figures + PDF from saved results
Rscript paper/run-full-workflow.R --skip-pdf          # sim + figures, no PDF
Rscript paper/run-full-workflow.R --from-step=mc_regression
```

Production defaults: **R=1000**, **n ∈ {50, 100}**, **parallel MC** (snowfall),
**`SPSURV_MC_BAYES_CORES=1`** per fit (avoids CPU oversubscription), Stan defaults.

Monitor overnight run:

```bash
tail -f paper/logs/workflow-*.log
wc -l paper/simulation/output/results-*.txt
ls -la paper/degree_llph_bpph.csv paper/bp-mcsim-results.rds
```

Run manifests: `paper/logs/manifests/run-*.rds`

## Build pipeline

From the package root:

```bash
# Illustrations only (001–006)
Rscript paper/render/render-paper-figures.R --figures=illustrations

# Monte Carlo figures (007–008) when summaries already exist under paper/
Rscript paper/render/render-paper-figures.R --figures=mc

# Full Monte Carlo path (Bayes uses spbp/rstan defaults: 4 chains, iter=2000)
Rscript -e 'devtools::load_all(".", quiet=TRUE); source("paper/simulation/monte_carlo_regression.R")'
Rscript paper/simulation/build-all-summaries.R
Rscript paper/render/render-paper-figures.R --figures=mc --require

# Everything (illustrations + MC figures)
Rscript paper/render/render-paper-figures.R --figures=all --require
```

Smoke test:

```bash
Rscript paper/render/render-paper-figures.R --run-smoke-sim --figures=all
```

## Data flow

| Figure | Kind | Simulation input | Summary file | Render script |
|--------|------|------------------|--------------|---------------|
| 001–005 | illustration | — | — | `render-tex-assets.R` |
| 006 | illustration | — | — | `render-larynx-degree-comparison.R` |
| 008 | simulation | (inline MC in build script) | `degree_llph_bpph.csv` | `render-degree-llph-bpph.R` |

Table `tab:bp-mcsim-abc` (main text) and appendix table `tab:degree-llph-bpph` are patched into
`spsurv.TeX` by `rebuild-mc-appendix-tables.R` from the same two summary files.

## LaTeX

```bash
cd paper && pdflatex spsurv.TeX
```

PDFs must sit alongside `spsurv.TeX` in `paper/` (not a separate `figures/` tree).
