# Canonical paths for the paper/spsurv.TeX build.
#
# Layout (all under paper/)
#   simulation/output/          raw Monte Carlo text output (gitignored)
#   simulation/*.R              Monte Carlo runners and aggregators
#   render/*.R                    figure/table render scripts
#   bp-mcsim-results.rds          aggregated MC summaries (Table tab:bp-mcsim-abc, appendix degree table)
#   degree_llph_bpph.csv          degree-sensitivity summaries (figure 008, appendix table)
#   *.pdf                         rendered figures 001--008 (gitignored)
#   spsurv.TeX                    manuscript
#   tex-fragments.rds             cached verbatim blocks (gitignored)

paper_figure_specs <- function() {
  list(
    fig_001 = list(
      id = "001",
      pdf = "001_larynx_survfit.pdf",
      tex_label = "fig:larynx-survfit-plot",
      kind = "illustration"
    ),
    fig_002 = list(
      id = "002",
      pdf = "002_larynx_bpxcox_residuals.pdf",
      tex_label = "fig:larynx-bpxcox",
      kind = "illustration"
    ),
    fig_003 = list(
      id = "003",
      pdf = "003_veteran_bpaft_bppo_survival.pdf",
      tex_label = "fig:bpaftxbppo",
      kind = "illustration"
    ),
    fig_004 = list(
      id = "004",
      pdf = "004_veteran_bpaft_bppo_cox.pdf",
      tex_label = "fig:bpaftxbppo-cox",
      kind = "illustration"
    ),
    fig_005 = list(
      id = "005",
      pdf = "005_veteran_martingale.pdf",
      tex_label = "fig:martingale",
      kind = "illustration"
    ),
    fig_006 = list(
      id = "006",
      pdf = "006_larynx_degree_comparison.pdf",
      tex_label = "fig:larynx-survfit-plot",
      kind = "illustration"
    ),
    fig_008 = list(
      id = "008",
      pdf = "008_degree_llph_bpph.pdf",
      tex_label = "fig:degree-llph-bpph-plot",
      kind = "simulation",
      summary = "degree_csv"
    )
  )
}

paper_figures_resolved <- function(paper_dir) {
  specs <- paper_figure_specs()
  summary_paths <- list(
    mcsim_rds = file.path(paper_dir, "bp-mcsim-results.rds"),
    degree_csv = file.path(paper_dir, "degree_llph_bpph.csv")
  )
  for (nm in names(specs)) {
    specs[[nm]]$path <- file.path(paper_dir, specs[[nm]]$pdf)
    if (!is.null(specs[[nm]]$summary)) {
      specs[[nm]]$summary_path <- summary_paths[[specs[[nm]]$summary]]
    }
  }
  specs
}

paper_figure_tex_graphics <- function(spec, width = "\\linewidth") {
  sprintf("\\includegraphics[width=%s]{%s}", width, spec$pdf)
}

paper_figures_by_kind <- function(figures, kind) {
  figures[vapply(figures, function(x) identical(x$kind, kind), logical(1))]
}

spsurv_pkg_root <- function() {
  env_root <- Sys.getenv("SPSURV_ROOT", unset = NA_character_)
  if (nzchar(env_root) && file.exists(file.path(env_root, "DESCRIPTION"))) {
    return(normalizePath(env_root, winslash = "/"))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg)) {
    script_dir <- normalizePath(dirname(file_arg[[1L]]), winslash = "/")
    for (candidate in c(
      normalizePath(file.path(script_dir, "..", ".."), winslash = "/"),
      normalizePath(file.path(script_dir, ".."), winslash = "/"),
      normalizePath(".", winslash = "/")
    )) {
      if (file.exists(file.path(candidate, "DESCRIPTION"))) {
        return(candidate)
      }
    }
  }

  for (candidate in c(
    normalizePath(".", winslash = "/"),
    normalizePath("..", winslash = "/")
  )) {
    if (file.exists(file.path(candidate, "DESCRIPTION"))) {
      return(candidate)
    }
  }

  stop("Could not locate package root (set SPSURV_ROOT).", call. = FALSE)
}

spsurv_paper_paths <- function(pkg_root = spsurv_pkg_root()) {
  paper_dir <- file.path(pkg_root, "paper")
  sim_dir <- file.path(paper_dir, "simulation")
  render_dir <- file.path(paper_dir, "render")
  list(
    pkg_root = pkg_root,
    paper_dir = paper_dir,
    render_dir = render_dir,
    paths_r = file.path(paper_dir, "paths.R"),
    tex_path = file.path(paper_dir, "spsurv.TeX"),
    tex_fragments = file.path(paper_dir, "tex-fragments.rds"),
    mcsim_rds = file.path(paper_dir, "bp-mcsim-results.rds"),
    degree_csv = file.path(paper_dir, "degree_llph_bpph.csv"),
    sim_dir = sim_dir,
    sim_output_dir = file.path(sim_dir, "output"),
    figures = paper_figures_resolved(paper_dir),
    figure_path = function(filename) file.path(paper_dir, filename)
  )
}

source_paper_paths <- function() {
  if (!exists("spsurv_paper_paths", mode = "function")) {
    candidates <- c(
      file.path(Sys.getenv("SPSURV_ROOT", unset = "."), "paper", "paths.R"),
      file.path("paper", "paths.R"),
      file.path("..", "paths.R"),
      file.path("..", "..", "paths.R"),
      "paths.R"
    )
    sourced <- FALSE
    for (src in candidates) {
      if (file.exists(src)) {
        source(src, local = FALSE)
        sourced <- TRUE
        break
      }
    }
    if (!sourced || !exists("spsurv_paper_paths", mode = "function")) {
      stop("Could not find paper/paths.R", call. = FALSE)
    }
  }
  spsurv_paper_paths()
}

load_paper_paths <- source_paper_paths
