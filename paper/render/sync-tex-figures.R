# Keep \\includegraphics{...} paths in paper/spsurv.TeX aligned with paper/paths.R.
#
# Usage (from package root):
#   Rscript paper/render/sync-tex-figures.R
#   Rscript paper/render/sync-tex-figures.R --verify-files

sync_tex_figure_paths <- function(
    paths,
    patch = TRUE,
    verify_files = FALSE) {
  tex_path <- paths$tex_path
  tex <- readLines(tex_path, warn = FALSE)
  figures <- paths$figures
  updated <- FALSE

  figure_for_pdf <- function(pdf) {
    for (spec in figures) {
      if (identical(spec$pdf, pdf)) {
        return(spec)
      }
    }
    id <- sub("^([0-9]{3})_.*", "\\1", pdf)
    key <- paste0("fig_", id)
    if (key %in% names(figures)) {
      return(figures[[key]])
    }
    NULL
  }

  ig_lines <- grep("\\\\includegraphics", tex)
  for (line_num in ig_lines) {
    line <- tex[[line_num]]
    m <- regexpr(
      "\\\\includegraphics(?:\\[[^]]*\\])?\\{([^}]+)\\}",
      line,
      perl = TRUE
    )
    if (m[[1L]] < 0L) {
      next
    }
    current_pdf <- regmatches(line, m)[[1L]]
    current_pdf <- sub(".*\\{([^}]+)\\}", "\\1", current_pdf)
    spec <- figure_for_pdf(current_pdf)
    if (is.null(spec)) {
      next
    }
    if (identical(current_pdf, spec$pdf)) {
      next
    }
    tex[[line_num]] <- sub(
      paste0("\\{", current_pdf, "\\}"),
      paste0("{", spec$pdf, "}"),
      line,
      fixed = TRUE
    )
    updated <- TRUE
    message(
      "Patched line ", line_num, ": ", current_pdf, " -> ", spec$pdf
    )
  }

  registry_pdfs <- vapply(figures, function(x) x$pdf, character(1))
  found_pdfs <- vapply(
    registry_pdfs,
    function(pdf) any(grepl(paste0("{", pdf, "}"), tex, fixed = TRUE)),
    logical(1)
  )
  missing_refs <- registry_pdfs[!found_pdfs]
  if (length(missing_refs)) {
    stop(
      "spsurv.TeX is missing \\includegraphics references to: ",
      paste(missing_refs, collapse = ", "),
      call. = FALSE
    )
  }

  if (verify_files) {
    missing_pdfs <- character()
    for (spec in figures) {
      if (!file.exists(spec$path)) {
        missing_pdfs <- c(missing_pdfs, spec$pdf)
      }
    }
    if (length(missing_pdfs)) {
      stop(
        "Missing figure PDF(s) under paper/: ",
        paste(missing_pdfs, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (patch && updated) {
    writeLines(tex, tex_path)
    message("Synced figure paths in ", tex_path)
  }

  invisible(updated)
}

if (identical(sys.nframe(), 0L)) {
  for (src in c("paper/paths.R", "../paths.R")) {
    if (file.exists(src)) {
      source(src, local = FALSE)
      break
    }
  }
  paths <- source_paper_paths()
  args <- commandArgs(trailingOnly = TRUE)
  verify_files <- "--verify-files" %in% args
  sync_tex_figure_paths(paths, patch = TRUE, verify_files = verify_files)
}
