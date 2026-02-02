#!/usr/bin/env Rscript

# Self-contained test runner for an R package using testthat:
# - Installs missing helper packages
# - Installs package deps (including Suggests)
# - Loads the package
# - Runs tests/testthat.R (or testthat::test_dir as fallback)
# - Exits non-zero if tests fail

args <- commandArgs(trailingOnly = TRUE)

# Usage:
#   Rscript run_tests.R [path_to_pkg]
pkg_path <- if (length(args) >= 1) args[1] else "."
pkg_path <- normalizePath(pkg_path, winslash = "/", mustWork = TRUE)

desc_path <- file.path(pkg_path, "DESCRIPTION")
if (!file.exists(desc_path)) {
  cat("ERROR: No DESCRIPTION found at:", desc_path, "\n")
  quit(status = 2)
}

cat("==> Package path:", pkg_path, "\n")

# ---- helper: install CRAN packages if missing
ensure_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_pkgs(c("remotes", "testthat", "pkgload", "desc"))

# ---- install dependencies (Imports + Suggests)
cat("\n==> Installing dependencies (Imports + Suggests)...\n")
tryCatch(
  remotes::install_deps(pkg_path, dependencies = TRUE, upgrade = "never"),
  error = function(e) {
    cat("ERROR installing dependencies:\n", conditionMessage(e), "\n")
    quit(status = 3)
  }
)

# ---- load package (without attaching)
pkg_name <- desc::desc(file = desc_path)$get("Package")
cat("\n==> Loading package:", pkg_name, "\n")

tryCatch(
  pkgload::load_all(pkg_path, export_all = FALSE, quiet = TRUE),
  error = function(e) {
    cat("ERROR loading package:\n", conditionMessage(e), "\n")
    quit(status = 4)
  }
)

# ---- run tests
cat("\n==> Running tests...\n")

testthat_r <- file.path(pkg_path, "tests", "testthat.R")
test_dir   <- file.path(pkg_path, "tests", "testthat")

status <- 0

tryCatch(
  {
    if (file.exists(testthat_r)) {
      # This is the canonical entry point created by usethis::use_testthat()
      sys.source(testthat_r, envir = new.env(parent = globalenv()))
    } else if (dir.exists(test_dir)) {
      # Fallback
      testthat::test_dir(test_dir, reporter = "summary")
    } else {
      cat("No tests found at tests/testthat.R or tests/testthat/\n")
      quit(status = 5)
    }
  },
  error = function(e) {
    cat("\n==> TESTS FAILED\n")
    cat(conditionMessage(e), "\n")
    status <<- 1
  }
)

if (status == 0) {
  cat("\n==> ALL TESTS PASSED\n")
} else {
  cat("\n==> SOME TESTS FAILED\n")
}

quit(status = status)
