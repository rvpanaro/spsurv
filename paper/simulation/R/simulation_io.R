#' Locate the newest timestamped Monte Carlo output file.
#'
#' @param prefix File stem prefix, e.g. \code{"results"}.
#' @param output_dir Directory containing \code{prefix-*.txt} files.
#' @return Normalized path, or \code{NULL} if none found.
find_latest_output <- function(prefix, output_dir = "output") {
  if (!dir.exists(output_dir)) {
    return(NULL)
  }
  pat <- paste0("^", prefix, "-[0-9]{8}-[0-9]{6}-[0-9]+\\.txt$")
  files <- list.files(output_dir, pattern = pat, full.names = TRUE)
  if (!length(files)) {
    return(NULL)
  }
  files[which.max(file.info(files)$mtime)]
}

#' Read a Monte Carlo results table with canonical column names.
read_results_table <- function(path) {
  d <- utils::read.table(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )
  if (ncol(d) != 13L) {
    stop(
      "Expected 13 columns in results file, found ", ncol(d), ": ", path,
      call. = FALSE
    )
  }
  names(d) <- c(
    "nsize", "par", "real", "estimate", "se", "RB", "lwr", "upr", "CP",
    "gdist", "approach", "model", "rep"
  )
  d$nsize <- as.integer(d$nsize)
  d$rep <- as.integer(d$rep)
  d$CP <- as.logical(d$CP)
  d
}

#' Read a Monte Carlo censoring table.
read_censoring_table <- function(path) {
  d <- utils::read.table(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )
  if (ncol(d) != 6L) {
    stop(
      "Expected 6 columns in censoring file, found ", ncol(d), ": ", path,
      call. = FALSE
    )
  }
  names(d) <- c(
    "nsize", "gdist", "approach", "model", "rep", "event_proportion"
  )
  d$nsize <- as.integer(d$nsize)
  d$rep <- as.integer(d$rep)
  d$event_proportion <- as.numeric(d$event_proportion)
  d
}
