#' Wrap lines for LaTeX verbatim blocks (default 65 columns).
#' @noRd
tex_verbatim_wrap_line <- function(line, width = 65L) {
  if (!nzchar(line) || nchar(line) <= width) {
    return(line)
  }

  prefix <- ""
  content <- line
  if (startsWith(line, "#> ")) {
    prefix <- "#> "
    content <- substr(line, 4L, nchar(line))
    content_width <- width - nchar(prefix)
  } else {
    content_width <- width
  }

  if (content_width < 8L) {
    stop("verbatim width is too small after accounting for prefixes", call. = FALSE)
  }

  chunks <- character()
  remaining <- content
  while (nchar(remaining) > content_width) {
    segment <- substr(remaining, 1L, content_width)
    brk <- gregexpr("[ ,;]", segment, perl = TRUE)[[1L]]
    brk <- brk[brk > 0L]
    if (!length(brk) || max(brk) < content_width %/% 2L) {
      cut <- content_width
    } else {
      cut <- max(brk)
    }
    chunks <- c(chunks, substr(remaining, 1L, cut))
    remaining <- sub("^\\s*", "", substr(remaining, cut + 1L, nchar(remaining)), perl = TRUE)
  }
  if (nzchar(remaining)) {
    chunks <- c(chunks, remaining)
  }

  paste0(prefix, chunks)
}

#' @noRd
tex_verbatim_wrap_lines <- function(x, width = 65L) {
  if (length(x) == 1L && grepl("\n", x, fixed = TRUE)) {
    lines <- strsplit(x, "\n", fixed = TRUE)[[1L]]
  } else {
    lines <- as.character(x)
  }
  unlist(lapply(lines, tex_verbatim_wrap_line, width = width), use.names = FALSE)
}
