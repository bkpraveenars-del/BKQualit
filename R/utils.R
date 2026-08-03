## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

## Convert an hclust object into plotting segments without depending on
## ggdendro. Heights are drawn on the y axis, leaves in hclust order on x.
.dendro_segments <- function(hc) {
  n <- length(hc$order)
  pos <- numeric(nrow(hc$merge))          # x position of each merge node
  leafpos <- match(seq_len(n), hc$order)  # x position of each original leaf
  segs <- list()
  for (i in seq_len(nrow(hc$merge))) {
    kids <- hc$merge[i, ]
    xy <- vapply(kids, function(k) {
      if (k < 0) c(leafpos[-k], 0) else c(pos[k], hc$height[k])
    }, numeric(2))
    xl <- xy[1, 1]; yl <- xy[2, 1]
    xr <- xy[1, 2]; yr <- xy[2, 2]
    h <- hc$height[i]
    segs[[length(segs) + 1L]] <- data.frame(
      x    = c(xl, xr, xl),
      y    = c(yl, yr, h),
      xend = c(xl, xr, xr),
      yend = c(h,  h,  h))
    pos[i] <- (xl + xr) / 2

  }
  out <- do.call(rbind, segs)
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------------------
## Exported helpers
## ---------------------------------------------------------------------------

#' Load a bundled example dataset
#'
#' Reads one of the demonstration datasets shipped with the package.
#'
#' @param name One of \code{"f2"} (single-cross segregation counts),
#'   \code{"families"} (segregation counts for many families),
#'   \code{"linkage"} (two-locus counts) or \code{"dus"} (qualitative
#'   descriptors scored on a set of accessions).
#' @return A data frame.
#' @examples
#' str(bq_data("f2"))
#' head(bq_data("dus"))
#' @export
bq_data <- function(name = c("f2", "families", "linkage", "dus")) {
  name <- match.arg(name)
  file <- switch(name,
    f2       = "f2_segregation.csv",
    families = "family_segregation.csv",
    linkage  = "linkage_data.csv",
    dus      = "dus_descriptors.csv")
  path <- system.file("extdata", file, package = "BKQualit")
  if (!nzchar(path))
    stop("Dataset '", name, "' not found. Is BKQualit installed correctly?")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Plot a BKQualit result
#'
#' A single plotting verb for every analysis in the package. The method chosen
#' depends on the class of \code{x}, and several results accept a \code{type}
#' argument selecting among alternative views.
#'
#' @param x An object returned by one of the \code{bq_*} analysis functions.
#' @param ... Arguments passed to the method, typically \code{type}.
#' @return A \code{ggplot} object.
#' @examples
#' res <- bq_segregate(bq_data("f2"), phenotype = "phenotype",
#'                     count = "count", trait = "trait")
#' bq_plot(res)
#' @export
bq_plot <- function(x, ...) UseMethod("bq_plot")

#' @export
bq_plot.default <- function(x, ...) {
  stop("bq_plot() has no method for objects of class '",
       paste(class(x), collapse = "/"), "'.")
}

#' Save every figure from a BKQualit result
#'
#' Convenience wrapper that writes all available views of a result to disk as
#' PNG files.
#'
#' @param x An object returned by one of the \code{bq_*} analysis functions.
#' @param path Directory to write into. Defaults to a temporary directory.
#' @param prefix File-name prefix.
#' @param width,height,dpi Passed to \code{ggplot2::ggsave}.
#' @return Character vector of the files written, invisibly.
#' @examples
#' res <- bq_segregate(bq_data("f2"), phenotype = "phenotype",
#'                     count = "count", trait = "trait")
#' bq_save(res, path = tempdir())
#' @export
bq_save <- function(x, path = tempdir(), prefix = NULL,
                    width = 8, height = 5.5, dpi = 300) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  cls <- class(x)[1]
  if (is.null(prefix)) prefix <- cls
  types <- switch(cls,
    bq_segregate    = NA_character_,
    bq_heterogeneity = NA_character_,
    bq_linkage      = c("profile", "estimates"),
    bq_diversity    = c("indices", "frequencies"),
    bq_mca          = c("biplot", "dendrogram", "scree"),
    NA_character_)
  out <- character(0)
  for (ty in types) {
    p <- if (is.na(ty)) bq_plot(x) else bq_plot(x, type = ty)
    f <- file.path(path, paste0(prefix,
      if (is.na(ty)) "" else paste0("_", ty), ".png"))
    ggplot2::ggsave(f, p, width = width, height = height, dpi = dpi)
    out <- c(out, f)
  }
  message("Wrote ", length(out), " figure(s) to ", normalizePath(path))
  invisible(out)
}
