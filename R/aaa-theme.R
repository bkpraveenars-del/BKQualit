#' BKQualit colour palettes
#'
#' Palettes for qualitative-trait figures. \code{"phenotype"} is a vivid
#' categorical set for phenotype classes; \code{"fit"} runs red to green for
#' goodness of fit; \code{"linkage"} is a sequential ramp for recombination
#' fraction; \code{"descriptor"} is a soft categorical set for descriptor
#' states.
#'
#' @param name Palette name: \code{"phenotype"} (default), \code{"fit"},
#'   \code{"linkage"}, \code{"descriptor"}, \code{"contrast"}.
#' @param n Number of colours; \code{NULL} returns the anchor palette.
#' @param reverse Logical; reverse the palette.
#' @return A character vector of hex colours.
#' @examples
#' bq_palette("phenotype")
#' bq_palette("linkage", 6)
#' @export
bq_palette <- function(name = c("phenotype", "fit", "linkage", "descriptor",
                                "contrast"), n = NULL, reverse = FALSE) {
  name <- match.arg(name)
  anchors <- list(
    phenotype  = c("#6A2C91", "#C2185B", "#F2A03D", "#4CAF50", "#1E88A8",
                   "#8D6E63", "#546E7A"),
    fit        = c("#B3261E", "#E8792B", "#F2C14E", "#8CC63F", "#1B7340"),
    linkage    = c("#0B3D91", "#3D7BC4", "#9EC5E8", "#F5C6C6", "#C62828"),
    descriptor = c("#A8DADC", "#457B9D", "#E9C46A", "#F4A261", "#E76F51",
                   "#B5838D", "#6D6875"),
    contrast   = c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51",
                   "#6A4C93", "#118AB2")
  )
  pal <- anchors[[name]]
  if (reverse) pal <- rev(pal)
  if (is.null(n)) return(pal)
  if (n <= length(pal) && name %in% c("phenotype", "descriptor", "contrast"))
    return(pal[seq_len(n)])
  grDevices::colorRampPalette(pal)(n)
}

#' A clean ggplot2 theme for BKQualit figures
#'
#' @param base_size Base font size in points.
#' @param base_family Font family.
#' @param grid Logical; draw light major grid lines.
#' @return A \code{ggplot2} theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(factor(cyl))) + geom_bar() + theme_bq()
#' @export
theme_bq <- function(base_size = 12, base_family = "", grid = TRUE) {
  ink <- "#1A1626"; panel <- "#FCFBFD"; subtle <- "#6C6880"
  gl <- if (grid) ggplot2::element_line(colour = "#E8E4EE", linewidth = 0.3)
        else ggplot2::element_blank()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = panel, colour = NA),
      panel.grid.major = gl,
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(fill = NA, colour = "#DCD6E4",
                                               linewidth = 0.4),
      axis.title = ggplot2::element_text(colour = ink, face = "bold"),
      axis.text  = ggplot2::element_text(colour = subtle),
      plot.title = ggplot2::element_text(colour = ink, face = "bold",
                     size = base_size * 1.3,
                     margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(colour = subtle,
                     size = base_size * 0.95,
                     margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(colour = subtle,
                     size = base_size * 0.75),
      legend.title = ggplot2::element_text(colour = ink, face = "bold"),
      legend.text  = ggplot2::element_text(colour = subtle),
      strip.background = ggplot2::element_rect(fill = "#6A2C91", colour = NA),
      strip.text = ggplot2::element_text(colour = "white", face = "bold",
                     margin = ggplot2::margin(3, 3, 3, 3))
    )
}

#' Discrete BKQualit colour and fill scales
#'
#' @param pal_name Palette name (see \code{\link{bq_palette}}).
#' @param reverse Logical; reverse the palette.
#' @param ... Passed to \code{ggplot2::discrete_scale}.
#' @return A ggplot2 scale.
#' @rdname bq_scales
#' @export
scale_colour_bq <- function(pal_name = "phenotype", reverse = FALSE, ...) {
  f <- function(n) bq_palette(pal_name, n = n, reverse = reverse)
  ggplot2::discrete_scale("colour", palette = f, ...)
}

#' @rdname bq_scales
#' @export
scale_color_bq <- scale_colour_bq

#' @rdname bq_scales
#' @export
scale_fill_bq <- function(pal_name = "phenotype", reverse = FALSE, ...) {
  f <- function(n) bq_palette(pal_name, n = n, reverse = reverse)
  ggplot2::discrete_scale("fill", palette = f, ...)
}
