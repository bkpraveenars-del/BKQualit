\name{bq_palette}
\alias{bq_palette}
\title{Colour palettes for qualitative-trait figures}
\usage{
bq_palette(name = c("phenotype", "fit", "linkage", "descriptor", "contrast"),
           n = NULL, reverse = FALSE)
}
\arguments{
  \item{name}{Palette name. \code{"phenotype"} is a vivid categorical set for phenotype classes; \code{"fit"} runs red to green for goodness of fit; \code{"linkage"} is a sequential ramp for recombination fraction; \code{"descriptor"} is a soft categorical set for descriptor states; \code{"contrast"} is a high-separation set for clusters.}
  \item{n}{Number of colours. If \code{NULL}, the whole palette is returned.}
  \item{reverse}{Logical; reverse the colour order.}
}
\value{A character vector of hexadecimal colours.}
\description{
Returns the colours used throughout the package. Categorical palettes are returned unchanged when enough colours are available and interpolated otherwise; sequential palettes are always interpolated.
}
\examples{
bq_palette("phenotype")
bq_palette("linkage", 6)
}
