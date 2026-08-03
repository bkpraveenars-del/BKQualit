\name{scale_colour_bq}
\alias{scale_colour_bq}
\alias{scale_color_bq}
\alias{scale_fill_bq}
\title{Discrete colour and fill scales for BKQualit palettes}
\usage{
scale_colour_bq(pal_name = "phenotype", reverse = FALSE, ...)

scale_fill_bq(pal_name = "phenotype", reverse = FALSE, ...)
}
\arguments{
  \item{pal_name}{Palette name; see \code{\link{bq_palette}}. This argument is called \code{pal_name} rather than \code{name} so that the \code{ggplot2} legend title, which is also called \code{name}, can still be passed through \code{...}.}
  \item{reverse}{Logical; reverse the colour order.}
  \item{...}{Passed to \code{ggplot2::discrete_scale}, for example \code{name} to set the legend title.}
}
\value{A \code{ggplot2} scale object.}
\description{
Discrete scales built on the package palettes. \code{scale_color_bq} is an alias of \code{scale_colour_bq}.
}
\seealso{\code{\link{bq_palette}}, \code{\link{theme_bq}}}
\examples{
library(ggplot2)
ggplot(bq_data("dus"), aes(growth_habit, fill = seed_shape)) +
  geom_bar() +
  scale_fill_bq("descriptor", name = "Seed shape") +
  theme_bq()
}
