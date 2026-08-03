\name{theme_bq}
\alias{theme_bq}
\title{The BKQualit ggplot2 theme}
\usage{
theme_bq(base_size = 12, base_family = "", grid = TRUE)
}
\arguments{
  \item{base_size}{Base font size in points.}
  \item{base_family}{Base font family.}
  \item{grid}{Logical; draw the background grid.}
}
\value{A \code{ggplot2} theme object.}
\description{
A clean theme with a soft background, muted grid and a coloured title, applied by every plotting method in the package. It can be added to any \code{ggplot} of your own.
}
\examples{
library(ggplot2)
ggplot(bq_data("dus"), aes(growth_habit)) +
  geom_bar(fill = bq_palette("phenotype", 1)) +
  theme_bq()
}
