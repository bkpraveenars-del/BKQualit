\name{bq_plot}
\alias{bq_plot}
\alias{bq_plot.default}
\title{Plot a BKQualit result}
\usage{
bq_plot(x, ...)
}
\arguments{
  \item{x}{An object returned by one of the \code{bq_*} analysis functions.}
  \item{...}{Arguments passed to the method, typically \code{type}.}
}
\value{A \code{ggplot} object.}
\description{
A single plotting verb for every analysis in the package. The method chosen depends on the class of \code{x}.
}
\details{
Several results accept a \code{type} argument selecting among alternative views: \code{bq_linkage} offers \code{"profile"} and \code{"estimates"}; \code{bq_diversity} offers \code{"indices"} and \code{"frequencies"}; \code{bq_mca} offers \code{"biplot"}, \code{"dendrogram"} and \code{"scree"}.
}
\seealso{\code{\link{bq_save}}}
\examples{
res <- bq_segregate(bq_data("f2"), phenotype = "phenotype",
                    count = "count", trait = "trait")
bq_plot(res)
}
