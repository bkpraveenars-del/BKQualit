\name{bq_save}
\alias{bq_save}
\title{Save every figure from a BKQualit result}
\usage{
bq_save(x, path = tempdir(), prefix = NULL, width = 8, height = 5.5,
        dpi = 300)
}
\arguments{
  \item{x}{An object returned by one of the \code{bq_*} analysis functions.}
  \item{path}{Directory to write into. Defaults to a temporary directory.}
  \item{prefix}{File-name prefix. Defaults to the class of \code{x}.}
  \item{width, height, dpi}{Passed to \code{ggplot2::ggsave}.}
}
\value{Character vector of the files written, invisibly.}
\description{
Convenience wrapper that writes all available views of a result to disk as PNG files. Nothing is written unless this function is called.
}
\examples{
res <- bq_segregate(bq_data("f2"), phenotype = "phenotype",
                    count = "count", trait = "trait")
bq_save(res, path = tempdir())
}
