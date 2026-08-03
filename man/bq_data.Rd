\name{bq_data}
\alias{bq_data}
\title{Load a bundled example dataset}
\usage{
bq_data(name = c("f2", "families", "linkage", "dus"))
}
\arguments{
  \item{name}{One of \code{"f2"} (single-cross segregation counts), \code{"families"} (segregation counts for many families), \code{"linkage"} (two-locus counts) or \code{"dus"} (qualitative descriptors scored on a set of accessions).}
}
\value{A data frame.}
\description{
Reads one of the demonstration datasets shipped with the package. The datasets are simulated from known parameters and exist for documentation, teaching and regression testing.
}
\examples{
str(bq_data("f2"))
head(bq_data("dus"))
}
