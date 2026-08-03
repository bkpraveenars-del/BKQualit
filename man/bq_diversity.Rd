\name{bq_diversity}
\alias{bq_diversity}
\title{Diversity of qualitative descriptors}
\usage{
bq_diversity(data, descriptors, accession = NULL, group = NULL)
}
\arguments{
  \item{data}{A data frame with one row per accession.}
  \item{descriptors}{Character vector naming the categorical columns.}
  \item{accession}{Character or \code{NULL}; the accession identifier column.}
  \item{group}{Character or \code{NULL}; an optional grouping column, giving diversity within each group as well as overall.}
}
\value{
An object of class \code{bq_diversity} containing per-descriptor diversity indices, the frequency of every state, and, if \code{group} is supplied, the same indices within each group.
}
\description{
Summarises a set of categorical descriptors, of the kind recorded on a distinctness, uniformity and stability (DUS) form, into per-descriptor diversity indices, and identifies which descriptors actually discriminate among accessions.
}
\details{
Three complementary indices are reported for each descriptor. The Shannon-Weaver index \eqn{H' = -\sum p_i \log p_i} is sensitive to rare states; the Simpson index \eqn{1 - \sum p_i^2} is dominated by common states; and evenness \eqn{J' = H'/\log s} rescales Shannon by the number of observed states \eqn{s} so that descriptors with different numbers of states can be compared directly. A descriptor with a single observed state is monomorphic, contributes nothing to discrimination, and is flagged.

A descriptor is useful only if it separates accessions. The discriminating power reported here is the probability that two randomly chosen accessions differ in that descriptor, which equals the Simpson index; descriptors are ranked on it so that the least informative can be dropped from a characterisation form.
}
\seealso{\code{\link{bq_mca}}}
\examples{
d <- bq_data("dus")
res <- bq_diversity(d,
  descriptors = c("growth_habit", "leaf_pubescence", "flower_colour",
                  "seed_colour", "seed_shape"),
  accession = "accession")
res
bq_plot(res)
bq_plot(res, type = "frequencies")
}
