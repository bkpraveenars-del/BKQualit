\name{bq_heterogeneity}
\alias{bq_heterogeneity}
\title{Pooled and heterogeneity chi-square across families}
\usage{
bq_heterogeneity(data, family, classes, trait = NULL, ratio = c(3, 1))
}
\arguments{
  \item{data}{A data frame with one row per family.}
  \item{family}{Character; the family identifier column.}
  \item{classes}{Character vector; the columns holding counts of each phenotype class.}
  \item{trait}{Character or \code{NULL}; an optional column identifying several traits, each partitioned separately.}
  \item{ratio}{Numeric vector giving the expected ratio, default \code{c(3, 1)}.}
}
\value{
An object of class \code{bq_heterogeneity} containing the per-family goodness of fit and the pooled, heterogeneity and total chi-square components.
}
\description{
Partitions the total chi-square computed over several segregating families into a pooled component and a heterogeneity component. This is the standard way to decide why a segregation hypothesis fails.
}
\details{
Writing \eqn{\chi^2_i} for the goodness-of-fit statistic of family \eqn{i}, the total is \eqn{\sum_i \chi^2_i} on \eqn{k(c-1)} degrees of freedom, the pooled statistic is computed from the summed counts on \eqn{c-1} degrees of freedom, and the heterogeneity statistic is the difference on \eqn{(k-1)(c-1)} degrees of freedom, where \eqn{k} is the number of families and \eqn{c} the number of phenotype classes. The partition is exact.

The two components answer different questions. A significant pooled chi-square with a non-significant heterogeneity means the families behave alike but the hypothesised ratio is wrong. A significant heterogeneity means the families do not behave alike, so pooling them is not justified and the pooled test should not be interpreted.
}
\seealso{\code{\link{bq_segregate}}}
\examples{
d <- bq_data("families")
res <- bq_heterogeneity(d, family = "family",
                        classes = c("dominant", "recessive"),
                        trait = "trait", ratio = c(3, 1))
res
bq_plot(res)
}
