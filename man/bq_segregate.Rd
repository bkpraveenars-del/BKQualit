\name{bq_segregate}
\alias{bq_segregate}
\title{Test observed segregation against Mendelian ratios}
\usage{
bq_segregate(data = NULL, phenotype = NULL, count = NULL,
             trait = NULL, ratio = NULL, counts = NULL,
             simulate = TRUE, B = 10000, seed = 1L)
}
\arguments{
  \item{data}{A data frame in long form, one row per phenotype class.}
  \item{phenotype}{Character; the phenotype-class column.}
  \item{count}{Character; the observed-count column.}
  \item{trait}{Character or \code{NULL}; a column identifying several traits, each tested separately.}
  \item{ratio}{Numeric vector giving a specific expected ratio. If \code{NULL}, every ratio in the library with the right number of classes is fitted and ranked.}
  \item{counts}{Named numeric vector, an alternative to supplying \code{data}.}
  \item{simulate}{Logical; also run a Monte Carlo exact test (default \code{TRUE}).}
  \item{B}{Monte Carlo replicates (default 10000).}
  \item{seed}{Random seed for the Monte Carlo test.}
}
\value{
An object of class \code{bq_segregate} holding the fitted ratios ranked by goodness of fit, the observed and expected counts, and the gene action implied by the best-fitting ratio.
}
\description{
Fits observed phenotype counts to the classical Mendelian expectations. When no ratio is stated, a library of standard mono-, di- and trihybrid ratios with the matching number of classes is fitted and ranked, and the gene action implied by the best fit is named.
}
\details{
The Pearson chi-square is reported with the Yates continuity correction applied only when there is a single degree of freedom, which is where the correction is defined; applying it more widely is a common error that makes tests conservative. A Monte Carlo exact test is also run, because the chi-square approximation is unreliable when an expected count falls below five, as it often does in the rarest class of a 63:1 or 9:3:3:1 table.

A high p value does not prove the hypothesis. It means the data are compatible with it. Where two ratios fit almost equally well, the printed output says so rather than declaring a winner.
}
\references{
Mather K (1951). \emph{The Measurement of Linkage in Heredity}, 2nd edition. Methuen, London.
}
\seealso{\code{\link{bq_heterogeneity}}, \code{\link{bq_plot}}}
\examples{
d <- bq_data("f2")
res <- bq_segregate(d, phenotype = "phenotype", count = "count",
                    trait = "trait")
res
bq_plot(res)

## a single trait against a stated hypothesis
bq_segregate(counts = c(purple = 312, white = 104), ratio = c(3, 1))
}
