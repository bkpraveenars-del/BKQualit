\name{bq_linkage}
\alias{bq_linkage}
\title{Estimate genetic linkage by maximum likelihood}
\usage{
bq_linkage(data, pair, classes, design = NULL, phase = NULL, conf = 0.95)
}
\arguments{
  \item{data}{A data frame with one row per locus pair.}
  \item{pair}{Character; a column naming the locus pair.}
  \item{classes}{Character vector of length four giving the count columns, in the order A_B_, A_bb, aaB_, aabb.}
  \item{design}{Character or \code{NULL}; a column holding \code{"F2"} or \code{"backcross"} for each pair. If \code{NULL}, \code{"F2"} is assumed.}
  \item{phase}{Character or \code{NULL}; a column holding \code{"coupling"} or \code{"repulsion"}. If \code{NULL}, coupling is assumed.}
  \item{conf}{Confidence level (default 0.95).}
}
\value{
An object of class \code{bq_linkage} containing the estimated recombination fraction, its confidence interval, the logarithm of odds (LOD) score, the map distance and a test against independent assortment, together with the likelihood profile for plotting.
}
\description{
Estimates the recombination fraction between two loci from two-locus second filial generation (F2) or backcross counts, and tests it against independent assortment.
}
\details{
For a backcross the recombination fraction is the proportion of recombinant progeny, a binomial parameter. For an F2 in coupling phase the four phenotype classes occur in the proportions \eqn{(2+\theta)/4}, \eqn{(1-\theta)/4}, \eqn{(1-\theta)/4} and \eqn{\theta/4} with \eqn{\theta=(1-r)^2}; in repulsion phase the outer and inner proportions are exchanged. The recombination fraction is obtained by maximising the multinomial likelihood numerically rather than from the older product-ratio approximation.

The LOD score compares the fitted model with independent assortment, and a value above 3 is the conventional threshold for declaring linkage. The confidence interval is obtained by likelihood-ratio inversion, retaining every value of the recombination fraction whose log-likelihood lies within 1.921 of the maximum; this is preferred to a symmetric standard-error interval because the likelihood is skewed near the boundary. Map distance is reported as 100 times the recombination fraction and is not corrected for double crossovers.
}
\references{
Allard RW (1956). Formulas and tables to facilitate the calculation of recombination values in heredity. \emph{Hilgardia} \strong{24}, 235--278. \doi{10.3733/hilg.v24n10p235}
}
\examples{
d <- bq_data("linkage")
res <- bq_linkage(d, pair = "pair",
                  classes = c("AB", "Ab", "aB", "ab"),
                  design = "design", phase = "phase")
res
bq_plot(res)
bq_plot(res, type = "estimates")
}
