\name{bq_mca}
\alias{bq_mca}
\title{Multiple correspondence analysis and Gower clustering of accessions}
\usage{
bq_mca(data, descriptors, accession = NULL, numeric_vars = NULL,
       k = 3, ndim = 2)
}
\arguments{
  \item{data}{A data frame with one row per accession.}
  \item{descriptors}{Character vector naming the categorical columns.}
  \item{accession}{Character or \code{NULL}; the accession identifier column.}
  \item{numeric_vars}{Character vector or \code{NULL}; numeric columns to add to the Gower distance. They do not enter the correspondence analysis.}
  \item{k}{Number of clusters to cut the dendrogram into (default 3).}
  \item{ndim}{Number of dimensions to retain (default 2).}
}
\value{
An object of class \code{bq_mca} containing the principal inertias with the Benzecri adjustment, row and column coordinates, the Gower distance matrix, the clustering, and a modal profile of each cluster.
}
\description{
Places accessions and descriptor states in a common low-dimensional space by multiple correspondence analysis (MCA), and groups the accessions by hierarchical clustering of the Gower distance.
}
\details{
MCA is computed from the singular value decomposition of the standardised indicator matrix. Writing \eqn{Z} for the \eqn{n \times J} indicator of \eqn{Q} descriptors, \eqn{P = Z/(nQ)}, and \eqn{r}, \eqn{c} for its row and column margins, the decomposition is applied to \eqn{D_r^{-1/2}(P - rc^{\top})D_c^{-1/2}}. Row and column principal coordinates then satisfy the transition relation \eqn{F = (Z/Q)G_{std}}, so an accession sits at the centroid of the descriptor states it carries, which is what makes an MCA biplot directly readable.

The raw principal inertias of an indicator MCA understate the structure badly, because coding a descriptor into several columns injects inertia that carries no information. The Benzecri adjustment is therefore also reported: axes with inertia above \eqn{1/Q} are rescaled to \eqn{((Q/(Q-1))(\lambda - 1/Q))^2}. Interpret the adjusted percentages, not the raw ones.

The Gower distance handles nominal and numeric descriptors together: nominal descriptors contribute 0 when two accessions match and 1 when they differ, numeric descriptors contribute the absolute difference divided by the range, and the contributions are averaged. Clustering uses Ward's minimum-variance criterion on that distance.
}
\references{
Gower JC (1971). A general coefficient of similarity and some of its properties. \emph{Biometrics} \strong{27}, 857--871. \doi{10.2307/2528823}
}
\seealso{\code{\link{bq_diversity}}}
\examples{
d <- bq_data("dus")
res <- bq_mca(d,
  descriptors = c("growth_habit", "leaf_pubescence", "flower_colour",
                  "seed_colour", "seed_shape", "pod_curvature"),
  accession = "accession", numeric_vars = "plant_height_cm", k = 3)
res
bq_plot(res)
bq_plot(res, type = "dendrogram")
bq_plot(res, type = "scree")
}
