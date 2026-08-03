#' BKQualit: Analysis of Qualitative Traits in Plant Breeding
#'
#' Tools for the genetic analysis of qualitative (discretely inherited)
#' characters: testing observed segregation against Mendelian expectations,
#' partitioning chi-square across families into pooled and heterogeneity
#' components, estimating linkage by maximum likelihood, and summarising and
#' clustering qualitative descriptors of the kind recorded on a distinctness,
#' uniformity and stability (DUS) form.
#'
#' @section Analysis functions:
#' \describe{
#'   \item{\code{\link{bq_segregate}}}{Fit observed counts to Mendelian ratios;
#'     ranks a library of classical ratios and names the gene action implied.}
#'   \item{\code{\link{bq_heterogeneity}}}{Partition the total chi-square over
#'     families into pooled and heterogeneity components.}
#'   \item{\code{\link{bq_linkage}}}{Estimate the recombination fraction by
#'     maximum likelihood with a logarithm of odds (LOD) score and a
#'     likelihood-ratio confidence interval.}
#'   \item{\code{\link{bq_diversity}}}{Shannon, Simpson and evenness indices
#'     for qualitative descriptors, with discriminating power.}
#'   \item{\code{\link{bq_mca}}}{Multiple correspondence analysis with the
#'     Benzecri adjustment, plus Gower-distance clustering of accessions.}
#' }
#'
#' @section Support:
#' \code{\link{bq_data}} loads the bundled example datasets,
#' \code{\link{bq_plot}} draws every result, \code{\link{bq_save}} writes all
#' figures for a result to disk, and \code{\link{bq_palette}},
#' \code{\link{theme_bq}}, \code{\link{scale_colour_bq}} and
#' \code{\link{scale_fill_bq}} provide the shared visual style.
#'
#' @section Design:
#' Every analysis function returns a classed list with a \code{print} method
#' that states the conclusion in words as well as numbers, and a
#' \code{bq_plot} method. Nothing is written to the user's file space unless
#' \code{\link{bq_save}} is called explicitly.
#'
#' @author Praveen Kumar B. K. \email{bkpraveenars@@gmail.com}
#' @references
#' Mather K (1951). \emph{The Measurement of Linkage in Heredity}, 2nd edition.
#' Methuen, London.
#'
#' Allard RW (1956). Formulas and tables to facilitate the calculation of
#' recombination values in heredity. \emph{Hilgardia} \strong{24}, 235--278.
#' \doi{10.3733/hilg.v24n10p235}
#'
#' Gower JC (1971). A general coefficient of similarity and some of its
#' properties. \emph{Biometrics} \strong{27}, 857--871. \doi{10.2307/2528823}
#' @keywords internal
"_PACKAGE"
