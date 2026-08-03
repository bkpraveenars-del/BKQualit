#' Estimate genetic linkage by maximum likelihood
#'
#' Estimates the recombination fraction between two loci from two-locus
#' second filial generation (F2) or backcross counts, and tests it against
#' independent assortment.
#'
#' For a backcross the recombination fraction is simply the proportion of
#' recombinant progeny, a binomial parameter. For an F2 in coupling phase the
#' four phenotype classes occur in the proportions
#' \deqn{\frac{2+\theta}{4},\ \frac{1-\theta}{4},\ \frac{1-\theta}{4},\
#'   \frac{\theta}{4}, \qquad \theta = (1-r)^2,}
#' and in repulsion phase the first and last proportions are exchanged with the
#' middle pair. The recombination fraction is obtained by maximising the
#' multinomial likelihood numerically rather than from the older product-ratio
#' approximation.
#'
#' Two quantities summarise the evidence. The logarithm of odds (LOD) score
#' compares the fitted model with independent assortment, and a value above 3
#' is the conventional threshold for declaring linkage. The confidence interval
#' is obtained by likelihood-ratio inversion, retaining every value of the
#' recombination fraction whose log-likelihood lies within 1.921 of the maximum;
#' this is preferred to a symmetric standard-error interval because the
#' likelihood is skewed near the boundary.
#'
#' @param data A data frame with one row per locus pair.
#' @param pair Character; a column naming the locus pair.
#' @param classes Character vector of length four giving the count columns, in
#'   the order A_B_, A_bb, aaB_, aabb.
#' @param design Character or \code{NULL}; a column holding \code{"F2"} or
#'   \code{"backcross"} for each pair. If \code{NULL}, \code{"F2"} is assumed.
#' @param phase Character or \code{NULL}; a column holding \code{"coupling"} or
#'   \code{"repulsion"}. If \code{NULL}, coupling is assumed.
#' @param conf Confidence level (default 0.95).
#' @return An object of class \code{bq_linkage}.
#' @references
#' Allard RW (1956). Formulas and tables to facilitate the calculation of
#' recombination values in heredity. \emph{Hilgardia} \strong{24}, 235--278.
#' \doi{10.3733/hilg.v24n10p235}
#' @examples
#' d <- bq_data("linkage")
#' res <- bq_linkage(d, pair = "pair",
#'                   classes = c("AB", "Ab", "aB", "ab"),
#'                   design = "design", phase = "phase")
#' res
#' bq_plot(res)
#' @export
bq_linkage <- function(data, pair, classes, design = NULL, phase = NULL,
                       conf = 0.95) {
  stopifnot(all(c(pair, classes, design, phase) %in% names(data)))
  if (length(classes) != 4)
    stop("'classes' must name exactly four count columns.")
  d <- data[stats::complete.cases(data[c(pair, classes)]), , drop = FALSE]
  cut <- stats::qchisq(conf, 1) / 2          # 1.920729 at 95 per cent

  ## expected proportions given r
  props <- function(r, des, ph) {
    if (des == "backcross") {
      if (ph == "repulsion") c(r, 1 - r, 1 - r, r) / 2
      else c(1 - r, r, r, 1 - r) / 2
    } else {
      th <- (1 - r)^2
      if (ph == "repulsion") {
        rr <- r^2
        c(2 + rr, 1 - rr, 1 - rr, rr) / 4
      } else {
        c(2 + th, 1 - th, 1 - th, th) / 4
      }
    }
  }
  loglik <- function(r, obs, des, ph) {
    p <- props(r, des, ph)
    if (any(p <= 0)) return(-Inf)
    sum(obs * log(p))
  }

  res <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    obs <- as.numeric(d[i, classes])
    n <- sum(obs)
    des <- if (is.null(design)) "F2" else as.character(d[[design]][i])
    ph  <- if (is.null(phase)) "coupling" else as.character(d[[phase]][i])

    opt <- stats::optimize(function(r) -loglik(r, obs, des, ph),
                           interval = c(1e-6, 0.5 - 1e-6))
    rhat <- opt$minimum
    ll1 <- -opt$objective
    ll0 <- loglik(0.5 - 1e-9, obs, des, ph)
    lod <- (ll1 - ll0) / log(10)

    grid <- seq(1e-6, 0.5 - 1e-6, length.out = 4000)
    lls <- vapply(grid, loglik, numeric(1), obs = obs, des = des, ph = ph)
    keep <- grid[lls >= ll1 - cut]
    lo <- if (length(keep)) min(keep) else NA_real_
    hi <- if (length(keep)) max(keep) else NA_real_

    ## test against independent assortment
    p0 <- props(0.5, des, ph); e0 <- n * p0
    chi <- sum((obs - e0)^2 / e0)
    pv <- stats::pchisq(chi, 3, lower.tail = FALSE)

    data.frame(pair = as.character(d[[pair]][i]), design = des, phase = ph,
      n = n, r = rhat, lower = lo, upper = hi, LOD = lod,
      map_distance_cM = 100 * rhat,
      chisq_independence = chi, p_independence = pv,
      linked = lod >= 3 & rhat < 0.5,
      row.names = NULL, stringsAsFactors = FALSE)
  }))
  res$Sig <- .stars(res$p_independence)

  ## likelihood profile for plotting
  prof <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    obs <- as.numeric(d[i, classes])
    des <- if (is.null(design)) "F2" else as.character(d[[design]][i])
    ph  <- if (is.null(phase)) "coupling" else as.character(d[[phase]][i])
    grid <- seq(1e-4, 0.5 - 1e-4, length.out = 300)
    lls <- vapply(grid, loglik, numeric(1), obs = obs, des = des, ph = ph)
    data.frame(pair = as.character(d[[pair]][i]), r = grid,
               lod = (lls - loglik(0.5 - 1e-9, obs, des, ph)) / log(10),
               row.names = NULL, stringsAsFactors = FALSE)
  }))

  structure(list(estimates = res, profile = prof, conf = conf),
            class = "bq_linkage")
}

#' @export
print.bq_linkage <- function(x, ...) {
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("  BKQualit | Linkage estimated by maximum likelihood\n")
  cat(strrep("=", 74), "\n\n", sep = "")
  e <- x$estimates
  print(data.frame(
    Pair = e$pair, Design = e$design, n = e$n,
    r = round(e$r, 4),
    CI = paste0("[", round(e$lower, 4), ", ", round(e$upper, 4), "]"),
    LOD = round(e$LOD, 3),
    cM = round(e$map_distance_cM, 2),
    Linked = ifelse(e$linked, "YES", "no"),
    stringsAsFactors = FALSE), row.names = FALSE)
  cat("\nTest against independent assortment (9:3:3:1 or 1:1:1:1):\n\n")
  print(data.frame(Pair = e$pair, ChiSq = round(e$chisq_independence, 3),
                   df = 3, p = format.pval(e$p_independence, digits = 4),
                   Sig = e$Sig, stringsAsFactors = FALSE), row.names = FALSE)
  cat("\nr is the recombination fraction: 0 = completely linked,",
      "0.5 = independent.\n")
  cat("LOD >= 3 is the conventional threshold for declaring linkage.\n")
  cat("Confidence intervals are obtained by likelihood-ratio inversion,\n")
  cat("which respects the skewness of the likelihood near the boundary.\n")
  cat("Map distance in centimorgans is reported as 100r (uncorrected).\n\n")
  invisible(x)
}

#' @export
bq_plot.bq_linkage <- function(x, type = c("profile", "estimates"), ...) {
  type <- match.arg(type)
  e <- x$estimates
  if (type == "profile") {
    ggplot2::ggplot(x$profile, ggplot2::aes(.data$r, .data$lod,
                                            colour = .data$pair)) +
      ggplot2::geom_hline(yintercept = 3, linetype = "dashed",
                          colour = "#B3261E", linewidth = 0.5) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(data = e,
        ggplot2::aes(x = .data$r, y = .data$LOD, colour = .data$pair),
        size = 3, inherit.aes = FALSE) +
      scale_colour_bq("phenotype", name = "Locus pair") +
      ggplot2::labs(title = "Linkage likelihood profiles",
        subtitle = "LOD against recombination fraction; dashed line is the LOD = 3 threshold",
        x = "recombination fraction (r)", y = "LOD score",
        caption = "BKQualit") +
      theme_bq()
  } else {
    ggplot2::ggplot(e, ggplot2::aes(x = stats::reorder(.data$pair, .data$r),
                                    y = .data$r, colour = .data$linked)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                          ymax = .data$upper),
                             width = 0.16, linewidth = 0.8) +
      ggplot2::geom_point(size = 3.6) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                          colour = "#6C6880") +
      ggplot2::coord_flip() +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#1B7340",
                                              `FALSE` = "#B3261E"),
        labels = c(`TRUE` = "linked", `FALSE` = "independent"), name = NULL) +
      ggplot2::labs(title = "Recombination fractions with confidence intervals",
        subtitle = paste0(round(100 * x$conf),
          "% likelihood-ratio intervals; 0.5 = independent assortment"),
        x = NULL, y = "recombination fraction (r)", caption = "BKQualit") +
      theme_bq()
  }
}
