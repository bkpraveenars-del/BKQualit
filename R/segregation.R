# ---- internal helpers -------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",
  ifelse(p < 0.10,  ".", "ns")))))
}

## The classical Mendelian expectations, indexed by number of phenotype classes.
.ratio_library <- function() {
  list(
    "3:1"              = c(3, 1),
    "1:1"              = c(1, 1),
    "1:2:1"            = c(1, 2, 1),
    "1:1:1:1"          = c(1, 1, 1, 1),
    "9:3:3:1"          = c(9, 3, 3, 1),
    "9:7"              = c(9, 7),
    "9:6:1"            = c(9, 6, 1),
    "9:3:4"            = c(9, 3, 4),
    "12:3:1"           = c(12, 3, 1),
    "13:3"             = c(13, 3),
    "15:1"             = c(15, 1),
    "27:9:9:9:3:3:3:1" = c(27, 9, 9, 9, 3, 3, 3, 1),
    "1:2:1:2:4:2:1:2:1" = c(1, 2, 1, 2, 4, 2, 1, 2, 1)
  )
}

## Gene-action interpretation of the common modified dihybrid ratios.
.ratio_meaning <- function(nm) {
  switch(nm,
    "3:1"       = "monohybrid, complete dominance",
    "1:1"       = "testcross / backcross segregation",
    "1:2:1"     = "monohybrid, incomplete dominance (codominance)",
    "1:1:1:1"   = "dihybrid testcross",
    "9:3:3:1"   = "dihybrid, two independent genes, complete dominance",
    "9:7"       = "duplicate recessive epistasis (complementary genes)",
    "9:6:1"     = "duplicate genes with cumulative effect",
    "9:3:4"     = "recessive epistasis",
    "12:3:1"    = "dominant epistasis",
    "13:3"      = "dominant suppression epistasis",
    "15:1"      = "duplicate dominant epistasis",
    "27:9:9:9:3:3:3:1" = "trihybrid, three independent genes",
    "1:2:1:2:4:2:1:2:1" = "dihybrid, incomplete dominance at both loci",
    "user-supplied expectation")
}

#' Test segregation against Mendelian expectations
#'
#' Compares observed phenotype counts with the ratio expected under a genetic
#' hypothesis. If no ratio is supplied, every classical Mendelian expectation
#' with the right number of classes is fitted and the candidates are ranked by
#' goodness of fit, together with the gene action each ratio implies.
#'
#' The chi-square goodness-of-fit statistic is reported with, and without, the
#' Yates continuity correction; the correction is applied only when there is one
#' degree of freedom, where it is appropriate. Because the chi-square
#' approximation is unreliable when expected counts are small, a Monte Carlo
#' exact test is also reported, and a warning is issued when any expected count
#' falls below five.
#'
#' @param data A data frame of counts, or \code{NULL} if \code{counts} is given
#'   directly.
#' @param phenotype Character; the phenotype column.
#' @param count Character; the count column.
#' @param trait Character or \code{NULL}; optional column identifying several
#'   traits, each tested separately.
#' @param ratio Numeric vector giving the expected ratio (for example
#'   \code{c(3, 1)}), or \code{NULL} to fit all classical ratios.
#' @param counts Optional named numeric vector of counts, used instead of
#'   \code{data}.
#' @param simulate Logical; compute a Monte Carlo exact p value.
#' @param B Number of Monte Carlo replicates (default 10000).
#' @param seed Random seed for the Monte Carlo test.
#' @return An object of class \code{bq_segregate}.
#' @references
#' Mather K (1951). \emph{The Measurement of Linkage in Heredity}, 2nd edition.
#' Methuen, London.
#' @examples
#' d <- bq_data("f2")
#' res <- bq_segregate(d, phenotype = "phenotype", count = "count",
#'                     trait = "trait")
#' res
#' bq_plot(res)
#'
#' ## a single trait against a stated hypothesis
#' bq_segregate(counts = c(purple = 312, white = 104), ratio = c(3, 1))
#' @export
bq_segregate <- function(data = NULL, phenotype = NULL, count = NULL,
                         trait = NULL, ratio = NULL, counts = NULL,
                         simulate = TRUE, B = 10000, seed = 1L) {
  if (is.null(data) && is.null(counts))
    stop("Supply either 'data' with phenotype/count columns, or 'counts'.")

  if (!is.null(counts)) {
    grp <- list(single = data.frame(
      phenotype = names(counts) %||% paste0("class", seq_along(counts)),
      count = as.numeric(counts), stringsAsFactors = FALSE))
  } else {
    stopifnot(all(c(phenotype, count) %in% names(data)))
    d <- data[stats::complete.cases(data[c(phenotype, count, trait)]), ]
    key <- if (is.null(trait)) rep("single", nrow(d)) else as.character(d[[trait]])
    grp <- split(data.frame(phenotype = as.character(d[[phenotype]]),
                            count = as.numeric(d[[count]]),
                            stringsAsFactors = FALSE), key)
  }

  ## Setting a seed makes the Monte Carlo test reproducible, but silently
  ## resetting the caller's random-number stream would corrupt any simulation
  ## this function is called from, so the previous state is restored on exit.
  if (simulate) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
              add = TRUE)
    }
    set.seed(seed)
  }
  lib <- .ratio_library()
  fit_one <- function(obs, r) {
    r <- as.numeric(r); p <- r / sum(r)
    n <- sum(obs); e <- n * p
    df <- length(obs) - 1
    chi <- sum((obs - e)^2 / e)
    chi_y <- if (df == 1) sum((abs(obs - e) - 0.5)^2 / e) else NA_real_
    pv <- stats::pchisq(chi, df, lower.tail = FALSE)
    pv_y <- if (df == 1) stats::pchisq(chi_y, df, lower.tail = FALSE) else NA_real_
    mc <- NA_real_
    if (simulate && n > 0) {
      sim <- stats::rmultinom(B, n, p)
      cs <- colSums((sim - e)^2 / e)
      mc <- (sum(cs >= chi) + 1) / (B + 1)
    }
    list(chi = chi, df = df, p = pv, chi_yates = chi_y, p_yates = pv_y,
         p_mc = mc, expected = e, min_exp = min(e))
  }

  out <- list(); details <- list()
  for (g in names(grp)) {
    dg <- grp[[g]]
    dg <- dg[order(-dg$count), , drop = FALSE]
    obs <- dg$count; k <- length(obs)

    if (!is.null(ratio)) {
      if (length(ratio) != k)
        stop("ratio has ", length(ratio), " parts but ", k,
             " phenotype classes were found for '", g, "'.")
      cands <- list(paste(ratio, collapse = ":"))
      names(cands) <- paste(ratio, collapse = ":")
      cands[[1]] <- ratio
    } else {
      cands <- lib[vapply(lib, length, integer(1)) == k]
      if (!length(cands))
        stop("No classical ratio has ", k, " classes; supply 'ratio'.")
    }

    tab <- do.call(rbind, lapply(names(cands), function(nm) {
      f <- fit_one(obs, cands[[nm]])
      data.frame(trait = g, ratio = nm, n = sum(obs),
                 chisq = f$chi, df = f$df, p = f$p,
                 chisq_yates = f$chi_yates, p_yates = f$p_yates,
                 p_montecarlo = f$p_mc, min_expected = f$min_exp,
                 gene_action = .ratio_meaning(nm),
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
    tab <- tab[order(-tab$p), , drop = FALSE]
    tab$Sig <- .stars(tab$p)
    tab$verdict <- ifelse(tab$p > 0.05, "fits", "rejected")
    rownames(tab) <- NULL
    out[[g]] <- tab

    best <- cands[[tab$ratio[1]]]
    e <- sum(obs) * as.numeric(best) / sum(as.numeric(best))
    details[[g]] <- data.frame(trait = g, phenotype = dg$phenotype,
      observed = obs, expected = e, deviation = obs - e,
      contribution = (obs - e)^2 / e,
      row.names = NULL, stringsAsFactors = FALSE)
  }

  fits <- do.call(rbind, out); rownames(fits) <- NULL
  best <- do.call(rbind, lapply(out, function(z) z[1, , drop = FALSE]))
  rownames(best) <- NULL
  obs_tab <- do.call(rbind, details); rownames(obs_tab) <- NULL

  if (any(fits$min_expected < 5))
    warning("Some expected counts are below 5; prefer the Monte Carlo p value.")

  structure(list(fits = fits, best = best, observed = obs_tab,
                 simulate = simulate, B = B, user_ratio = ratio),
            class = "bq_segregate")
}

#' @export
print.bq_segregate <- function(x, ...) {
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("  BKQualit | Segregation analysis against Mendelian expectations\n")
  cat(strrep("=", 74), "\n\n", sep = "")
  b <- x$best
  cat("BEST-FITTING RATIO for each trait:\n\n")
  tb <- data.frame(Trait = b$trait, n = b$n, Ratio = b$ratio,
                   ChiSq = round(b$chisq, 4), df = b$df,
                   p = round(b$p, 4), Verdict = b$verdict,
                   stringsAsFactors = FALSE)
  if (x$simulate) tb$p_MC <- round(b$p_montecarlo, 4)
  print(tb, row.names = FALSE)
  cat("\nImplied gene action:\n")
  for (i in seq_len(nrow(b)))
    cat("  ", b$trait[i], ": ", b$ratio[i], " - ", b$gene_action[i], "\n", sep = "")

  cat("\nObserved vs expected (under the best ratio):\n\n")
  o <- x$observed
  o$expected <- round(o$expected, 2); o$deviation <- round(o$deviation, 2)
  o$contribution <- round(o$contribution, 3)
  print(o, row.names = FALSE)

  if (is.null(x$user_ratio) && nrow(x$fits) > nrow(b)) {
    cat("\nAll candidate ratios tested (ranked by fit):\n\n")
    f <- x$fits
    print(data.frame(Trait = f$trait, Ratio = f$ratio,
                     ChiSq = round(f$chisq, 3), df = f$df,
                     p = round(f$p, 4), Verdict = f$verdict,
                     stringsAsFactors = FALSE), row.names = FALSE)
  }
  cat("\np > 0.05 means the observed counts are CONSISTENT with that ratio.\n")
  cat("Yates correction is applied only at 1 df, where it is appropriate.\n")
  if (any(x$fits$min_expected < 5))
    cat("Some expected counts < 5: rely on the Monte Carlo p value.\n")
  cat("\n")
  invisible(x)
}

#' @export
bq_plot.bq_segregate <- function(x, ...) {
  o <- x$observed
  long <- rbind(
    data.frame(trait = o$trait, phenotype = o$phenotype,
               kind = "observed", value = o$observed,
               stringsAsFactors = FALSE),
    data.frame(trait = o$trait, phenotype = o$phenotype,
               kind = "expected", value = o$expected,
               stringsAsFactors = FALSE))
  lab <- stats::setNames(paste0(x$best$trait, "\n", x$best$ratio,
                                "  (p = ", round(x$best$p, 3), ")"),
                         x$best$trait)
  long$facet <- lab[long$trait]
  ggplot2::ggplot(long, ggplot2::aes(x = .data$phenotype, y = .data$value,
                                     fill = .data$kind)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72),
                      width = 0.66) +
    ggplot2::facet_wrap(~ .data$facet, scales = "free") +
    ggplot2::scale_fill_manual(values = c(observed = "#6A2C91",
                                          expected = "#E9C46A"), name = NULL) +
    ggplot2::labs(title = "Observed and expected segregation",
      subtitle = "bars matching closely indicate the ratio fits the data",
      x = NULL, y = "number of plants", caption = "BKQualit") +
    theme_bq() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}
