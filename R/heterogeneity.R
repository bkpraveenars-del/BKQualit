#' Pooled and heterogeneity chi-square across families
#'
#' Partitions the total chi-square computed over several segregating families
#' into a pooled component and a heterogeneity component. This is the standard
#' way to decide \emph{why} a segregation hypothesis fails.
#'
#' Writing \eqn{\chi^2_i} for the goodness-of-fit statistic of family \eqn{i},
#' the total is \eqn{\sum_i \chi^2_i} on \eqn{k(c-1)} degrees of freedom, the
#' pooled statistic is computed from the summed counts on \eqn{c-1} degrees of
#' freedom, and the heterogeneity statistic is the difference on
#' \eqn{(k-1)(c-1)} degrees of freedom, where \eqn{k} is the number of families
#' and \eqn{c} the number of phenotype classes. The partition is exact.
#'
#' The two components answer different questions. A significant \emph{pooled}
#' chi-square with a non-significant heterogeneity means the families behave
#' alike but the hypothesised ratio is wrong. A significant \emph{heterogeneity}
#' means the families do not behave alike, so pooling them is not justified and
#' the pooled test should not be interpreted.
#'
#' @param data A data frame with one row per family.
#' @param family Character; the family identifier column.
#' @param classes Character vector; the columns holding counts of each
#'   phenotype class.
#' @param trait Character or \code{NULL}; optional column identifying several
#'   traits, each partitioned separately.
#' @param ratio Numeric vector giving the expected ratio, default
#'   \code{c(3, 1)}.
#' @return An object of class \code{bq_heterogeneity}.
#' @examples
#' d <- bq_data("families")
#' res <- bq_heterogeneity(d, family = "family",
#'                         classes = c("dominant", "recessive"),
#'                         trait = "trait", ratio = c(3, 1))
#' res
#' bq_plot(res)
#' @export
bq_heterogeneity <- function(data, family, classes, trait = NULL,
                             ratio = c(3, 1)) {
  stopifnot(all(c(family, classes, trait) %in% names(data)))
  if (length(ratio) != length(classes))
    stop("ratio has ", length(ratio), " parts but ", length(classes),
         " class columns were given.")
  d <- data[stats::complete.cases(data[c(family, classes, trait)]), ]
  p <- as.numeric(ratio) / sum(as.numeric(ratio))
  cc <- length(classes)

  key <- if (is.null(trait)) rep("all", nrow(d)) else as.character(d[[trait]])
  per <- list(); summ <- list()

  for (g in unique(key)) {
    dg <- d[key == g, , drop = FALSE]
    M <- as.matrix(dg[classes])
    n_i <- rowSums(M)
    E <- outer(n_i, p)
    chi_i <- rowSums((M - E)^2 / E)
    df_i <- cc - 1

    total <- sum(chi_i); df_total <- nrow(M) * df_i
    pooled_obs <- colSums(M)
    e_pool <- sum(pooled_obs) * p
    chi_pool <- sum((pooled_obs - e_pool)^2 / e_pool); df_pool <- df_i
    chi_het <- total - chi_pool; df_het <- (nrow(M) - 1) * df_i

    per[[g]] <- data.frame(trait = g, family = as.character(dg[[family]]),
      n = n_i, chisq = chi_i, df = df_i,
      p = stats::pchisq(chi_i, df_i, lower.tail = FALSE),
      row.names = NULL, stringsAsFactors = FALSE)
    per[[g]]$Sig <- .stars(per[[g]]$p)

    summ[[g]] <- data.frame(
      trait = g,
      source = c("Pooled", "Heterogeneity", "Total"),
      chisq = c(chi_pool, chi_het, total),
      df = c(df_pool, df_het, df_total),
      p = c(stats::pchisq(chi_pool, df_pool, lower.tail = FALSE),
            stats::pchisq(chi_het, df_het, lower.tail = FALSE),
            stats::pchisq(total, df_total, lower.tail = FALSE)),
      row.names = NULL, stringsAsFactors = FALSE)
    summ[[g]]$Sig <- .stars(summ[[g]]$p)
  }

  pt <- do.call(rbind, per); rownames(pt) <- NULL
  st <- do.call(rbind, summ); rownames(st) <- NULL

  structure(list(per_family = pt, partition = st, ratio = ratio,
                 classes = classes),
            class = "bq_heterogeneity")
}

#' @export
print.bq_heterogeneity <- function(x, ...) {
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("  BKQualit | Pooled and heterogeneity chi-square  (ratio ",
      paste(x$ratio, collapse = ":"), ")\n", sep = "")
  cat(strrep("=", 74), "\n\n", sep = "")
  s <- x$partition
  for (g in unique(s$trait)) {
    sg <- s[s$trait == g, ]
    cat("Trait: ", g, "\n", sep = "")
    print(data.frame(Source = sg$source, ChiSq = round(sg$chisq, 4),
                     df = sg$df, p = round(sg$p, 4), Sig = sg$Sig,
                     stringsAsFactors = FALSE), row.names = FALSE)
    pooled <- sg[sg$source == "Pooled", ]
    het <- sg[sg$source == "Heterogeneity", ]
    cat("  Interpretation: ")
    if (het$p < 0.05) {
      cat("families are HETEROGENEOUS (p = ", round(het$p, 4), ").\n",
          "  Do not pool them; examine the families separately.\n", sep = "")
    } else if (pooled$p < 0.05) {
      cat("families agree, but the ratio ", paste(x$ratio, collapse = ":"),
          " is REJECTED (p = ", round(pooled$p, 4), ").\n",
          "  Consider a different genetic hypothesis.\n", sep = "")
    } else {
      cat("families agree AND the ratio fits.\n",
          "  The hypothesis is supported across all families.\n", sep = "")
    }
    cat("\n")
  }
  cat("Per-family goodness of fit:\n\n")
  pf <- x$per_family
  pf$chisq <- round(pf$chisq, 4); pf$p <- round(pf$p, 4)
  print(pf, row.names = FALSE)
  cat("\nPartition: Pooled + Heterogeneity = Total (exact, by construction).\n\n")
  invisible(x)
}

#' @export
bq_plot.bq_heterogeneity <- function(x, ...) {
  pf <- x$per_family
  crit <- stats::qchisq(0.95, pf$df[1])
  ggplot2::ggplot(pf, ggplot2::aes(x = .data$family, y = .data$chisq,
                                   fill = .data$p > 0.05)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_hline(yintercept = crit, linetype = "dashed",
                        colour = "#B3261E", linewidth = 0.5) +
    ggplot2::annotate("text", x = 1, y = crit, vjust = -0.6, hjust = 0,
                      size = 3, colour = "#B3261E",
                      label = paste0("5% critical value = ", round(crit, 2))) +
    ggplot2::facet_wrap(~ .data$trait, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#1B7340",
                                          `FALSE` = "#B3261E"),
      labels = c(`TRUE` = "fits ratio", `FALSE` = "deviates"), name = NULL) +
    ggplot2::labs(title = "Goodness of fit family by family",
      subtitle = paste0("expected ratio ", paste(x$ratio, collapse = ":"),
                        "; bars above the line deviate significantly"),
      x = NULL, y = expression(chi^2), caption = "BKQualit") +
    theme_bq() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
