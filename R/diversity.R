#' Diversity of qualitative descriptors
#'
#' Summarises a set of categorical descriptors (the kind recorded on a
#' distinctness, uniformity and stability, or DUS, form) into per-descriptor
#' diversity indices and identifies which descriptors actually discriminate
#' among accessions.
#'
#' Three complementary indices are reported for each descriptor. The
#' Shannon-Weaver index \eqn{H' = -\sum p_i \log p_i} is sensitive to rare
#' states; the Simpson index \eqn{1 - \sum p_i^2} is dominated by common states;
#' and evenness \eqn{J' = H'/\log s} rescales Shannon by the number of observed
#' states \eqn{s} so that descriptors with different numbers of states can be
#' compared directly. A descriptor with a single observed state is
#' monomorphic, contributes nothing to discrimination, and is flagged.
#'
#' A descriptor is also useful only if it separates accessions. The
#' \emph{discriminating power} reported here is the probability that two
#' randomly chosen accessions differ in that descriptor, which equals the
#' Simpson index; descriptors are ranked on it so that the least informative
#' can be dropped from a characterisation form.
#'
#' @param data A data frame with one row per accession.
#' @param descriptors Character vector naming the categorical columns.
#' @param accession Character or \code{NULL}; the accession identifier column.
#' @param group Character or \code{NULL}; an optional grouping column, giving
#'   diversity within each group as well as overall.
#' @return An object of class \code{bq_diversity}.
#' @examples
#' d <- bq_data("dus")
#' res <- bq_diversity(d,
#'   descriptors = c("growth_habit", "leaf_pubescence", "flower_colour",
#'                   "seed_colour", "seed_shape"),
#'   accession = "accession")
#' res
#' bq_plot(res)
#' @export
bq_diversity <- function(data, descriptors, accession = NULL, group = NULL) {
  stopifnot(all(descriptors %in% names(data)))
  d <- data[stats::complete.cases(data[descriptors]), , drop = FALSE]

  one <- function(v) {
    tb <- table(as.character(v))
    p <- as.numeric(tb) / sum(tb)
    s <- length(p)
    H <- -sum(p * log(p))
    D <- 1 - sum(p^2)
    J <- if (s > 1) H / log(s) else 0
    c(states = s, shannon = H, simpson = D, evenness = J,
      max_freq = max(p), n = sum(tb))
  }

  ov <- as.data.frame(t(vapply(descriptors, function(v) one(d[[v]]),
                               numeric(6))))
  ov$descriptor <- rownames(ov); rownames(ov) <- NULL
  ov <- ov[c("descriptor", "n", "states", "shannon", "simpson", "evenness",
             "max_freq")]
  ov$discriminating_power <- ov$simpson
  ov$status <- ifelse(ov$states < 2, "monomorphic",
               ifelse(ov$evenness >= 0.75, "highly informative",
               ifelse(ov$evenness >= 0.40, "informative", "weak")))
  ov <- ov[order(-ov$discriminating_power), ]
  rownames(ov) <- NULL

  ## frequency table of every state
  freq <- do.call(rbind, lapply(descriptors, function(v) {
    tb <- table(as.character(d[[v]]))
    data.frame(descriptor = v, state = names(tb), count = as.numeric(tb),
               frequency = as.numeric(tb) / sum(tb),
               row.names = NULL, stringsAsFactors = FALSE)
  }))

  by_group <- NULL
  if (!is.null(group)) {
    by_group <- do.call(rbind, lapply(unique(as.character(d[[group]])),
      function(g) {
        dg <- d[as.character(d[[group]]) == g, , drop = FALSE]
        m <- as.data.frame(t(vapply(descriptors,
               function(v) one(dg[[v]]), numeric(6))))
        m$descriptor <- rownames(m); m$group <- g
        rownames(m) <- NULL
        m[c("group", "descriptor", "n", "states", "shannon", "simpson",
            "evenness")]
      }))
  }

  ## overall mean diversity per accession is not defined; report totals
  structure(list(indices = ov, frequencies = freq, by_group = by_group,
                 n_accessions = nrow(d),
                 mean_shannon = mean(ov$shannon),
                 mean_evenness = mean(ov$evenness),
                 accession = accession),
            class = "bq_diversity")
}

#' @export
print.bq_diversity <- function(x, ...) {
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("  BKQualit | Diversity of qualitative descriptors\n")
  cat(strrep("=", 74), "\n\n", sep = "")
  cat("Accessions scored : ", x$n_accessions, "\n", sep = "")
  cat("Descriptors       : ", nrow(x$indices), "\n", sep = "")
  cat("Mean Shannon H'   : ", round(x$mean_shannon, 4), "\n", sep = "")
  cat("Mean evenness J'  : ", round(x$mean_evenness, 4), "\n\n", sep = "")
  i <- x$indices
  print(data.frame(Descriptor = i$descriptor, States = i$states,
    `Shannon H` = round(i$shannon, 4), `Simpson D` = round(i$simpson, 4),
    `Evenness J` = round(i$evenness, 4),
    `Discrim. power` = round(i$discriminating_power, 4),
    Status = i$status, check.names = FALSE, stringsAsFactors = FALSE),
    row.names = FALSE)
  mono <- i$descriptor[i$states < 2]
  cat("\n")
  if (length(mono))
    cat("Monomorphic (no discriminating value, consider dropping): ",
        paste(mono, collapse = ", "), "\n", sep = "")
  best <- i$descriptor[1]
  cat("Most discriminating descriptor: ", best, " (", round(i$simpson[1], 3),
      " chance that two random accessions differ)\n", sep = "")
  cat("\nH' rewards rare states; D and discriminating power reward balance;\n")
  cat("J' rescales H' by the number of states so descriptors are comparable.\n\n")
  invisible(x)
}

#' @export
bq_plot.bq_diversity <- function(x, type = c("indices", "frequencies"), ...) {
  type <- match.arg(type)
  if (type == "indices") {
    i <- x$indices
    long <- rbind(
      data.frame(descriptor = i$descriptor, index = "Shannon H'",
                 value = i$shannon, stringsAsFactors = FALSE),
      data.frame(descriptor = i$descriptor, index = "Simpson D",
                 value = i$simpson, stringsAsFactors = FALSE),
      data.frame(descriptor = i$descriptor, index = "Evenness J'",
                 value = i$evenness, stringsAsFactors = FALSE))
    long$descriptor <- factor(long$descriptor, levels = rev(i$descriptor))
    ggplot2::ggplot(long, ggplot2::aes(.data$descriptor, .data$value,
                                       fill = .data$index)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78),
                        width = 0.72) +
      ggplot2::coord_flip() +
      scale_fill_bq("descriptor", name = NULL) +
      ggplot2::labs(title = "Diversity of qualitative descriptors",
        subtitle = "descriptors ordered by discriminating power",
        x = NULL, y = "index value", caption = "BKQualit") +
      theme_bq()
  } else {
    f <- x$frequencies
    ggplot2::ggplot(f, ggplot2::aes(x = .data$descriptor, y = .data$frequency,
                                    fill = .data$state)) +
      ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = .data$state),
        position = ggplot2::position_stack(vjust = 0.5), size = 2.7,
        colour = "white") +
      ggplot2::coord_flip() +
      scale_fill_bq("phenotype", name = NULL) +
      ggplot2::labs(title = "State frequencies within each descriptor",
        subtitle = "wide, even bands indicate an informative descriptor",
        x = NULL, y = "proportion of accessions", caption = "BKQualit") +
      theme_bq() +
      ggplot2::theme(legend.position = "none")
  }
}
