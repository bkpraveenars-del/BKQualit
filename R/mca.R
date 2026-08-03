#' Multiple correspondence analysis and Gower clustering of accessions
#'
#' Places accessions and descriptor states in a common low-dimensional space by
#' multiple correspondence analysis (MCA), and groups the accessions by
#' hierarchical clustering of the Gower distance.
#'
#' MCA is computed from the singular value decomposition of the standardised
#' indicator matrix. Writing \eqn{Z} for the \eqn{n \times J} indicator of
#' \eqn{Q} descriptors, \eqn{P = Z/(nQ)}, and \eqn{r}, \eqn{c} for its row and
#' column margins, the decomposition is applied to
#' \eqn{D_r^{-1/2}(P - rc^{\top})D_c^{-1/2}}. Row and column principal
#' coordinates then satisfy the transition relation \eqn{F = (Z/Q)\,G_{std}},
#' so an accession sits at the centroid of the descriptor states it carries,
#' which is what makes an MCA biplot directly readable.
#'
#' The raw principal inertias of an indicator MCA are known to understate the
#' structure badly, because coding a descriptor into several columns injects
#' inertia that carries no information. The Benzecri adjustment is therefore
#' also reported: axes with inertia above \eqn{1/Q} are rescaled to
#' \eqn{\left(\frac{Q}{Q-1}\left(\lambda - \frac{1}{Q}\right)\right)^{2}}.
#' Interpret the adjusted percentages, not the raw ones.
#'
#' The Gower distance handles nominal and numeric descriptors together: nominal
#' descriptors contribute 0 when two accessions match and 1 when they differ,
#' numeric descriptors contribute the absolute difference divided by the range,
#' and the contributions are averaged. Clustering uses Ward's minimum-variance
#' criterion on that distance.
#'
#' @param data A data frame with one row per accession.
#' @param descriptors Character vector naming the categorical columns.
#' @param accession Character or \code{NULL}; the accession identifier column.
#' @param numeric_vars Character vector or \code{NULL}; numeric columns to add
#'   to the Gower distance (they do not enter the MCA).
#' @param k Number of clusters to cut the dendrogram into (default 3).
#' @param ndim Number of MCA dimensions to retain (default 2).
#' @return An object of class \code{bq_mca}.
#' @references
#' Gower JC (1971). A general coefficient of similarity and some of its
#' properties. \emph{Biometrics} \strong{27}, 857--871. \doi{10.2307/2528823}
#' @examples
#' d <- bq_data("dus")
#' res <- bq_mca(d,
#'   descriptors = c("growth_habit", "leaf_pubescence", "flower_colour",
#'                   "seed_colour", "seed_shape", "pod_curvature"),
#'   accession = "accession", numeric_vars = "plant_height_cm", k = 3)
#' res
#' bq_plot(res)
#' @export
bq_mca <- function(data, descriptors, accession = NULL, numeric_vars = NULL,
                   k = 3, ndim = 2) {
  stopifnot(all(descriptors %in% names(data)))
  keepcols <- c(descriptors, numeric_vars)
  d <- data[stats::complete.cases(data[keepcols]), , drop = FALSE]
  n <- nrow(d); Q <- length(descriptors)
  if (Q < 2) stop("At least two descriptors are needed for MCA.")
  ids <- if (is.null(accession)) paste0("A", seq_len(n))
         else as.character(d[[accession]])

  ## ---- indicator matrix ----------------------------------------------------
  levs <- lapply(descriptors, function(v) sort(unique(as.character(d[[v]]))))
  names(levs) <- descriptors
  colnm <- unlist(lapply(descriptors, function(v)
    paste0(v, "=", levs[[v]])), use.names = FALSE)
  Z <- matrix(0, n, length(colnm), dimnames = list(ids, colnm))
  for (v in descriptors) {
    idx <- match(paste0(v, "=", as.character(d[[v]])), colnm)
    Z[cbind(seq_len(n), idx)] <- 1
  }
  J <- ncol(Z)

  P <- Z / (n * Q)
  r <- rowSums(P); cm <- colSums(P)
  S <- diag(r^-0.5) %*% (P - outer(r, cm)) %*% diag(cm^-0.5)
  sv <- svd(S)
  keep <- sv$d > 1e-9
  dd <- sv$d[keep]; U <- sv$u[, keep, drop = FALSE]
  V <- sv$v[, keep, drop = FALSE]
  lam <- dd^2

  ## Benzecri adjustment
  thr <- 1 / Q
  sel <- lam > thr
  adj <- rep(NA_real_, length(lam))
  if (any(sel)) adj[sel] <- ((Q / (Q - 1)) * (lam[sel] - thr))^2

  ndim <- min(ndim, length(lam))
  F <- diag(r^-0.5) %*% U %*% diag(dd, nrow = length(dd))   # row principal
  G <- diag(cm^-0.5) %*% V %*% diag(dd, nrow = length(dd))  # col principal

  eig <- data.frame(dim = seq_along(lam), inertia = lam,
    pct_raw = 100 * lam / sum(lam),
    inertia_benzecri = adj,
    pct_benzecri = if (any(sel)) 100 * adj / sum(adj, na.rm = TRUE)
                   else rep(NA_real_, length(lam)),
    row.names = NULL)

  rows <- data.frame(accession = ids,
                     F[, seq_len(ndim), drop = FALSE],
                     row.names = NULL, stringsAsFactors = FALSE)
  names(rows)[-1] <- paste0("Dim", seq_len(ndim))
  cols <- data.frame(state = colnm,
                     descriptor = sub("=.*$", "", colnm),
                     G[, seq_len(ndim), drop = FALSE],
                     row.names = NULL, stringsAsFactors = FALSE)
  names(cols)[-(1:2)] <- paste0("Dim", seq_len(ndim))

  ## ---- Gower distance and clustering ---------------------------------------
  Dm <- matrix(0, n, n)
  used <- 0L
  for (v in descriptors) {
    x <- as.character(d[[v]])
    Dm <- Dm + outer(x, x, function(a, b) as.numeric(a != b))
    used <- used + 1L
  }
  for (v in numeric_vars) {
    x <- as.numeric(d[[v]]); rg <- diff(range(x))
    if (rg > 0) {
      Dm <- Dm + abs(outer(x, x, "-")) / rg
      used <- used + 1L
    }
  }
  Dm <- Dm / used
  dimnames(Dm) <- list(ids, ids)
  hc <- stats::hclust(stats::as.dist(Dm), method = "ward.D2")
  k <- max(1L, min(as.integer(k), n - 1L))
  grp <- stats::cutree(hc, k = k)
  rows$cluster <- factor(paste0("C", grp[ids]))

  ## cluster profile: modal state of every descriptor
  prof <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
    dg <- d[grp[ids] == g, , drop = FALSE]
    data.frame(cluster = paste0("C", g), size = nrow(dg),
      t(vapply(descriptors, function(v) {
        tb <- table(as.character(dg[[v]]))
        paste0(names(tb)[which.max(tb)], " (",
               round(100 * max(tb) / sum(tb)), "%)")
      }, character(1))),
      row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
  }))

  structure(list(eigen = eig, rows = rows, cols = cols, profile = prof,
                 distance = Dm, hclust = hc, k = k, Q = Q, J = J,
                 total_inertia = sum(lam), ndim = ndim),
            class = "bq_mca")
}

#' @export
print.bq_mca <- function(x, ...) {
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("  BKQualit | Multiple correspondence analysis and Gower clustering\n")
  cat(strrep("=", 74), "\n\n", sep = "")
  cat("Accessions ", nrow(x$rows), " | descriptors ", x$Q,
      " | states ", x$J, " | clusters ", x$k, "\n", sep = "")
  cat("Total inertia ", round(x$total_inertia, 4),
      "  (theoretical (J - Q)/Q = ", round((x$J - x$Q) / x$Q, 4), ")\n\n",
      sep = "")
  e <- utils::head(x$eigen, 6)
  print(data.frame(Dim = e$dim, Inertia = round(e$inertia, 4),
    `% raw` = round(e$pct_raw, 2),
    `% Benzecri` = ifelse(is.na(e$pct_benzecri), "-",
                          format(round(e$pct_benzecri, 2))),
    check.names = FALSE, stringsAsFactors = FALSE), row.names = FALSE)
  b <- x$eigen$pct_benzecri
  if (any(!is.na(b)))
    cat("\nFirst two axes carry ", round(sum(b[1:2], na.rm = TRUE), 1),
        "% of the adjusted inertia.\n", sep = "")
  cat("Use the Benzecri column: raw indicator inertias understate structure.\n")
  cat("\nCluster profiles (modal state of each descriptor):\n\n")
  print(x$profile, row.names = FALSE)
  cat("\nGower distance range: ", round(min(x$distance), 3), " to ",
      round(max(x$distance), 3),
      "  (0 = identical, 1 = differ in every descriptor)\n\n", sep = "")
  invisible(x)
}

#' @export
bq_plot.bq_mca <- function(x, type = c("biplot", "dendrogram", "scree"), ...) {
  type <- match.arg(type)
  if (type == "biplot") {
    ## stat_ellipse warns when a group has fewer than four points, so add it
    ## only when every cluster is large enough.
    big_enough <- min(table(x$rows$cluster)) >= 4L
    p <- ggplot2::ggplot(x$rows, ggplot2::aes(.data$Dim1, .data$Dim2)) +
      ggplot2::geom_hline(yintercept = 0, colour = "#D8D4E2") +
      ggplot2::geom_vline(xintercept = 0, colour = "#D8D4E2")
    if (big_enough)
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(colour = .data$cluster),
                                     type = "norm", linewidth = 0.5,
                                     alpha = 0.7)
    p +
      ggplot2::geom_point(ggplot2::aes(colour = .data$cluster), size = 2.6,
                          alpha = 0.85) +
      ggplot2::geom_point(data = x$cols, shape = 23, size = 2.8,
        fill = "#F2B705", colour = "#3A3450", inherit.aes = FALSE,
        ggplot2::aes(.data$Dim1, .data$Dim2)) +
      ggplot2::geom_text(data = x$cols, ggplot2::aes(.data$Dim1, .data$Dim2,
        label = sub("^.*=", "", .data$state)), inherit.aes = FALSE,
        size = 2.7, vjust = -0.9, colour = "#3A3450") +
      scale_colour_bq("contrast", name = "Cluster") +
      ggplot2::labs(title = "Multiple correspondence analysis biplot",
        subtitle = paste0("dots = accessions, diamonds = descriptor states; ",
          "Dim1 ", round(x$eigen$pct_benzecri[1], 1), "%, Dim2 ",
          round(x$eigen$pct_benzecri[2], 1), "% adjusted inertia"),
        x = "Dimension 1", y = "Dimension 2", caption = "BKQualit") +
      theme_bq()
  } else if (type == "scree") {
    e <- x$eigen
    e$pct <- ifelse(is.na(e$pct_benzecri), 0, e$pct_benzecri)
    ggplot2::ggplot(e, ggplot2::aes(factor(.data$dim), .data$pct)) +
      ggplot2::geom_col(fill = "#5B3E96", width = 0.65) +
      ggplot2::geom_line(ggplot2::aes(group = 1), colour = "#F2B705",
                         linewidth = 0.9) +
      ggplot2::geom_point(colour = "#F2B705", size = 2.4) +
      ggplot2::labs(title = "Scree plot of adjusted inertia",
        subtitle = "Benzecri-adjusted percentages",
        x = "dimension", y = "% of adjusted inertia", caption = "BKQualit") +
      theme_bq()
  } else {
    hc <- x$hclust
    ord <- hc$order
    labs <- hc$labels[ord]
    seg <- .dendro_segments(hc)
    grp <- stats::cutree(hc, k = x$k)
    lab <- data.frame(x = seq_along(ord), y = 0, label = labs,
      cluster = factor(paste0("C", grp[labs])), stringsAsFactors = FALSE)
    ggplot2::ggplot() +
      ggplot2::geom_segment(data = seg,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                     yend = .data$yend), colour = "#6C6880",
        linewidth = 0.45) +
      ggplot2::geom_text(data = lab,
        ggplot2::aes(.data$x, .data$y, label = .data$label,
                     colour = .data$cluster),
        angle = 90, hjust = 1.1, size = 2.5) +
      ggplot2::expand_limits(y = -max(hc$height) * 0.45) +
      scale_colour_bq("contrast", name = "Cluster") +
      ggplot2::labs(title = "Gower distance dendrogram",
        subtitle = paste0("Ward's minimum-variance linkage, cut into ",
                          x$k, " clusters"),
        x = NULL, y = "Gower distance", caption = "BKQualit") +
      theme_bq() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  }
}
