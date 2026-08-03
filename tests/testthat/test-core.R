## ---------------------------------------------------------------------------
## Every expected value below was computed independently before being written
## here, so these tests check the package against arithmetic, not against
## itself.
## ---------------------------------------------------------------------------

test_that("bq_segregate reproduces hand-computed chi-square values", {
  ## 312 : 104 against 3:1 gives expected 312 and 104 exactly, so chi2 = 0
  res <- bq_segregate(counts = c(purple = 312, white = 104), ratio = c(3, 1),
                      simulate = FALSE)
  expect_equal(res$fits$chisq[1], 0, tolerance = 1e-10)
  expect_equal(res$fits$df[1], 1)
  expect_equal(res$fits$ratio[1], "3:1")

  ## 290 : 110 against 3:1. E = 300, 100.
  ## chi2 = 100/300 + 100/100 = 4/3
  r2 <- bq_segregate(counts = c(A = 290, B = 110), ratio = c(3, 1),
                     simulate = FALSE)
  expect_equal(r2$fits$chisq[1], 4 / 3, tolerance = 1e-9)

  ## 315:108:101:32 against 9:3:3:1 gives chi2 = 0.47002 on 3 df
  r3 <- bq_segregate(counts = c(AB = 315, Ab = 108, aB = 101, ab = 32),
                     ratio = c(9, 3, 3, 1), simulate = FALSE)
  expect_equal(r3$fits$chisq[1], 0.47002, tolerance = 1e-4)
  expect_equal(r3$fits$df[1], 3)
})

test_that("Yates correction is applied only at one degree of freedom", {
  ## At 1 df the correction exists and must shrink the statistic.
  ## (|290-300|-0.5)^2/300 + (|110-100|-0.5)^2/100 = 1.203333
  one_df <- bq_segregate(counts = c(A = 290, B = 110), ratio = c(3, 1),
                         simulate = FALSE)$fits
  expect_false(is.na(one_df$chisq_yates[1]))
  expect_equal(one_df$chisq_yates[1], 1.2033333, tolerance = 1e-6)
  expect_lt(one_df$chisq_yates[1], one_df$chisq[1])

  ## Above 1 df the correction is undefined and must not be reported.
  three_df <- bq_segregate(counts = c(AB = 315, Ab = 108, aB = 101, ab = 32),
                           ratio = c(9, 3, 3, 1), simulate = FALSE)$fits
  expect_true(is.na(three_df$chisq_yates[1]))
  expect_true(is.na(three_df$p_yates[1]))
})

test_that("candidate ratios are ranked and the gene action is named", {
  d <- bq_data("f2")
  res <- bq_segregate(d, phenotype = "phenotype", count = "count",
                      trait = "trait", simulate = FALSE)
  ## every trait gets exactly one best row
  expect_equal(nrow(res$best), length(unique(d$trait)))
  ## the fits table is sorted by p within each trait
  for (g in unique(res$fits$trait)) {
    p <- res$fits$p[res$fits$trait == g]
    expect_false(is.unsorted(rev(p)))
  }
  expect_true(all(nzchar(res$best$gene_action)))
  ## the simulated grain_type data were generated on a 9:6:1 expectation
  gt <- res$best[res$best$trait == "grain_type", ]
  expect_equal(gt$ratio, "9:6:1")
})

test_that("chi-square partition is exact: pooled + heterogeneity = total", {
  d <- bq_data("families")
  res <- bq_heterogeneity(d, family = "family",
                          classes = c("dominant", "recessive"),
                          trait = "trait", ratio = c(3, 1))
  s <- res$partition
  for (g in unique(s$trait)) {
    sg <- s[s$trait == g, ]
    pooled <- sg$chisq[sg$source == "Pooled"]
    het    <- sg$chisq[sg$source == "Heterogeneity"]
    total  <- sg$chisq[sg$source == "Total"]
    expect_equal(pooled + het, total, tolerance = 1e-10)
    expect_equal(sg$df[sg$source == "Pooled"] +
                   sg$df[sg$source == "Heterogeneity"],
                 sg$df[sg$source == "Total"])
  }
  ## the total must also equal the sum of the per-family statistics
  for (g in unique(s$trait)) {
    tot <- s$chisq[s$trait == g & s$source == "Total"]
    per <- sum(res$per_family$chisq[res$per_family$trait == g])
    expect_equal(tot, per, tolerance = 1e-10)
  }
})

test_that("linkage recovers a known recombination fraction", {
  ## Counts generated from r = 0.18 in coupling phase, n = 200000
  d <- data.frame(pair = "X-Y", design = "F2", phase = "coupling",
                  AB = 133620, Ab = 16380, aB = 16380, ab = 33620,
                  stringsAsFactors = FALSE)
  res <- bq_linkage(d, pair = "pair", classes = c("AB", "Ab", "aB", "ab"),
                    design = "design", phase = "phase")
  e <- res$estimates
  expect_equal(e$r[1], 0.18, tolerance = 1e-3)
  expect_lt(e$lower[1], 0.18)
  expect_gt(e$upper[1], 0.18)
  expect_gt(e$LOD[1], 3)
  expect_true(e$linked[1])
  expect_equal(e$map_distance_cM[1], 100 * e$r[1], tolerance = 1e-10)
})

test_that("independent assortment gives r near 0.5 and a LOD near zero", {
  d <- data.frame(pair = "I-J", design = "F2", phase = "coupling",
                  AB = 900, Ab = 300, aB = 300, ab = 100,
                  stringsAsFactors = FALSE)
  res <- bq_linkage(d, pair = "pair", classes = c("AB", "Ab", "aB", "ab"),
                    design = "design", phase = "phase")
  e <- res$estimates
  expect_gt(e$r[1], 0.45)
  expect_lt(e$LOD[1], 1)
  expect_false(e$linked[1])
  expect_gt(e$p_independence[1], 0.05)
})

test_that("diversity indices behave correctly at the extremes", {
  d <- data.frame(a = rep("one", 20),                # monomorphic
                  b = rep(c("x", "y"), each = 10),   # perfectly even
                  stringsAsFactors = FALSE)
  i <- bq_diversity(d, descriptors = c("a", "b"))$indices
  mono <- i[i$descriptor == "a", ]
  even <- i[i$descriptor == "b", ]
  expect_equal(mono$shannon, 0, tolerance = 1e-12)
  expect_equal(mono$simpson, 0, tolerance = 1e-12)
  expect_equal(mono$evenness, 0, tolerance = 1e-12)
  expect_equal(mono$status, "monomorphic")
  expect_equal(even$shannon, log(2), tolerance = 1e-10)
  expect_equal(even$simpson, 0.5, tolerance = 1e-12)
  expect_equal(even$evenness, 1, tolerance = 1e-10)
  ## discriminating power is the Simpson index by construction
  expect_equal(i$discriminating_power, i$simpson, tolerance = 1e-12)
})

test_that("MCA reproduces its theoretical total inertia and is centred", {
  d <- bq_data("dus")
  desc <- c("growth_habit", "leaf_pubescence", "flower_colour",
            "seed_colour", "seed_shape", "pod_curvature")
  res <- bq_mca(d, descriptors = desc, accession = "accession",
                numeric_vars = "plant_height_cm", k = 3)
  ## total inertia of an indicator MCA is exactly (J - Q) / Q
  expect_equal(res$total_inertia, (res$J - res$Q) / res$Q, tolerance = 1e-9)
  ## rows carry equal mass, so their coordinates must average to zero
  expect_equal(mean(res$rows$Dim1), 0, tolerance = 1e-8)
  expect_equal(mean(res$rows$Dim2), 0, tolerance = 1e-8)
  ## the Benzecri adjustment concentrates inertia on the leading axes
  expect_gte(res$eigen$pct_benzecri[1], res$eigen$pct_raw[1])
  expect_equal(nrow(res$rows), nrow(d))
  expect_equal(nlevels(res$rows$cluster), 3L)
})

test_that("Gower distance is a proper distance on [0, 1]", {
  d <- bq_data("dus")
  res <- bq_mca(d, descriptors = c("growth_habit", "flower_colour",
                                   "seed_shape"),
                accession = "accession", k = 2)
  D <- res$distance
  expect_equal(D, t(D), tolerance = 1e-12)
  expect_true(all(diag(D) == 0))
  expect_true(all(D >= 0 & D <= 1))
})

test_that("every result prints and plots without error", {
  f2 <- bq_segregate(bq_data("f2"), phenotype = "phenotype",
                     count = "count", trait = "trait", simulate = FALSE)
  he <- bq_heterogeneity(bq_data("families"), family = "family",
                         classes = c("dominant", "recessive"),
                         trait = "trait")
  li <- bq_linkage(bq_data("linkage"), pair = "pair",
                   classes = c("AB", "Ab", "aB", "ab"),
                   design = "design", phase = "phase")
  dv <- bq_diversity(bq_data("dus"),
                     descriptors = c("growth_habit", "flower_colour",
                                     "seed_shape"))
  mc <- bq_mca(bq_data("dus"),
               descriptors = c("growth_habit", "flower_colour",
                               "seed_shape", "seed_colour"),
               accession = "accession", k = 3)
  for (obj in list(f2, he, li, dv, mc)) {
    expect_output(print(obj))
    expect_s3_class(bq_plot(obj), "ggplot")
  }
  expect_s3_class(bq_plot(li, type = "estimates"), "ggplot")
  expect_s3_class(bq_plot(dv, type = "frequencies"), "ggplot")
  expect_s3_class(bq_plot(mc, type = "dendrogram"), "ggplot")
  expect_s3_class(bq_plot(mc, type = "scree"), "ggplot")
  expect_error(bq_plot(structure(list(), class = "not_a_bq_object")))
})

test_that("bq_data loads all four datasets", {
  for (nm in c("f2", "families", "linkage", "dus")) {
    d <- bq_data(nm)
    expect_s3_class(d, "data.frame")
    expect_gt(nrow(d), 0)
  }
})

test_that("bq_save writes files only where it is told to", {
  res <- bq_segregate(bq_data("f2"), phenotype = "phenotype",
                      count = "count", trait = "trait", simulate = FALSE)
  dir <- file.path(tempdir(), "bq_test_figs")
  files <- suppressMessages(bq_save(res, path = dir))
  expect_true(all(file.exists(files)))
  unlink(dir, recursive = TRUE)
})

test_that("input errors are caught with a clear message", {
  expect_error(bq_segregate(), "Supply either")
  expect_error(bq_segregate(counts = c(a = 10, b = 5, c = 5),
                            ratio = c(3, 1)), "phenotype classes")
  expect_error(bq_heterogeneity(bq_data("families"), family = "family",
                                classes = c("dominant", "recessive"),
                                ratio = c(9, 3, 3, 1)), "class columns")
  expect_error(bq_linkage(bq_data("linkage"), pair = "pair",
                          classes = c("AB", "Ab")), "exactly four")
})
