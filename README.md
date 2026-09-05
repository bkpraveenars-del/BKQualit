# BKQualit

**Analysis of Qualitative Traits, Segregation and Genetic Linkage**

[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-276DC3.svg)](https://cran.r-project.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22324861.svg)](https://doi.org/10.5281/zenodo.22324861)

Qualitative traits are usually analysed with a hand-typed chi-square and left
there. BKQualit does the rest of the job: it names the gene action, says
*why* a hypothesis failed, estimates linkage properly, and draws the figure.

Written for teaching and research in genetics and plant breeding by
**Dr. Praveen Kumar B. K.**, Department of Genetics and Plant Breeding,
Agriculture University, Jodhpur.

## Citation

Kumar, P. (2026). *BKQualit: Analysis of Qualitative Traits, Segregation and
Genetic Linkage*. Zenodo. https://doi.org/10.5281/zenodo.22324861

---

## Install

```r
install.packages("BKQualit")            # once on CRAN
# development version
remotes::install_github("bkpraveenars-del/BKQualit")
```

## Five functions, one plotting verb

| Function | Question it answers |
|---|---|
| `bq_segregate()` | Which Mendelian ratio do these counts fit, and what gene action does it imply? |
| `bq_heterogeneity()` | Is the ratio wrong, or do the families simply differ from one another? |
| `bq_linkage()` | How far apart are these two loci, and is the linkage real? |
| `bq_diversity()` | Which of my descriptors actually discriminate between accessions? |
| `bq_mca()` | How do accessions group, and which descriptor states drive the grouping? |

Every result prints a plain-English verdict and plots with `bq_plot()`.

---

## Sixty seconds with the package

```r
library(BKQualit)

## 1. Which ratio fits? Do not guess - fit them all.
res <- bq_segregate(bq_data("f2"), phenotype = "phenotype",
                    count = "count", trait = "trait")
res
bq_plot(res)
```

It fits every classical ratio with the right number of classes, ranks them by
goodness of fit, and tells you what each one means (output abridged; the full
printout also carries the Monte Carlo p value, the observed-versus-expected
table, and every candidate ratio that was tried):

```
BEST-FITTING RATIO for each trait:

        Trait   n  Ratio  ChiSq df      p Verdict
 flower_colour 478    3:1 0.0251  1 0.8741    fits
     seed_coat 405    9:7 0.0004  1 0.9850    fits
    grain_type 462  9:6:1 0.0067  2 0.9966    fits
      awn_type 592   13:3 0.7096  1 0.3996    fits

Implied gene action:
   seed_coat : 9:7 - duplicate recessive epistasis (complementary genes)
   grain_type : 9:6:1 - duplicate genes with cumulative effect
   awn_type : 13:3 - dominant suppression epistasis
```

```r
## 2. The ratio failed. Is the hypothesis wrong, or are the families unalike?
het <- bq_heterogeneity(bq_data("families"), family = "family",
                        classes = c("dominant", "recessive"),
                        trait = "trait", ratio = c(3, 1))
het
```

```
Trait: disease_reaction
        Source   ChiSq df      p Sig
        Pooled 25.5670  1 0.0000 ***
 Heterogeneity 58.8975 11 0.0000 ***
         Total 84.4645 12 0.0000 ***
  Interpretation: families are HETEROGENEOUS (p = 0.0000).
  Do not pool them; examine the families separately.
```

That distinction is the whole point. A single pooled chi-square would have
called the 3:1 hypothesis wrong; the partition shows the real problem is that
the families disagree, so pooling them was never valid.

```r
## 3. Linkage by maximum likelihood, not by the old product-ratio table
lk <- bq_linkage(bq_data("linkage"), pair = "pair",
                 classes = c("AB", "Ab", "aB", "ab"),
                 design = "design", phase = "phase")
lk
bq_plot(lk)                    # likelihood profiles with the LOD = 3 line
bq_plot(lk, type = "estimates")

## 4. Which descriptors earn their place on the DUS form?
dv <- bq_diversity(bq_data("dus"),
        descriptors = c("growth_habit", "leaf_pubescence", "flower_colour",
                        "seed_colour", "seed_shape"))
dv
bq_plot(dv)

## 5. Group the accessions, and see what drives the grouping
mc <- bq_mca(bq_data("dus"),
        descriptors = c("growth_habit", "leaf_pubescence", "flower_colour",
                        "seed_colour", "seed_shape", "pod_curvature"),
        accession = "accession", numeric_vars = "plant_height_cm", k = 3)
mc
bq_plot(mc)                     # accessions and descriptor states together
bq_plot(mc, type = "dendrogram")
```

---

## What is done differently here

**The chi-square is partitioned, not just reported.** Pooled plus heterogeneity
equals total, exactly. A poor fit is attributed either to the wrong hypothesis
or to variation between families. Most software gives you one number and leaves
the diagnosis to you.

**Yates' correction is applied only at one degree of freedom.** That is where
it is defined. Applying it to a 9:3:3:1 table, which is common in practice,
makes the test conservative for no reason.

**Sparse tables get a Monte Carlo exact test.** The chi-square approximation is
unreliable when an expected count falls below five, which happens routinely in
the rarest class of a 15:1 or 63:1 table.

**Linkage is estimated by maximising the likelihood.** The recombination
fraction comes from the multinomial likelihood, not the product-ratio
approximation, and the confidence interval comes from likelihood-ratio
inversion rather than a symmetric standard error - which matters, because the
likelihood is skewed as r approaches its boundary at 0.5. Coupling and
repulsion phase, F2 and backcross, are all handled.

**MCA inertias are Benzecri-adjusted.** Raw indicator inertias badly understate
the structure because splitting a descriptor into several columns injects
inertia that carries no information. On the bundled data the first two axes
carry 30% of the raw inertia but 86% of the adjusted inertia. The adjusted
figure is the honest one.

**Gower distance handles mixed data.** Nominal and numeric descriptors enter
the same distance, so plant height can sit alongside flower colour.

---

## Companion packages

| Package | Scope |
|---|---|
| [**BKBreed**](https://github.com/bkpraveenars-del/BKBreed) | Experimental designs, variability, correlation and path analysis, D2 diversity, stability and multi-location trials, combining ability |
| [**BKMutate**](https://github.com/bkpraveenars-del/BKMutate) | Induced mutagenesis: dose-response and LD50, mutagenic effectiveness and efficiency, chlorophyll mutation spectrum |
| **BKQualit** | Qualitative traits: segregation, heterogeneity, linkage, descriptor diversity |

---

## A note on the bundled data

The four example datasets are **simulated from known parameters**, so the
correct answer is known in advance and the package can be tested against
arithmetic rather than against itself. The generating parameters, and the
values the package recovers from them, are tabulated in
[`data-raw/PROVENANCE.md`](data-raw/PROVENANCE.md). They are not field
observations and should not be cited as such.

---

## Citation

```r
citation("BKQualit")
```

> B. K., Praveen Kumar (2026). *BKQualit: Analysis of Qualitative Traits,
> Segregation and Genetic Linkage*. R package version 0.1.0.
> https://github.com/bkpraveenars-del/BKQualit

---

## References

Mather K (1951). *The Measurement of Linkage in Heredity*, 2nd ed. Methuen, London.

Allard RW (1956). Formulas and tables to facilitate the calculation of
recombination values in heredity. *Hilgardia* **24**, 235-278.
[doi:10.3733/hilg.v24n10p235](https://doi.org/10.3733/hilg.v24n10p235)

Gower JC (1971). A general coefficient of similarity and some of its
properties. *Biometrics* **27**, 857-871.
[doi:10.2307/2528823](https://doi.org/10.2307/2528823)

---

## License

GPL-3
