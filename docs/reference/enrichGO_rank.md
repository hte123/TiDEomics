# GO enrichment with ranked gene list

Gene ontology enrichment analysis of a ranked gene list with
clusterProfiler. Genes can be ranked by variance decomposition results.

## Usage

``` r
enrichGO_rank(
  rank_table,
  gene_rank_by,
  OrgDb,
  keyType = "SYMBOL",
  go_rank_by = "p.adjust",
  category = NULL,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  simplify = FALSE,
  simplify_cutoff = 0.7,
  simplify_by = "p.adjust",
  simplify_select_fun = min,
  simplify_measure = "Wang",
  ...
)
```

## Arguments

- rank_table:

  A data frame with at least two columns: 'Feature' for gene names and
  one or more columns of variables for ranking the genes, e.g. output of
  [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md)
  containing variance decomposition results

- gene_rank_by:

  Variable in `rank_table` to rank the genes by, e.g. "Time", "Group" in
  the output of
  [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md)

- OrgDb:

  Organism database, e.g. org.Hs.eg.db, org.Mm.eg.db

- keyType:

  Available options are `AnnotationDbi::keytypes(OrgDb)` (default is
  "SYMBOL")

- go_rank_by:

  Variable in the GO enrichment result to rank the GO terms by (default
  is "p.adjust", other options include "pvalue", "qvalue", "NES",
  "setSize", "enrichmentScore", etc.)

- category:

  GO category to analyze (default is all three of BP, MF, CC)

- pvalueCutoff:

  Parameter of
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)
  (default is 0.05)

- pAdjustMethod:

  Parameter of
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)
  (default is "BH")

- simplify:

  Whether to simplify the GO terms by removing redundant terms with
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html)
  (default is FALSE)

- simplify_cutoff:

  Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  cutoff for similarity when simplifying GO terms (default is 0.7)

- simplify_by:

  Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  method to choose the representative term (default is "p.adjust")

- simplify_select_fun:

  Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  function to select the representative term (default is `min`)

- simplify_measure:

  Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  method to calculate similarity when simplifying GO terms (default is
  "Wang")

- ...:

  additional arguments passed to
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)

## Value

A list containing `all`, a list of `gseaResult` objects with the GSEA
results. If `simplify = TRUE`, the output also contains `simplified`, a
list of simplified `gseaResult` objects.

## Examples

``` r
if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(org.Mm.eg.db)
    data(example_obj)

    var_decomp <- decomp_variance(example_obj,
        features = rownames(example_obj)[1:100], assay = 1)
    example_go_rank <- enrichGO_rank(var_decomp, gene_rank_by = "Time",
        OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
        category = "BP")
}
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
#> 100 features evaluated: 71 successful; 0 errors; 29 warnings; 0 skipped before fitting for insufficient observations or variance.
#> Warning: There were 1 pathways for which P-values were not calculated properly; pvalue, NES and log2err are set to NA for them. The usual cause is an unbalanced gene-level statistic (no genes on one side of the ranking), which makes the permutation null degenerate; a gene set covering the whole ranked list is excluded before this point. You can try to increase nPermSimple, or pass reportNA = TRUE to keep these rows.
#> Warning: Invalid p-values detected (NA, non-finite, <0, or >1). qvalue will be computed on valid p-values only.
#> Warning: 1 of 151 gene set(s) have no estimable p-value and were dropped (pvalue, NES and log2err are NA for them). Pass reportNA = TRUE to keep them, or raise nPermSimple.
#> No significant GO terms found for category: BP
#> No significant GO terms found in any category.
```
