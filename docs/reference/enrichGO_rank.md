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

  (Optional) Available options are `AnnotationDbi::keytypes(OrgDb)`
  (default is "SYMBOL")

- go_rank_by:

  (Optional) Variable in the GO enrichment result to rank the GO terms
  by (default is "p.adjust", other options include "pvalue", "qvalue",
  "NES", "setSize", "enrichmentScore", etc.)

- category:

  (Optional) GO category to analyze (default is all three of BP, MF, CC)

- pvalueCutoff:

  (Optional) Parameter of
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)
  (default is 0.05)

- pAdjustMethod:

  (Optional) Parameter of
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)
  (default is "BH")

- ...:

  additional arguments passed to
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)

## Value

A list of `gseaResult` objects containing the GSEA results

## Examples

``` r
if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(org.Mm.eg.db)
    data(example_obj)
    example_obj <- normalise_to_start(example_obj)

    var_decomp <- decomp_variance(example_obj,
        features = rownames(example_obj)[1:100], assay = 1)
    example_go_rank <- enrichGO_rank(var_decomp, gene_rank_by = "Time",
        OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
        category = "BP")
}
#> Normalising to group baseline at each feature's first non-NA time point.
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
#> Warning: There were 4907 pathways for which P-values were not calculated properly due to unbalanced gene-level statistic values. For such pathways pvalue, NES and log2err are set to NA. You can try to increase nPermSimple.
#> Warning: Invalid p-values detected (NA, non-finite, <0, or >1). qvalue will be computed on valid p-values only.
#> Warning: NA values detected in gene set IDs. Replacing with string 'NA'.
#> Warning: Duplicate gene set IDs detected: NA... (Total 1). Unique suffixes added.
#> Removing NA ID gene sets for BP.
```
