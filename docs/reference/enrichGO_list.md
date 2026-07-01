# GO enrichment with gene sets

Gene ontology enrichment analysis of multiple gene sets with
clusterProfiler, allowing for using different or same background genes
for each gene set

## Usage

``` r
enrichGO_list(
  gene_list,
  keyType = "SYMBOL",
  OrgDb,
  universe = NULL,
  universe_list = NULL,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  category = NULL,
  simplify = FALSE,
  simplify_cutoff = 0.7,
  simplify_by = "p.adjust",
  simplify_select_fun = min,
  simplify_measure = "Wang",
  ...
)
```

## Arguments

- gene_list:

  A named list of gene vectors, or a data.frame with `Feature` and
  `Module` columns from
  [`WGCNA_module()`](https://hte123.github.io/TiDEomics/reference/WGCNA_module.md).

- keyType:

  (Optional) Available options are `AnnotationDbi::keytypes(OrgDb)`
  (default is "SYMBOL")

- OrgDb:

  Organism database, e.g. org.Hs.eg.db, org.Mm.eg.db

- universe:

  Background genes for all input gene sets, used if `universe_list` is
  not provided

- universe_list:

  Background genes for each input gene set, a list of gene vectors with
  the same names as `gene_list`

- pAdjustMethod:

  (Optional) Parameter of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)
  (default is "BH")

- pvalueCutoff:

  (Optional) Parameter of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)
  (default is 0.05)

- qvalueCutoff:

  (Optional) Parameter of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)
  (default is 0.05)

- category:

  (Optional) GO category to analyze (default is all three of BP, MF, CC)

- simplify:

  (Optional) Whether to simplify the GO terms by removing redundant
  terms with
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html)
  function. (default is FALSE)

- simplify_cutoff:

  (Optional) Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  cutoff for similarity when simplifying GO terms (default is 0.7)

- simplify_by:

  (Optional) Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  method to choose representative term when simplifying GO terms
  (default is "p.adjust")

- simplify_select_fun:

  (Optional) Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  function to select representative term when simplifying GO terms
  (default is `min`)

- simplify_measure:

  (Optional) Parameter of
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html),
  method to calculate similarity when simplifying GO terms (default is
  "Wang")

- ...:

  additional arguments passed to
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)

## Value

A nested list of GO enrichment results with sublists: 'all' and
'simplified' (if `simplify = TRUE`) containing merged results of all or
simplified terms across gene sets for each GO category; 'unmerged_all'
and 'unmerged_simplified' (if `simplify = TRUE`) including all or
simplified terms for each gene set and GO category.

## Examples

``` r
if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(org.Mm.eg.db)
    library(clusterProfiler)
    data(example_net)
    # select two modules for demonstration
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))
    # set cutoff to 1 to show all results for demonstration
    example_go_list = enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)
}
#> Warning: replacing previous import 'BiocGenerics::transform' by 'S4Vectors::transform' when loading 'AnnotationDbi'
#> Warning: replacing previous import 'utils::data' by 'BiocGenerics::data' when loading 'Biostrings'
#> Warning: replacing previous import 'BiocGenerics::transform' by 'S4Vectors::transform' when loading 'Biostrings'
#> 
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> Warning: package 'generics' was built under R version 4.6.1
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from 'package:utils':
#> 
#>     data
#> The following objects are masked from 'package:base':
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     I, expand.grid, unname
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:grDevices':
#> 
#>     windows
#> 
#> clusterProfiler v4.21.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
#> 
#> Please cite:
#> 
#> S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang,
#> W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize
#> multiomics data. Nature Protocols. 2024, 19(11):3292-3320
#> 
#> Attaching package: 'clusterProfiler'
#> The following object is masked from 'package:AnnotationDbi':
#> 
#>     select
#> The following object is masked from 'package:IRanges':
#> 
#>     slice
#> The following object is masked from 'package:S4Vectors':
#> 
#>     rename
#> The following object is masked from 'package:stats':
#> 
#>     filter
#> Performing GO enrichment for category: BP
#> Processing gene list: 1
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> Processing gene list: 2
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> Merging GO enrichment results across gene lists for each category.
```
