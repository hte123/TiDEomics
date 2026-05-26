# Gene-set enrichment via enrichR

Enrichment analysis of multiple gene sets against databases from enrichR
(e.g. DSigDB for drug signatures, ChEA for TF targets, KEGG for
pathways, DrugBank for drug targets).

## Usage

``` r
enrichR_list(
  gene_list,
  databases = c("DSigDB", "DrugMatrix"),
  site = "Enrichr",
  universe = NULL,
  universe_list = NULL,
  pvalueCutoff = 0.05
)
```

## Arguments

- gene_list:

  A named list of gene vectors, or a data.frame with `Feature` and
  `Module` columns from
  [`WGCNA_module()`](https://hte123.github.io/TiDEomics/reference/WGCNA_module.md).
  Gene identifiers must match the convention of the chosen site (e.g.
  human gene symbols for Enrichr, fly symbols for FlyEnrichr).

- databases:

  Character vector of enrichR database names to query. Available
  databases can be checked with
  [`enrichR::listEnrichrDbs()`](https://rdrr.io/pkg/enrichR/man/listEnrichrDbs.html).
  (default: `c("DSigDB", "DrugMatrix")`)

- site:

  enrichR site to query. See
  [`enrichR::listEnrichrSites()`](https://rdrr.io/pkg/enrichR/man/listEnrichrSites.html).
  (default: `"Enrichr"`)

- universe:

  Background genes for all input gene sets, used if `universe_list` is
  not provided.

- universe_list:

  Background genes for each input gene set, a list of gene vectors with
  the same names as `gene_list`.

- pvalueCutoff:

  Adjusted p-value cutoff for filtering enriched terms (default: 0.05).

## Value

A named list of data.frames, one per database. Each data.frame has
columns `Cluster`, `Description`, `p.adjust` (Adjusted.P.value from
enrichR output), `Odds.Ratio`, `Combined.Score`, `Genes`, and any
additional columns returned by the enrichR API. Compatible with
`plot_modules_h(enrich_list = result, enrich_category = "DSigDB")`.

## Details

**Multiple testing:** P-values are adjusted (Benjamini–Hochberg)
independently for each module within each database. No correction is
applied across databases or across modules.

Requires an internet connection.

## Examples

``` r
data(example_net)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:AnnotationDbi':
#> 
#>     select
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from 'package:generics':
#> 
#>     explain
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
example_module <- WGCNA_module(example_net) %>%
    dplyr::filter(Module %in% c("1", "2"))
# Use high pvalueCutoff for demonstration
# enrichr_out <- enrichR_list(example_module, 
#     databases = c("KEGG_2019_Mouse"),
#     universe = example_module$Feature, pvalueCutoff = 0.5)
```
