# Flatten nested enrichment results

Flattens the nested enrichment result structure returned by
[`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md),
[`enrichR_list()`](https://hte123.github.io/TiDEomics/reference/enrichR_list.md),
or
[`enrich_msigdb()`](https://hte123.github.io/TiDEomics/reference/enrich_msigdb.md)
into a flat named list of data frames or `enrichResult` objects. This is
useful for exporting results to other tools such as `DeeDeeExperiment`.

## Usage

``` r
flatten_enrich(enrich_list)
```

## Arguments

- enrich_list:

  Output of
  [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md),
  [`enrichR_list()`](https://hte123.github.io/TiDEomics/reference/enrichR_list.md),
  or
  [`enrich_msigdb()`](https://hte123.github.io/TiDEomics/reference/enrich_msigdb.md);
  or a named list combining several such outputs (e.g.
  `list(GO = go_res, MSigDB = msigdb_res)`).

## Value

A flat named list of data frames or `enrichResult` objects with
`DeeDeeExperiment`-compatible names and columns (e.g.
"clusterProfiler_BP_1", "enrichr_DSigDB"). Returns an empty list if
`enrich_list` is empty.

## Examples

``` r
data(example_go)
enrich_flat <- flatten_enrich(example_go)
str(enrich_flat, max.level = 1)
#> List of 4
#>  $ clusterProfiler_BP_1:Formal class 'enrichResult' [package "enrichit"] with 15 slots
#>  $ clusterProfiler_BP_2:Formal class 'enrichResult' [package "enrichit"] with 15 slots
#>  $ clusterProfiler_CC_1:Formal class 'enrichResult' [package "enrichit"] with 15 slots
#>  $ clusterProfiler_CC_2:Formal class 'enrichResult' [package "enrichit"] with 15 slots
```
