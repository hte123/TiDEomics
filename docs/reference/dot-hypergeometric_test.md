# Hypergeometric enrichment test

One-tailed hypergeometric (Fisher's exact) test of query genes against
named gene sets.

## Usage

``` r
.hypergeometric_test(query, gene_sets, universe, min_overlap = 1L)
```

## Arguments

- query:

  Character vector of query genes.

- gene_sets:

  Named list of gene vectors (gene set name -\> genes, already filtered
  to universe).

- universe:

  Character vector of all background genes.

- min_overlap:

  Minimum number of overlapping genes to report (default: 1).

## Value

A data.frame with columns `Term`, `N`, `m`, `k`, `q`, `pvalue`, or NULL
if no gene set passes filters.
