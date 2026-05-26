# Prepare input gene list for enrichment functions

Convert a WGCNA_module() output data.frame to a named gene list, or pass
through an existing named list unchanged.

## Usage

``` r
.prepare_gene_list(x)
```

## Arguments

- x:

  A data.frame with 'Feature' and 'Module' columns (as returned by
  WGCNA_module()), or a named list of gene vectors.

## Value

A named list of gene vectors, where names correspond to module names.
