# Calculate consensus and dataset-specific kME scores

Runs
[`WGCNA::consensusKME()`](https://rdrr.io/pkg/WGCNA/man/consensusKME.html)
and extracts assigned-module scores for each dataset and for the
consensus network.

## Usage

``` r
.consensus_kme_scores(consensus, exclude_grey = TRUE)
```

## Arguments

- consensus:

  Output from
  [`run_consensus_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_consensus_WGCNA.md).

- exclude_grey:

  Logical; whether to exclude the grey module.

## Value

A list containing module assignments, score vectors, and labels.
