# Prepare WGCNA input data

Prepare input data for WGCNA, including transposing the data to have
samples in rows and features in columns, and removing bad samples and
features. Features can be pre-filtered, e.g. by residual variance
calculated by
[`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md),
to remove noisy features before preparing the data for WGCNA.

## Usage

``` r
WGCNA_input(se_obj, assay)
```

## Arguments

- se_obj:

  A SummarizedExperiment object. Data normalised to time point 0 can be
  in the second assay slot, created by
  [`normalise_to_start()`](https://hte123.github.io/TiDEomics/reference/normalise_to_start.md).

- assay:

  Which assay slot of the SummarizedExperiment object to use for WGCNA
  input.

## Value

A data frame with samples in rows and features in columns, filtered to
remove bad samples and features, ready for use in
[`prepare_WGCNA()`](https://hte123.github.io/TiDEomics/reference/prepare_WGCNA.md).
