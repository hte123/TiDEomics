# Convert a TiDEomics result to DeeDeeExperiment

Wraps the output of
[`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md)
into a `DeeDeeExperiment` object.

## Usage

``` r
as_DeeDeeExperiment(tide)
```

## Arguments

- tide:

  A named list returned by
  [`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md).

## Value

A `DeeDeeExperiment` object if the `DeeDeeExperiment` package is
available, otherwise an error.

## Examples

``` r
data(tutorial_data)
data(tutorial_sample_info)
tide <- prepare_tide(tutorial_data, tutorial_sample_info,
    keep = "threshold", residual_threshold = 100)
#> No Subject column specified. Samples treated as independent. For repeated-measures designs, set subject_col to the column identifying biological subjects.
#> Converting 'Group' column to factor. Default order is alphabetical.
#> Converting 'Replicate' column to factor. Default order is numerical.
#> Converting 'Batch' column to factor. Default order is numerical.
#> Normalising to group baseline at each feature's first non-NA time point.
#> Preparing TiDEomics input: 500 features, 40 samples, 4 groups
#> Filtering criteria: >=50% values >0 in >=2 of groups: IFNbeta, IFNgamma, LPS, untreated
#> Cross-group filter (Exp_ratio >= 0.5 in >= 2 groups): kept 340 of 500 features (68%)
#> --- assay: orig ---
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
#> Residual filter (threshold): kept 334 of 340 features (98.2%)
#> --- assay: norm ---
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
#> Residual filter (threshold): kept 333 of 340 features (97.9%)
if (requireNamespace("DeeDeeExperiment", quietly = TRUE)) {
    dde <- as_DeeDeeExperiment(tide)
}
#> Warning: replacing previous import 'BiocGenerics::transform' by 'IRanges::transform' when loading 'DESeq2'
```
