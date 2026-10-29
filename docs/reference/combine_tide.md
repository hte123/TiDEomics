# Combine multiple `prepare_tide()` outputs

Combine a list of
[`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md)
output objects for multi-dataset or multi-omics analysis.

## Usage

``` r
combine_tide(tide_list, assay = 1)
```

## Arguments

- tide_list:

  A non-empty named list of results from
  [`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md).

- assay:

  The assay to use: a numeric index or character name, e.g. 1 or "orig"
  for original data, 2 or "norm" for time 0 normalised data.

## Value

A list containing the filtered and merged SE objects (all features and
shared features), variance decomposition, shared feature names, the
resolved assay name, and a per-dataset filtering summary.

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
#> 340 features evaluated: 267 successful; 0 errors; 73 warnings; 0 skipped before fitting for insufficient observations or variance.
#> Residual filter (threshold): kept 334 of 340 features (98.2%)
#> --- assay: norm ---
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
#> 340 features evaluated: 267 successful; 0 errors; 73 warnings; 0 skipped before fitting for insufficient observations or variance.
#> Residual filter (threshold): kept 333 of 340 features (97.9%)
combined <- combine_tide(list(RNA = tide, Protein = tide), assay = "norm")
#> Combining 2 prepare_tide() results using assay 'norm'
#> Shared filtered features: 333 (RNA=333, Protein=333)
combined$filter_summary
#>   Dataset n_total n_kept pct_kept residual_cutoff
#> 1     RNA     500    333     66.6             100
#> 2 Protein     500    333     66.6             100
```
