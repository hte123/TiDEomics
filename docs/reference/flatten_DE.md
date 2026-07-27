# Flatten nested differential expression results

Flattens the nested DE result structure returned by
[`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md)
or
[`DE_between_time()`](https://hte123.github.io/TiDEomics/reference/DE_between_time.md)
into a named list of data frames, one entry per contrast-time
combination. This is useful for exporting results to other tools such as
`DeeDeeExperiment`.

## Usage

``` r
flatten_DE(de_list)
```

## Arguments

- de_list:

  A list of DE results. Accepts individual output of
  [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md)
  or
  [`DE_between_time()`](https://hte123.github.io/TiDEomics/reference/DE_between_time.md)
  with an `all_list` element; or a named list (e.g.
  `list(between_group = ..., between_time = ...)`), each containing an
  `all_list` element.

## Value

A named list of data frames with columns renamed for `DeeDeeExperiment`
compatibility (`log2FoldChange`, `pvalue`, `padj`). Returns an empty
list if `de_list` is empty.

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
tide$DE <- DE_between_group(tide$se, assay = "norm",
    filter = 1, trend = TRUE)
#> Comparing group IFNgamma to IFNbeta at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to IFNbeta at Time 2: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to IFNbeta at Time 4: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to IFNbeta at Time 6: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to IFNbeta at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to IFNbeta at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNbeta at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNbeta at Time 2: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNbeta at Time 4: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNbeta at Time 6: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNbeta at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNbeta at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing untreated vs IFNbeta: time points only in IFNbeta: 2, 4, 6; only in untreated: none
#> Comparing group untreated to IFNbeta at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated to IFNbeta at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated to IFNbeta at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to IFNgamma at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to IFNgamma at Time 2: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to IFNgamma at Time 4: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to IFNgamma at Time 6: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to IFNgamma at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to IFNgamma at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNgamma at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNgamma at Time 2: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNgamma at Time 4: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNgamma at Time 6: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNgamma at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to IFNgamma at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing untreated vs IFNgamma: time points only in IFNgamma: 2, 4, 6; only in untreated: none
#> Comparing group untreated to IFNgamma at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated to IFNgamma at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated to IFNgamma at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to LPS at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to LPS at Time 2: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to LPS at Time 4: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to LPS at Time 6: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to LPS at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to LPS at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to LPS at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to LPS at Time 2: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to LPS at Time 4: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to LPS at Time 6: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to LPS at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to LPS at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing untreated vs LPS: time points only in LPS: 2, 4, 6; only in untreated: none
#> Comparing group untreated to LPS at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated to LPS at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated to LPS at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing IFNbeta vs untreated: time points only in untreated: none; only in IFNbeta: 2, 4, 6
#> Comparing group IFNbeta to untreated at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to untreated at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta to untreated at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing IFNgamma vs untreated: time points only in untreated: none; only in IFNgamma: 2, 4, 6
#> Comparing group IFNgamma to untreated at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to untreated at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma to untreated at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing LPS vs untreated: time points only in untreated: none; only in LPS: 2, 4, 6
#> Comparing group LPS to untreated at Time 0: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to untreated at Time 8: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS to untreated at Time 24: keeping 500 of 500 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
de_flat <- flatten_DE(tide$DE)
str(de_flat, max.level = 1)
#> List of 54
#>  $ IFNgamma-IFNbeta_T0   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-IFNbeta_T2   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-IFNbeta_T4   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-IFNbeta_T6   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-IFNbeta_T8   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-IFNbeta_T24  :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNbeta_T0        :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNbeta_T2        :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNbeta_T4        :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNbeta_T6        :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNbeta_T8        :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNbeta_T24       :'data.frame':  500 obs. of  11 variables:
#>  $ untreated-IFNbeta_T0  :'data.frame':  500 obs. of  12 variables:
#>  $ untreated-IFNbeta_T8  :'data.frame':  500 obs. of  11 variables:
#>  $ untreated-IFNbeta_T24 :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-IFNgamma_T0   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-IFNgamma_T2   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-IFNgamma_T4   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-IFNgamma_T6   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-IFNgamma_T8   :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-IFNgamma_T24  :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNgamma_T0       :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNgamma_T2       :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNgamma_T4       :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNgamma_T6       :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNgamma_T8       :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-IFNgamma_T24      :'data.frame':  500 obs. of  11 variables:
#>  $ untreated-IFNgamma_T0 :'data.frame':  500 obs. of  12 variables:
#>  $ untreated-IFNgamma_T8 :'data.frame':  500 obs. of  11 variables:
#>  $ untreated-IFNgamma_T24:'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-LPS_T0        :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-LPS_T2        :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-LPS_T4        :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-LPS_T6        :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-LPS_T8        :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-LPS_T24       :'data.frame':  500 obs. of  11 variables:
#>  $ IFNgamma-LPS_T0       :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-LPS_T2       :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-LPS_T4       :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-LPS_T6       :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-LPS_T8       :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-LPS_T24      :'data.frame':  500 obs. of  11 variables:
#>  $ untreated-LPS_T0      :'data.frame':  500 obs. of  12 variables:
#>  $ untreated-LPS_T8      :'data.frame':  500 obs. of  11 variables:
#>  $ untreated-LPS_T24     :'data.frame':  500 obs. of  11 variables:
#>  $ IFNbeta-untreated_T0  :'data.frame':  500 obs. of  12 variables:
#>  $ IFNbeta-untreated_T8  :'data.frame':  500 obs. of  11 variables:
#>  $ IFNbeta-untreated_T24 :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-untreated_T0 :'data.frame':  500 obs. of  12 variables:
#>  $ IFNgamma-untreated_T8 :'data.frame':  500 obs. of  11 variables:
#>  $ IFNgamma-untreated_T24:'data.frame':  500 obs. of  12 variables:
#>  $ LPS-untreated_T0      :'data.frame':  500 obs. of  12 variables:
#>  $ LPS-untreated_T8      :'data.frame':  500 obs. of  11 variables:
#>  $ LPS-untreated_T24     :'data.frame':  500 obs. of  11 variables:
```
