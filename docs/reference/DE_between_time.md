# DE between time points

Differential expression analysis between time points within each group
by limma

## Usage

``` r
DE_between_time(
  se_obj,
  group = NULL,
  filter = NULL,
  assay,
  adjP_thres = 0.05,
  logFC_thres = 1,
  trend = FALSE
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)

- group:

  A character vector specifying which groups to analyse. If NULL, all
  groups in the 'Group' column will be used. (default is NULL)

- filter:

  Minimum number of replicates required in both conditions for a feature
  to be tested. If NULL, the minimum number of replicates across all
  groups and time points will be used. (default is NULL)

- assay:

  The assay to use in the SummarizedExperiment object: a numeric index
  or character name, e.g. 1 or "orig" for original data, 2 or "norm" for
  time 0 normalised data. The selected assay should contain
  log-transformed, normalised values. A warning is issued if the data
  appears to be un-logged raw counts.

- adjP_thres:

  Threshold for adjusted p-value to consider a feature as differentially
  expressed (default is 0.05)

- logFC_thres:

  Threshold for log2 fold change to consider a feature as differentially
  expressed (default is 1)

- trend:

  Logical, passed to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html). Use
  `TRUE` only when the mean-variance relationship is expected to vary
  with abundance. The default `FALSE` is usually appropriate for
  log-transformed, normalised data unless diagnostics support modelling
  a mean-variance trend.

## Value

A list with: `all_list` (nested list of DE results per group and time
comparison, all features); `de_list` (nested list of significant
features only); `fit_list` (nested list of limma `MArrayLM` fit objects,
for use with
[`limma::plotSA()`](https://rdrr.io/pkg/limma/man/plotSA.html));
`time_series` (vector of all available time points). The output can be
passed to
[`plot_DE_between_time()`](https://hte123.github.io/TiDEomics/reference/plot_DE_between_time.md)
for visualisation.

## Examples

``` r
data(example_obj)
DE_between_time_out <- DE_between_time(example_obj, assay = 1)
#> Non-NA replicate number filter not specified. Using minimum number of replicates across all groups and time points: 0
#> Comparing group IFNbeta time 2 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 4 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 6 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 8 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 24 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 4 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 6 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 8 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 24 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 6 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 8 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 24 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 8 to 6: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 24 to 6: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNbeta time 24 to 8: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 2 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 4 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 6 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 8 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 24 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 4 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 6 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 8 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 24 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 6 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 8 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 24 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 8 to 6: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 24 to 6: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group IFNgamma time 24 to 8: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 2 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 4 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 6 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 8 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 24 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 4 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 6 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 8 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 24 to 2: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 6 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 8 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 24 to 4: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 8 to 6: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 24 to 6: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group LPS time 24 to 8: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated time 8 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated time 24 to 0: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Comparing group untreated time 24 to 8: keeping 100 of 100 features (100.0%)
#> Warning: Zero sample variances detected, have been offset away from zero
```
