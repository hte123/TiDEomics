# Summarise WGCNA module metrics

Computes per-module metrics for a WGCNA network: module size, proportion
of features assigned, and module membership (kME). Unassigned module is
excluded.

## Usage

``` r
summarise_module_metrics(net)
```

## Arguments

- net:

  A WGCNA network object from
  [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md).

## Value

A data frame with columns: `Module`, `Size`, `Proportion`, `MeanKME`,
`MeanKME2`, `MedianKME`, `SDKME`, `MinKME`, `MaxKME`.

## Details

`MeanKME` and `MeanKME2` are the mean and mean squared module membership
(kME = cor(feature, module eigengene)). `MeanKME2` = mean(kME^2) is the
average proportion of per-feature variance explained by the module
eigengene.

## Examples

``` r
data(example_net)
summarise_module_metrics(example_net)
#>   Module Size Proportion   MeanKME  MeanKME2 MedianKME     SDKME    MinKME
#> 1     M1   52       0.52 0.6208447 0.4157510 0.6069187 0.1757754 0.2908221
#> 2     M2   16       0.16 0.6309279 0.4346083 0.6944016 0.1974188 0.3559238
#>      MaxKME
#> 1 0.9580086
#> 2 0.8938520
```
