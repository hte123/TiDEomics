# Summarise WGCNA module metrics

Computes per-module metrics for a WGCNA network: module size, proportion
of features assigned, and module membership (kME). Unassigned module is
excluded. For consensus WGCNA output, metrics are returned for the
consensus and for each individual dataset.

## Usage

``` r
summarise_module_metrics(net)
```

## Arguments

- net:

  A WGCNA network object from
  [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)
  or
  [`run_consensus_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_consensus_WGCNA.md).

## Value

For an individual WGCNA network, a data frame with columns: `Module`,
`Size`, `Proportion`, `MeanKME`, `MeanKME2`, `MedianKME`, `SDKME`,
`Q25KME`, `Q75KME`, `MinKME`, and `MaxKME`. For consensus WGCNA, a named
list of these data frames, with `Consensus` followed by one element for
each individual dataset.

## Details

`MeanKME` is mean module membership
(`kME = cor(feature, module eigengene)`). `MeanKME2` is the mean squared
module membership (`mean(kME^2)`), equivalent to the mean per-feature
proportion of variance explained by the module eigengene. `Q25KME` and
`Q75KME` are the 25th and 75th kME percentiles. For consensus output, it
is calculated from consensus or individual-dataset kME scores,
respectively.

## Examples

``` r
data(example_net)
summarise_module_metrics(example_net)
#>   Module Size Proportion   MeanKME  MeanKME2 MedianKME     SDKME    Q25KME
#> 1     M1   52       0.52 0.6208447 0.4157510 0.6069187 0.1757754 0.4920063
#> 2     M2   16       0.16 0.6309279 0.4346083 0.6944016 0.1974188 0.4398500
#>      Q75KME    MinKME    MaxKME
#> 1 0.7716648 0.2908221 0.9580086
#> 2 0.7849205 0.3559238 0.8938520
data(example_cons)
summarise_module_metrics(example_cons)
#> $Consensus
#>   Module Size Proportion   MeanKME  MeanKME2 MedianKME     SDKME    Q25KME
#> 1     M1   51  0.5368421 0.5318533 0.3158869 0.5215082 0.1835195 0.3905486
#> 2     M2   12  0.1263158 0.6145119 0.4261337 0.6934061 0.2300405 0.3964295
#> 3     M3   11  0.1157895 0.5466791 0.3348633 0.4965516 0.1990121 0.4169232
#>      Q75KME    MinKME    MaxKME
#> 1 0.6747267 0.2091264 0.8369927
#> 2 0.7921272 0.2301306 0.8997157
#> 3 0.6628851 0.2559315 0.8820006
#> 
#> $Exp1
#>   Module Size Proportion   MeanKME  MeanKME2 MedianKME     SDKME    Q25KME
#> 1     M1   51  0.5368421 0.5991889 0.3954120 0.6009237 0.1926457 0.4685013
#> 2     M2   12  0.1263158 0.6714406 0.5004950 0.7758978 0.2327602 0.3964295
#> 3     M3   11  0.1157895 0.6001163 0.3933021 0.6369312 0.1909944 0.4542118
#>      Q75KME    MinKME    MaxKME
#> 1 0.7637917 0.2508767 0.9476520
#> 2 0.8465950 0.3111618 0.9316782
#> 3 0.7387240 0.3290971 0.8820006
#> 
#> $Exp2
#>   Module Size Proportion   MeanKME  MeanKME2 MedianKME     SDKME    Q25KME
#> 1     M1   51  0.5368421 0.5437077 0.3278069 0.5425846 0.1811978 0.3974038
#> 2     M2   12  0.1263158 0.6306196 0.4422440 0.6934061 0.2204859 0.5196878
#> 3     M3   11  0.1157895 0.5655531 0.3543330 0.4965516 0.1947587 0.4390390
#>      Q75KME    MinKME    MaxKME
#> 1 0.6819027 0.2091264 0.8467566
#> 2 0.7921272 0.2301306 0.8997157
#> 3 0.6930078 0.2559315 0.8882523
#> 
```
