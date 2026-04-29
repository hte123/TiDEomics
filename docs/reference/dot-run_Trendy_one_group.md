# Segmented regression analysis (for one group)

Run segmented regression analysis with Trendy on imputed data for one
group of samples in a SummarizedExperiment object. For use in the
[`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)
function.

If less than 2 \* minNumInSeg time points are available, Trendy analysis
will not be performed and NULL will be returned.

## Usage

``` r
.run_Trendy_one_group(
  se_obj_imp,
  minExp = 0.5,
  feature = NULL,
  maxK = 1,
  meanCut = 0,
  minNumInSeg = 3,
  NCores = 2,
  ...
)
```

## Arguments

- se_obj_imp:

  A SummarizedExperiment object with no missing values

- minExp:

  Minimum expression ratio for a feature to be included in the analysis
  (default is 0.5)

- feature:

  A character vector of feature names to be included in the analysis. If
  NULL, all features with expression ratio (Exp_ratio) \>= minExp will
  be used, based on results of
  [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
  before imputation. (default is NULL)

- maxK:

  Parameter of
  [`Trendy::trendy()`](https://rdrr.io/pkg/Trendy/man/trendy.html),
  maximum number of breakpoints allowed in the segmented regression
  model (default is 1)

- meanCut:

  Parameter of
  [`Trendy::trendy()`](https://rdrr.io/pkg/Trendy/man/trendy.html),
  minimum mean expression required for a feature to be included in the
  analysis (default is 0)

- minNumInSeg:

  Parameter of
  [`Trendy::trendy()`](https://rdrr.io/pkg/Trendy/man/trendy.html),
  minimum number of samples required in each segment (default is 3)

- NCores:

  Number of cores to use for parallel processing (default is 2)

- ...:

  Additional arguments to be passed to the
  [`Trendy::trendy()`](https://rdrr.io/pkg/Trendy/man/trendy.html)

## Value

Trendy analysis result, including the fitted model parameters and
statistics for each feature.
