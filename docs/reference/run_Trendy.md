# Segmented regression analysis

Run segmented regression analysis with Trendy on imputed data for
multiple groups of samples in a list of SummarizedExperiment objects.
The results will be returned in a list format. Each element in the list
corresponds to one group of samples and contains the Trendy analysis
results for that group.

If the `feature` parameter is not specified, all features expressed in
at least `minExp` time points in each group will be included in the
analysis for the group.
[`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
should be run before imputation to calculate the expression ratio for
each feature in each group. Missing values can be imputed after
[`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
using the
[`impute_groups()`](https://hte123.github.io/TiDEomics/reference/impute_groups.md)
function.

If a group has less than 4 time points are available, Trendy analysis
cannot be performed and NULL will be returned for the group.

## Usage

``` r
run_Trendy(
  se_obj_imp_list,
  group = NULL,
  feature = NULL,
  minExp = 0.5,
  maxK = 1,
  meanCut = 0,
  minNumInSeg = 3,
  NCores = 1,
  ...
)
```

## Arguments

- se_obj_imp_list:

  A list of SummarizedExperiment objects with no missing values

- group:

  A character vector of group names to be analysed. If NULL, all groups
  in the input list will be used (default is NULL)

- feature:

  A character vector of feature names to be included in the analysis. If
  NULL, all features with expression ratio (Exp_ratio) \>= minExp will
  be used, based on results of
  [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
  before imputation. (default is NULL)

- minExp:

  Minimum expression ratio for a feature to be included in the analysis
  (default is 0.5)

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

  Number of cores to use for parallel processing (default is 1)

- ...:

  Additional arguments to be passed to the
  [`Trendy::trendy()`](https://rdrr.io/pkg/Trendy/man/trendy.html)

## Value

A list of Trendy analysis results, including the fitted model parameters
and statistics for each feature in each group.

## Examples

``` r
data(example_obj)
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged_list <- calc_feature_property(example_obj_merged_list,
    threshold = 0)
# no missing value in the example dataset, so imputation is not necessary
example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
#> Group IFNbeta: no missing values.
#> Group IFNgamma: no missing values.
#> Group LPS: no missing values.
#> Group untreated: no missing values.

# "untreated" group has only 3 time points, so Trendy analysis will not be
# performed for this group
example_res_list <- run_Trendy(example_obj_merged_imp_list, maxK = 1,
    minNumInSeg = 2, meanCut = 0,
    feature = rownames(example_obj_merged_imp_list[[1]])[1:10])
#> Max number of breakpoints: 1
#> Min mean expression: 0
#> Min number of samples in each segment: 2
#> Running Trendy for group: IFNbeta
#> Using 10 specified features present in the data.
#> Warning: MulticoreParam() not supported on Windows, use SnowParam()
#> breakpoint estimate(s): 9.519174 
#> Running Trendy for group: IFNgamma
#> Using 10 specified features present in the data.
#> Warning: MulticoreParam() not supported on Windows, use SnowParam()
#> breakpoint estimate(s): 9.519174 
#> Running Trendy for group: LPS
#> Using 10 specified features present in the data.
#> Warning: MulticoreParam() not supported on Windows, use SnowParam()
#> breakpoint estimate(s): 8.885903 
#> Running Trendy for group: untreated
#> Trendy analysis is not performed for group untreated: number of time points (3) less than required ((maxK + 1) * minNumInSeg = 4
```
