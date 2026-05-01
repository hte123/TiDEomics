# Plot segmented regression

Plot segmented regression results for specified features in each group,
using the
[`Trendy::plotFeature()`](https://rdrr.io/pkg/Trendy/man/plotFeature.html)
function. The input is a list of SummarizedExperiment objects with no
missing values and a list of Trendy analysis results for each group.

## Usage

``` r
plot_segments(
  se_obj_imp_list,
  res_list,
  feature,
  group = NULL,
  nrow = NULL,
  ...
)
```

## Arguments

- se_obj_imp_list:

  A list of SummarizedExperiment objects with no missing values

- res_list:

  A list of the Trendy analysis results, output of
  [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)

- feature:

  A character vector of feature names to be plotted

- group:

  A character vector of group names to be analysed. If NULL, all groups
  in the input list will be used (default is NULL)

- nrow:

  Number of rows in the plot layout. If NULL and multiple features are
  provided, it will be set to half the number of features (rounded up)
  (default is NULL)

- ...:

  Additional arguments to be passed to the
  [`Trendy::plotFeature()`](https://rdrr.io/pkg/Trendy/man/plotFeature.html)

## Value

Plots of the segmented regression fits for the specified features in
each group.

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged_list <-
    calc_feature_property(example_obj_merged_list, threshold = 0)
# no missing value in the example dataset, so imputation is not necessary
example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)

# "untreated" group has only 3 time points, so Trendy analysis will not be
# performed for this group
# example_res_list <- run_Trendy(example_obj_merged_imp_list, maxK = 1,
#     minNumInSeg = 2, meanCut = 0)
data("example_res_list")

plot_segments(example_obj_merged_imp_list, example_res_list,
    feature = c("Mctp1"))
#> Plotting segmented regression for group: IFNbeta

#> Plotting segmented regression for group: IFNgamma

#> Plotting segmented regression for group: LPS

plot_breakpoints(example_res_list)
#> Warning: number of columns of result is not a multiple of vector length (arg 9)
#> Warning: number of columns of result is not a multiple of vector length (arg 3)
#> Warning: number of columns of result is not a multiple of vector length (arg 5)

trendy_summary <- summarise_Trendy(example_res_list)
#> Warning: number of columns of result is not a multiple of vector length (arg 9)
#> Warning: number of columns of result is not a multiple of vector length (arg 3)
#> Warning: number of columns of result is not a multiple of vector length (arg 5)
trendy_list <- extract_segment_trends(trendy_summary)
```
