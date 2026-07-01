# Extract feature trends

Extract features with different segment trends from trendy segmented
regression results summary

## Usage

``` r
extract_segment_trends(trendy.summary)
```

## Arguments

- trendy.summary:

  A data frame, output of
  [`summarise_Trendy()`](https://hte123.github.io/TiDEomics/reference/summarise_Trendy.md),
  containing the segmented regression results for all features in all
  groups.

## Value

A nested list of features grouped by their segment trend patterns within
each group.

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

data("example_res_list")
trendy_summary <- summarise_Trendy(example_res_list)
#> Warning: number of columns of result is not a multiple of vector length (arg 25)
#> Warning: number of columns of result is not a multiple of vector length (arg 19)
#> Warning: number of columns of result is not a multiple of vector length (arg 16)
trendy_list <- extract_segment_trends(trendy_summary)
```
