# Plot breakpoint distribution

Plot breakpoint distribution among time points for each group

## Usage

``` r
plot_breakpoints(res_list, group = NULL, fontsize = 8, ...)
```

## Arguments

- res_list:

  A list of the Trendy analysis results, output of
  [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)

- group:

  (Optional) A character vector of group names to be plotted. If NULL,
  all groups in the input list will be used (default is NULL)

- fontsize:

  (Optional) Font size for the plot (default is 8)

- ...:

  Additional arguments to be passed to
  [`Trendy::topTrendy()`](https://rdrr.io/pkg/Trendy/man/topTrendy.html)

## Value

A plot showing the distribution of breakpoints over time for each
specified group

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
