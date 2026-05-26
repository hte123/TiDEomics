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
data("example_res_list")
plot_breakpoints(example_res_list)
#> Warning: number of columns of result is not a multiple of vector length (arg 9)
#> Warning: number of columns of result is not a multiple of vector length (arg 3)
#> Warning: number of columns of result is not a multiple of vector length (arg 5)
```
