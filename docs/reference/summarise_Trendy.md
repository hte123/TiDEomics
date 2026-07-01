# Summarise Trendy results

Summarise Trendy results into data frame

## Usage

``` r
summarise_Trendy(res_list, ...)
```

## Arguments

- res_list:

  A list of Trendy analysis results, output of
  [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)

- ...:

  Additional arguments to be passed to the
  [`Trendy::topTrendy()`](https://rdrr.io/pkg/Trendy/man/topTrendy.html)

## Value

A data frame of summary results of Trendy, including breakpoints,
segment slopes and p-values for each fitted feature in each group. A
"Pattern" column is added to show the combination of segment trends (up,
down, stable) for each feature.

## Examples

``` r
data("example_res_list")
trendy_summary <- summarise_Trendy(example_res_list)
#> Warning: number of columns of result is not a multiple of vector length (arg 25)
#> Warning: number of columns of result is not a multiple of vector length (arg 19)
#> Warning: number of columns of result is not a multiple of vector length (arg 16)
```
