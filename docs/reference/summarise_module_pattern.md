# Summarise module patterns

Summarise module patterns by integrating WGCNA output with Trendy
results

## Usage

``` r
summarise_module_pattern(module, trendy_summary)
```

## Arguments

- module:

  A data frame with columns "Feature" and "Module"

- trendy_summary:

  A data frame of trendy results, output of
  [`summarise_Trendy()`](https://hte123.github.io/TiDEomics/reference/summarise_Trendy.md)

## Value

A list of data frames, each data frame shows the pattern counts in a
module

## Examples

``` r
data("example_res_list")
trendy_summary <- summarise_Trendy(example_res_list)
#> Warning: number of columns of result is not a multiple of vector length (arg 25)
#> Warning: number of columns of result is not a multiple of vector length (arg 19)
#> Warning: number of columns of result is not a multiple of vector length (arg 16)

data(example_net)
example_module <- WGCNA_module(example_net)
summarise_module_pattern(example_module, trendy_summary)
#> Most common pattern in each module (IFNbeta, IFNgamma, LPS):
#> Module 1: stable_up, NA, NA
#> Module 2: NA, up_stable, down_stable
#> $`1`
#> # A tibble: 36 × 2
#>    Pattern                        Count
#>    <chr>                          <int>
#>  1 stable_up, NA, NA                  3
#>  2 down_stable, NA, stable_stable     2
#>  3 down_up, NA, NA                    2
#>  4 down_up, down_stable, down_up      2
#>  5 NA, NA, stable_stable              1
#>  6 NA, NA, stable_up                  1
#>  7 NA, down, NA                       1
#>  8 NA, down_stable, down              1
#>  9 NA, down_stable, stable_stable     1
#> 10 NA, down_up, stable_stable         1
#> # ℹ 26 more rows
#> 
#> $`2`
#> # A tibble: 13 × 2
#>    Pattern                               Count
#>    <chr>                                 <int>
#>  1 NA, up_stable, down_stable                1
#>  2 stable_stable, NA, NA                     1
#>  3 up, NA, stable_stable                     1
#>  4 up, stable_up, NA                         1
#>  5 up_down, NA, NA                           1
#>  6 up_down, NA, stable_down                  1
#>  7 up_down, down, NA                         1
#>  8 up_down, stable_down, NA                  1
#>  9 up_down, stable_stable, stable_stable     1
#> 10 up_down, up_down, stable_stable           1
#> 11 up_stable, NA, stable_stable              1
#> 12 up_stable, up_stable, NA                  1
#> 13 up_stable, up_stable, down_up             1
#> 
```
