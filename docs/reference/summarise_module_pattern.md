# Summarise module patterns

Summarise module patterns by integrating WGCNA output with Trendy
results

## Usage

``` r
summarise_module_pattern(module, trendy_summary, print_top_n = TRUE, top_n = 5)
```

## Arguments

- module:

  A data frame with columns "Feature" and "Module"

- trendy_summary:

  A data frame of trendy results, output of
  [`summarise_Trendy()`](https://hte123.github.io/TiDEomics/reference/summarise_Trendy.md)

- print_top_n:

  Whether to output the top n patterns per module as text (default is
  TRUE)

- top_n:

  Number of top patterns to show per module when `print_top_n` is TRUE
  (default is 5)

## Value

A list of data frames, each data frame shows the pattern counts in a
module

## Examples

``` r
data("example_res_list")
trendy_summary <- summarise_Trendy(example_res_list)
#> Warning: number of columns of result is not a multiple of vector length (arg 9)
#> Warning: number of columns of result is not a multiple of vector length (arg 3)
#> Warning: number of columns of result is not a multiple of vector length (arg 5)

data(example_net)
example_module <- WGCNA_module(example_net) 
summarise_module_pattern(example_module, trendy_summary)
#> Most common pattern in each module (IFNbeta, IFNgamma, LPS):
#> Module 1: down_up, NA, NA
#> Module 2: NA, NA, stable_stable
#> Module 3: NA, NA, stable_stable
#> Module 1 top 5 patterns:
#>       IFNbeta, IFNgamma, LPS Count
#> 1            down_up, NA, NA     5
#> 2 down_up, NA, stable_stable     2
#> 3   down_up, down_stable, NA     2
#> 4            NA, down_up, NA     1
#> 5          NA, up, stable_up     1
#> Module 2 top 5 patterns:
#>             IFNbeta, IFNgamma, LPS Count
#> 1            NA, NA, stable_stable     2
#> 2            stable_stable, NA, NA     2
#> 3            NA, stable_stable, NA     1
#> 4                     down, NA, NA     1
#> 5 stable_stable, down_stable, down     1
#> Module 3 top 5 patterns:
#>         IFNbeta, IFNgamma, LPS Count
#> 1        NA, NA, stable_stable     1
#> 2                 NA, down, NA     1
#> 3 NA, stable_up, stable_stable     1
#> 4        down, NA, down_stable     1
#> 5      down, NA, stable_stable     1
#> $`1`
#>                     IFNbeta, IFNgamma, LPS Count
#> 1                          down_up, NA, NA     5
#> 2               down_up, NA, stable_stable     2
#> 3                 down_up, down_stable, NA     2
#> 4                          NA, down_up, NA     1
#> 5                        NA, up, stable_up     1
#> 6         down_stable, down, stable_stable     1
#> 7    down_stable, down_stable, down_stable     1
#> 8  down_stable, stable_stable, down_stable     1
#> 9            down_up, down_up, down_stable     1
#> 10              stable_stable, down_up, NA     1
#> 11            stable_stable, stable_up, NA     1
#> 12                       stable_up, NA, NA     1
#> 13              stable_up, NA, down_stable     1
#> 14            stable_up, stable_stable, NA     1
#> 
#> $`2`
#>             IFNbeta, IFNgamma, LPS Count
#> 1            NA, NA, stable_stable     2
#> 2            stable_stable, NA, NA     2
#> 3            NA, stable_stable, NA     1
#> 4                     down, NA, NA     1
#> 5 stable_stable, down_stable, down     1
#> 
#> $`3`
#>               IFNbeta, IFNgamma, LPS Count
#> 1              NA, NA, stable_stable     1
#> 2                       NA, down, NA     1
#> 3       NA, stable_up, stable_stable     1
#> 4              down, NA, down_stable     1
#> 5            down, NA, stable_stable     1
#> 6                   down, down, down     1
#> 7            down, down, down_stable     1
#> 8         down, down_down, down_down     1
#> 9         down_down, down, down_down     1
#> 10 down_stable, down_stable, down_up     1
#> 11      stable_down, NA, down_stable     1
#> 12  stable_stable, NA, stable_stable     1
#> 
```
