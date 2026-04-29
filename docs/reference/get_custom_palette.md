# Get custom color palette

Get the custom color palette set by `set_custom_palette`. If no custom
palette has been set, the function will return a default color palette
based on the number of groups specified in the `groups` argument.

## Usage

``` r
get_custom_palette(groups)
```

## Arguments

- groups:

  A character vector of group names for which to retrieve the colors
  from the custom palette.

## Value

A named vector of colors corresponding to the specified group names.

## Examples

``` r
get_custom_palette(c("untreated", "IFNbeta"))
#> untreated   IFNbeta 
#> "#F8766D" "#00BFC4" 
```
