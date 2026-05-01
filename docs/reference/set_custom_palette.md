# Set custom color palette

Set a custom color palette for groups to use in all subsequent plotting
functions that support custom palettes.

## Usage

``` r
set_custom_palette(palette)
```

## Arguments

- palette:

  A named vector where names correspond to group names and values are
  the corresponding colors (e.g.,
  `c("Group1" = "#FF0000", "Group2" = "#00FF00")`).

## Value

The function does not return a value but sets the custom palette as an
option that can be accessed by other plotting functions. The palette
will be stored in the R session options under the name "custom.palette".

## Examples

``` r
custom_palette <- c("untreated" = "#1b9e77", "IFNbeta" = "#d95f02",
    "IFNgamma" = "#7570b3", "LPS" = "#e7298a")
set_custom_palette(custom_palette)
#> Custom palette has been set successfully for groups untreated, IFNbeta, IFNgamma, LPS.
```
