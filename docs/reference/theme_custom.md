# Custom ggplot2 theme

A minimal theme used by all TiDEomics plotting functions, built on
`theme_minimal()`. Exported so users can apply it to their own plots for
consistent styling.

## Usage

``` r
theme_custom(base_size = 8, panel_border = FALSE, legend_position = "right")
```

## Arguments

- base_size:

  Base font size (default: 8).

- panel_border:

  Logical, whether to draw a border around the panel (default: FALSE).

- legend_position:

  Legend position (default: "right").

## Value

A ggplot2 theme object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
    geom_point() +
    theme_custom(base_size = 10)
```
