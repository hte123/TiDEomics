# Plot PCA by group (one object)

Plot PCA for each group separately, with optional circles around time
points and arrows indicating trajectory over time.

## Usage

``` r
plot_pca_by_group(
  se_obj,
  circle = TRUE,
  arrow = TRUE,
  nrow = 1,
  fontsize = 8,
  assay = 1
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)

- circle:

  Logical, whether to draw circles (ellipses) around samples of each
  time point (default is TRUE)

- arrow:

  Logical, whether to draw arrows indicating the trajectory over time
  (default is TRUE)

- nrow:

  Number of rows for arranging the PCA plots (default is 1)

- fontsize:

  Font size for the PCA plots (default is 8)

- assay:

  Assay index to use, where 1 is the original data and 2 is normalised
  to time 0 (if available) (default is 1)

## Value

A series of PCA plots showing the distribution of samples in each group,
coloured by Time.

## Examples

``` r
data("example")
plot_pca_by_group(example_obj, circle = TRUE, arrow = TRUE)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_segment()`).
```
