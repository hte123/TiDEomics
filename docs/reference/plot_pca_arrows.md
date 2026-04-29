# Plot PCA with arrows

Plot PCA with arrows indicating trajectory over time

## Usage

``` r
plot_pca_arrows(se_obj, circle = TRUE, arrow = TRUE, fontsize = 8, assay = 1)
```

## Arguments

- se_obj:

  A SummarizedExperiment object

- circle:

  Logical, whether to draw circles (ellipses) around samples of each
  time point (default is TRUE)

- arrow:

  Logical, whether to draw arrows indicating the trajectory over time
  (default is TRUE)

- fontsize:

  Font size for the PCA plot (default is 8)

- assay:

  Assay index to use, where 1 is the original data and 2 is normalised
  to time 0 (if available) (default is 1)

## Value

A PCA plot showing the distribution of samples, coloured by Time, with
arrows indicating the trajectory over time.

## Examples

``` r
data("example")
plot_pca_arrows(example_obj[, example_obj$Group == "IFNbeta"])
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_segment()`).
```
