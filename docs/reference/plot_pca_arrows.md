# Plot PCA with arrows

Plot PCA with arrows indicating trajectory over time, for a single
group.

## Usage

``` r
plot_pca_arrows(
  se_obj,
  circle = TRUE,
  arrow = TRUE,
  pc1 = 1,
  pc2 = 2,
  fontsize = 8,
  assay = 1
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object (single group).

- circle:

  Logical, whether to draw ellipses around samples of each time point
  (default: TRUE).

- arrow:

  Logical, whether to draw arrows indicating the trajectory over time
  (default: TRUE).

- pc1:

  Principal component for the x-axis (default: 1).

- pc2:

  Principal component for the y-axis (default: 2).

- fontsize:

  Font size for the plot (default is 8)

- assay:

  The assay to use in the SummarizedExperiment object: a numeric index
  or character name, e.g. 1 or "orig" for original data, 2 or "norm" for
  time 0 normalised data. (Default: 1.)

## Value

A PCA plot coloured by Time, with optional ellipses and trajectory
arrows.

## Examples

``` r
data(example_obj)
plot_pca_arrows(example_obj[, example_obj$Group == "IFNbeta"])
#> Warning: Using size for a discrete variable is not advised.
```
