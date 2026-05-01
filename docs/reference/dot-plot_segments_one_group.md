# Plot segmented regression (one group)

Plot segmented regression results for specified features

## Usage

``` r
.plot_segments_one_group(se_obj_imp, res, feature, group, nrow = NULL, ...)
```

## Arguments

- se_obj_imp:

  A SummarizedExperiment object with no missing values

- res:

  Result of the Trendy analysis, output of
  [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)

- feature:

  A character vector of feature names to be plotted

- group:

  A character string specifying the group name (for title purposes)

- nrow:

  Number of rows in the plot layout. If NULL and multiple features are
  provided, it will be set to half the number of features (rounded up)
  (default is NULL)

- ...:

  Additional arguments to be passed to the
  [`Trendy::plotFeature()`](https://rdrr.io/pkg/Trendy/man/plotFeature.html)

## Value

A plot of the segmented regression fits for the specified features.
