# Plot feature abundance over time

Plot trend of feature abundances / expression over time by mean and
standard deviation (SD) for each group. Accepts either a
SummarizedExperiment object (computing mean/SD internally via
[`calc_mean_sd()`](https://hte123.github.io/TiDEomics/reference/calc_mean_sd.md))
or a pre-computed table from
[`calc_mean_sd()`](https://hte123.github.io/TiDEomics/reference/calc_mean_sd.md).

## Usage

``` r
plot_trend(
  se_obj,
  assay = 1,
  groups = NULL,
  features,
  title = "Feature",
  ylab = "Abundance",
  errorbar = TRUE,
  fontsize = 8
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object, or a data.frame from
  [`calc_mean_sd()`](https://hte123.github.io/TiDEomics/reference/calc_mean_sd.md)
  with columns: Feature, Time, Group, Mean, SD.

- assay:

  Assay index when `se_obj` is a SummarizedExperiment (default: 1 =
  original, 2 = time-0 normalised).

- groups:

  Groups to be plotted, if NULL, all groups will be used (default is
  NULL)

- features:

  Features to be plotted, if NULL, an error will be raised

- title:

  Title of the plot (default is "Feature")

- ylab:

  Y axis label of the plot (default is "Abundance")

- errorbar:

  Whether to plot error bars (default is TRUE)

- fontsize:

  Font size for the plot (default is 8)

## Value

Plot of feature abundances over time by mean and SD

## Examples

``` r
data(example_obj)
plot_trend(example_obj,
    features = sample(rownames(example_obj), 4))
#> Group not specified. Plotting all groups: IFNbeta, IFNgamma, LPS, untreated
```
