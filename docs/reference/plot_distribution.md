# Abundance distribution plot

This function generates density plots to visualise the distribution of
abundance values across samples. The user can choose to facet the plot
by "Group", "Time", or "Sample" to compare distributions across
different conditions or samples. If no faceting variable is specified, a
single density plot will be generated for all samples combined.

## Usage

``` r
plot_distribution(se_obj, facet_by = NULL, fontsize = 8)
```

## Arguments

- se_obj:

  A SummarizedExperiment object, produced by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)
  function, containing the abundance data and associated sample
  information.

- facet_by:

  (Optional) A character string specifying the variable to facet the
  plot by. It can be one of "Group", "Time", or "Sample". If NULL, no
  faceting will be applied and a single density plot will be generated
  for all samples (default is NULL).

- fontsize:

  (Optional) An integer specifying the font size for the plot (default
  is 8).

## Value

A plot showing the density distribution of abundance values, optionally
faceted by the specified variable.

## Examples

``` r
data("example")
plot_distribution(example_obj)

plot_distribution(example_obj, facet_by = "Group")
#> Picking joint bandwidth of 1.19
#> Picking joint bandwidth of 1.21
#> Picking joint bandwidth of 1.23
#> Picking joint bandwidth of 1.24
```
