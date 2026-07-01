# Prepare for plotting WGCNA modules

Prepare the input data for plotting WGCNA modules

Here, different from running WGCNA, the input data should have the
replicates merged, instead of having multiple samples per group, time
and feature (gene).

If certain time points are missing in some groups, NA values are added.

## Usage

``` r
.plot_modules_input(module, se_obj_merged, assay, scale)
```

## Arguments

- module:

  A data frame with columns "Feature" and "Module"

- se_obj_merged:

  A SummarizedExperiment object, with one value for each feature at each
  time point in each group (replicates merged). The colData of the
  object should contain columns "Sample", "Group", and "Time". The
  object can be produced by
  [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md),
  [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md)
  and
  [`merge_groups()`](https://hte123.github.io/TiDEomics/reference/merge_groups.md).

- assay:

  The assay index in the SummarizedExperiment object to use

- scale:

  Whether to scale the data (z-score) across samples for each feature

## Value

A long-format data frame suitable for ggplot2, with columns "Feature",
"Module", "Abundance", "Group", and "Time"
