# Plot correlation matrix

Plot correlation matrix between samples as a heatmap

## Usage

``` r
plot_cor_matrix(
  se_obj,
  assay = 1,
  use = "pairwise.complete.obs",
  method = c("spearman", "pearson", "kendall"),
  label_group = TRUE,
  label_time = TRUE,
  label_rep = TRUE,
  label_batch = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 8,
  cellwidth = 1,
  cellheight = 1,
  title = NULL,
  shared_only = FALSE,
  names = NULL,
  plot_cophenetic = FALSE,
  ...
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md),
  or a named list of SummarizedExperiment objects (one per dataset).

- assay:

  The assay to use in the SummarizedExperiment object: a numeric index
  or character name, e.g. 1 or "orig" for original data, 2 or "norm" for
  time 0 normalised data. (Default: 1.)

- use:

  Parameter of [`stats::cor()`](https://rdrr.io/r/stats/cor.html)
  (default is "pairwise.complete.obs")

- method:

  Parameter of [`stats::cor()`](https://rdrr.io/r/stats/cor.html)
  (default is "spearman")

- label_group:

  Whether to label Group (default is TRUE)

- label_time:

  Whether to label Time (default is TRUE)

- label_rep:

  Whether to label Replicate (default is TRUE)

- label_batch:

  Whether to label Batch (default is TRUE)

- show_rownames:

  Whether to show row names in the heatmap (default is FALSE)

- show_colnames:

  Whether to show column names in the heatmap (default is FALSE)

- fontsize:

  Font size for the plot (default is 8)

- cellwidth:

  Cell width for the heatmap (default is 1)

- cellheight:

  Cell height for the heatmap (default is 1)

- title:

  Title of the heatmap. If NULL, auto-generated as "Sample correlation
  ()" with the method capitalised (e.g. "Sample correlation (Spearman)")
  (default: NULL)

- shared_only:

  Logical, whether to subset each dataset to the features shared by all
  datasets before computing correlations (default is FALSE; only applies
  when `se_obj` is a list of multiple datasets)

- names:

  Optional character vector of dataset names for the per-heatmap titles.
  Default: NULL, i.e. the names of `se_obj` (when a list). Only applies
  when `se_obj` is a list of multiple datasets.

- plot_cophenetic:

  Logical, whether to include the pairwise Pearson correlations of
  cophenetic distances in the combined plot title (default is FALSE) in
  addition to reporting by text.

- ...:

  Additional arguments to be passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## Value

A heatmap showing the correlation between samples. With a named list of
datasets, all heatmaps share the same colour scale, including pairwise
correlations of cophenetic distances in the combined title when
`plot_cophenetic` is TRUE, and are returned invisibly as a named list of
Heatmap objects.

## Details

With a named list of datasets that share the same samples, the
per-dataset dendrograms are compared pairwise using Pearson correlations
of their cophenetic distances (see
[`?stats::cophenetic`](https://rdrr.io/r/stats/cophenetic.html)).

## Examples

``` r
data(example_obj)
data(example_obj_2)
plot_cor_matrix(example_obj)

plot_cor_matrix(list(Exp1 = example_obj, Exp2 = example_obj_2))
#> Cophenetic correlation between datasets (Pearson): 0.83
```
