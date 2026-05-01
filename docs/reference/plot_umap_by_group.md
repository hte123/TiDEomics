# UMAP by group (one object)

Plot UMAP by group, input is a SummarizedExperiment object with a
"Group" column in the colData.

## Usage

``` r
plot_umap_by_group(
  se_obj,
  seed = 1234,
  nrow = 1,
  umap_neighbors = NULL,
  fontsize = 8,
  assay = 1
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)

- seed:

  Random seed for UMAP (default is 1234)

- nrow:

  Number of rows for arranging the UMAP plots (default is 1)

- umap_neighbors:

  UMAP n_neighbors parameter (default is selected by
  `umap_n_neighbors()` function based on the number of samples)

- fontsize:

  Font size for the plot (default is 8)

- assay:

  Assay index to use, where 1 is the original data and 2 is normalised
  to time 0 (if available) (default is 1)

## Value

A series of UMAP plots showing the distribution of samples in each
group, coloured by Time.

## Examples

``` r
data("example")
plot_umap_by_group(example_obj)
#> Using n_neighbors = 5
#> Using n_neighbors = 5
#> Using n_neighbors = 5
#> Using n_neighbors = 4
```
