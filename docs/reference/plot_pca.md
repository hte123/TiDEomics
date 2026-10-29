# Plot PCA

Plot principal component analysis (PCA), using features without missing
values

The samples are coloured by Group and sized by Time. Ellipses are drawn
for each Group.

## Usage

``` r
plot_pca(
  se_obj,
  plot_screeplot = TRUE,
  plot_loadings = TRUE,
  plot_morepc = TRUE,
  morepc = seq(1, 5),
  pc1 = 1,
  pc2 = 2,
  circle = FALSE,
  xlim_min = NULL,
  xlim_max = NULL,
  ylim_min = NULL,
  ylim_max = NULL,
  fontsize = 8,
  assay = 1,
  shared_only = FALSE,
  names = NULL
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md),
  or a named list of SummarizedExperiment objects (one per dataset)

- plot_screeplot:

  Logical, whether to plot screeplot with
  [`PCAtools::screeplot()`](https://rdrr.io/pkg/PCAtools/man/screeplot.html)
  (default is TRUE)

- plot_loadings:

  Logical, whether to plot loadings with
  [`PCAtools::plotloadings()`](https://rdrr.io/pkg/PCAtools/man/plotloadings.html)
  (default is TRUE)

- plot_morepc:

  Logical, whether to plot pairs of more PCs with
  [`PCAtools::pairsplot()`](https://rdrr.io/pkg/PCAtools/man/pairsplot.html)
  (default is TRUE)

- morepc:

  A numeric vector of PCs to plot in the pairs plot (default is 1:5)

- pc1:

  Numeric, which PC to use for the x-axis (default is 1)

- pc2:

  Numeric, which PC to use for the y-axis (default is 2)

- circle:

  Logical, whether to draw circles (ellipses) around samples of each
  group (default is FALSE)

- xlim_min:

  Minimum x-axis limit when drawing ellipses (default: NULL,
  auto-computed with 0.35x padding around the PC1 range)

- xlim_max:

  Maximum x-axis limit when drawing ellipses (default: NULL,
  auto-computed with 0.35x padding around the PC1 range)

- ylim_min:

  Minimum y-axis limit when drawing ellipses (default: NULL,
  auto-computed with 0.35x padding around the PC2 range)

- ylim_max:

  Maximum y-axis limit when drawing ellipses (default: NULL,
  auto-computed with 0.35x padding around the PC2 range)

- fontsize:

  Font size for the PCA plot (default is 8)

- assay:

  The assay to use in the SummarizedExperiment object: a numeric index
  or character name, e.g. 1 or "orig" for original data, 2 or "norm" for
  time 0 normalised data. (Default: 1.)

- shared_only:

  Logical, whether to subset each dataset to the features shared by all
  datasets before PCA (default is FALSE; only applies when `se_obj` is a
  list of multiple datasets)

- names:

  Optional character vector of dataset names for the plot titles.
  Default: NULL, i.e. the names of `se_obj` (when a list).

## Value

A list with `p_list` (a named list of ggplot objects) and `pca` (a
PCAtools output object for custom plotting). When `se_obj` is a named
list of datasets, `p_list` and `pca` are named lists (one per dataset)
and the combined PCA scatter plots are returned as `plot`.

## Examples

``` r
data(example_obj)
data(example_obj_2)
PC = plot_pca(example_obj, morepc = seq(1, 3))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the PCAtools package.
#>   Please report the issue to the authors.
#> -- variables retained:
#> Angpt1, Lpin1, Tmigd3, Pkib, Kitl, Efna2, Gprc5b, Slc13a3, Zbtb47, Tns4, Serpinf1
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the PCAtools package.
#>   Please report the issue to the authors.
#> Warning: Using size for a discrete variable is not advised.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
PC$p_list
#> $p3

#> 
#> $p4

#> 
#> $p1

#> 
#> $p5

#> 
PC_multi = plot_pca(list(Exp1 = example_obj, Exp2 = example_obj_2))
#> -- variables retained:
#> Angpt1, Lpin1, Tmigd3, Pkib, Kitl, Efna2, Gprc5b, Slc13a3, Zbtb47, Tns4, Serpinf1
#> Warning: Using size for a discrete variable is not advised.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> -- variables retained:
#> Angpt1, Lpin1, Efna2, Tmigd3_Exp2, Kitl, Gprc5b, Slc13a3, Zbtb47, Tns4, Serpinf1_Exp2, Clcn2
#> Warning: Using size for a discrete variable is not advised.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
PC_multi$plot
```
