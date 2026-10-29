# Plot WGCNA results

Plot WGCNA module eigengenes and module-group/time correlations based on
output of
[`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)

## Usage

``` r
plot_WGCNA(net, fontsize = 8)
```

## Arguments

- net:

  WGCNA network object output by
  [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)

- fontsize:

  Font size for plots (default is 8)

## Value

Plots of WGCNA module dendrogram, module eigengenes, pairwise
scatterplots of eigengenes, clustering of module eigengenes, and
module-trait correlation heatmap

## References

https://github.com/edo98811/WGCNA_official_documentation/blob/main/FemaleLiver-03-relateModsToExt.R

## Examples

``` r
data("example_net")
plot_WGCNA(example_net, fontsize = 8)

#> Warning: argument 1 does not name a graphical parameter
#> Warning: argument 1 does not name a graphical parameter
#> Warning: argument 1 does not name a graphical parameter
#> Warning: argument 1 does not name a graphical parameter
#> Warning: argument 1 does not name a graphical parameter
#> Warning: argument 1 does not name a graphical parameter



```
