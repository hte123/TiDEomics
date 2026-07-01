# Plot GO enrichment

Visualization of GO enrichment analysis with dotplot, cnetplot, or
emapplot in clusterProfiler

## Usage

``` r
plot_GO(
  go_list,
  plot_dotplot = TRUE,
  plot_cnetplot = FALSE,
  plot_emapplot = FALSE,
  showCategory_dotplot = 5,
  showCategory_cnetplot = 5,
  showCategory_emapplot = 5,
  fontsize = 8,
  label = "genes",
  ...
)
```

## Arguments

- go_list:

  A list of enriched GO results, output from
  [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)

- plot_dotplot:

  Whether to plot dotplot (default is TRUE)

- plot_cnetplot:

  Whether to plot cnetplot (default is FALSE)

- plot_emapplot:

  Whether to plot emapplot (default is FALSE)

- showCategory_dotplot:

  showCategory parameter of
  [`clusterProfiler::dotplot()`](https://rdrr.io/pkg/enrichplot/man/dotplot.html)
  (default is 5)

- showCategory_cnetplot:

  showCategory parameter of
  [`clusterProfiler::cnetplot()`](https://rdrr.io/pkg/ggtangle/man/cnetplot.html)
  (default is 5)

- showCategory_emapplot:

  showCategory parameter of
  [`clusterProfiler::emapplot()`](https://rdrr.io/pkg/enrichplot/man/emapplot.html)
  (default is 5)

- fontsize:

  Font size for the plots (default is 8)

- label:

  Title for the plots (default is "genes")

- ...:

  Additional parameters to pass to the plotting functions
  [`clusterProfiler::dotplot()`](https://rdrr.io/pkg/enrichplot/man/dotplot.html),
  [`clusterProfiler::cnetplot()`](https://rdrr.io/pkg/ggtangle/man/cnetplot.html),
  or
  [`clusterProfiler::emapplot()`](https://rdrr.io/pkg/enrichplot/man/emapplot.html)

## Value

A list of plots visualizing the GO enrichment results.

## Examples

``` r
if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(org.Mm.eg.db)
    data(example_net)
    # select two modules for demonstration
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))
    # set cutoff to 1 to show all results for demonstration
    example_go_list = enrichGO_list(example_module, OrgDb = org.Mm.eg.db,
        universe = example_module$Feature,
        pvalueCutoff = 1, qvalueCutoff = 1,
        category = "BP", simplify = FALSE)
    plot_GO(example_go_list$all, plot_dotplot = TRUE,
        plot_emapplot = FALSE, plot_cnetplot = FALSE)
}
#> Performing GO enrichment for category: BP
#> Processing gene list: 1
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> Processing gene list: 2
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> Merging GO enrichment results across gene lists for each category.
#> $dotplot_BP

#> 
```
