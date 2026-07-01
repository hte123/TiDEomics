# Convert WGCNA output to feature-module data frame

Converts the module assignments from
[`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)
output into two-column data.frame with `Feature` and `Module` columns,
expected by
[`plot_modules_v()`](https://hte123.github.io/TiDEomics/reference/plot_modules_v.md),
[`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md),
[`summarise_module_pattern()`](https://hte123.github.io/TiDEomics/reference/summarise_module_pattern.md),
and all enrichment functions.

## Usage

``` r
WGCNA_module(net, exclude_grey = FALSE)
```

## Arguments

- net:

  The output of
  [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)
  (a list containing `$colors`)

- exclude_grey:

  Logical. If `TRUE`, features assigned to module `0` (grey /
  unassigned) are removed. Default: FALSE

## Value

A data.frame with columns `Feature` (character) and `Module` (factor
ordered by decreasing module size). When `exclude_grey = TRUE`,
grey/unassigned features are excluded.

## Examples

``` r
data(example_net)
module <- WGCNA_module(example_net)
```
