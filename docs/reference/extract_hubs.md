# Extract hub features from WGCNA modules

Identifies hub features within each module using module membership
(kME). For unsigned networks, absolute kME values are used to reproduce
the behaviour of
[`WGCNA::blockwiseModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseModules.html).
For signed and signed hybrid networks, signed kME values are retained.

## Usage

``` r
extract_hubs(net, top_n = 3, exclude_grey = TRUE)
```

## Arguments

- net:

  An output object from
  [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)

- top_n:

  Number of hub features to return per module. Default: 3.

- exclude_grey:

  Logical; if `TRUE`, exclude the grey/unassigned module. Default:
  `TRUE`.

## Value

A character vector containing the top `top_n` features from each module
concatenated together.

## Examples

``` r
data(example_net)
hubs <- extract_hubs(example_net, top_n = 3)
```
