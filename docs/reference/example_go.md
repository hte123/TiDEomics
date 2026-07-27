# `enrichGO_list()` output object for runnable examples

Code for producing the data is in `inst/script/generate_example_go.R`

## Usage

``` r
data(example_go)
```

## Format

A nested list with GO enrichment results for modules "1" and "2" from
'example_net', containing categories "BP" and "CC", with simplify = TRUE
applied:

- all:

  Merged results across gene sets for each GO category

- simplified:

  Simplified merged results

- unmerged_all:

  Per-gene-set results for each category

- unmerged_simplified:

  Simplified per-gene-set results

## Source

GSE263759
