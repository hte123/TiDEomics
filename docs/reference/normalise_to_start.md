# Normalise to time 0

Normalise to starting time point, to make the mean of starting time
point replicates 0

## Usage

``` r
normalise_to_start(se_obj)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)

## Value

A SummarizedExperiment object with normalised data in the second assay
slot.

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged <- merge_groups(example_obj_merged_list)
```
