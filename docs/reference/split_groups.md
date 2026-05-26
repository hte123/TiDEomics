# Split groups

Split SummarizedExperiment object by groups

## Usage

``` r
split_groups(se_obj)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)

## Value

A list of SummarizedExperiment objects for each group in the input
object. Each object in the list corresponds to one group of samples.

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged <- merge_groups(example_obj_merged_list)
```
