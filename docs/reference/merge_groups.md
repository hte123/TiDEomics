# Merge groups into one object

Merge SummarizedExperiment object of different groups into one

## Usage

``` r
merge_groups(se_obj_list)
```

## Arguments

- se_obj_list:

  A list of SummarizedExperiment objects, such as output of
  [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md)
  and
  [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md).
  Names of the list components are the group names.

## Value

A SummarizedExperiment object containing all samples from the input
list. The 'Group' column in the colData will indicate the group of each
sample. New sample names will be prefixed with the group name.

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged <- merge_groups(example_obj_merged_list)
```
