# Impute missing values

Impute missing values for each group of samples in a list of
SummarizedExperiment objects, with the minimum value in the group.

The input samples can contain replicates, or merged replicates (mean of
replicates).

## Usage

``` r
impute_groups(se_obj_list)
```

## Arguments

- se_obj_list:

  A list of SummarizedExperiment objects, created by
  [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md)
  or
  [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md),
  each corresponds to one group of samples.

## Value

A list of SummarizedExperiment objects with the missing values imputed
for each group. Each object in the list corresponds to one group of
samples.

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
```
