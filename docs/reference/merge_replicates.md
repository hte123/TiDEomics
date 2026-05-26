# Merge replicates

Calculate mean of replicates for each feature at each time point for
each group.

When a Subject column is present, merging is done in two stages:

1.  Average replicates within each Group-Subject-Time combination.

2.  Average across subjects within each Group-Time combination. This
    gives equal weight to each subject regardless of replicate count.

Without a Subject column, all samples at the same Group-Time are
averaged together in a single step.

## Usage

``` r
merge_replicates(se_obj_list)
```

## Arguments

- se_obj_list:

  A list of SummarizedExperiment objects created by
  [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md),
  each corresponds to one group of samples.

## Value

A list of SummarizedExperiment objects containing the mean of replicates
for each feature at each time point for each group. Each object in the
list corresponds to one group of samples.

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged <- merge_groups(example_obj_merged_list)
```
