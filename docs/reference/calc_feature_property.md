# Calculate feature property

Calculate randomness p-value, maximum fold change along time course, and
ratio of expressed time points for each feature in each group, based on
the mean of replicates at each time point.

## Usage

``` r
calc_feature_property(se_obj_merged_list, threshold = NULL)
```

## Arguments

- se_obj_merged_list:

  A list of SummarizedExperiment objects created by
  [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md),
  each containing the mean of replicates for each feature at each time
  point for one group.

- threshold:

  A numeric value to determine whether a feature is "expressed" in
  samples. (default is NULL, meaning any non-NA value is considered
  expressed).

## Value

A list of SummarizedExperiment objects, each corresponding to one group,
with updated rowData containing the following columns: Feature, Group,
Exp_threshold, T_exp (number of time points with values \>
Exp_threshold), T_total (total number of time points), P_trend (p-value
from Bartels' test for randomness against a trend), Max_FC (maximum fold
change across time points), Max_FC_time (difference between the time
points at which maximum and minimum expression occur; positive if max
occurs after min, negative otherwise), Exp_ratio (T_exp / T_total).

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)

example_obj_merged_list <-
    calc_feature_property(example_obj_merged_list, threshold = 0)
property_random_fc <- summarise_feature_property(example_obj_merged_list)
```
