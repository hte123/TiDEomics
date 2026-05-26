# Summarise feature properties

This function takes a list of merged `SummarizedExperiment` objects,
which are the output of the
[`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
function, and summarizes the feature properties across all groups. It
extracts the row data from each `SummarizedExperiment` object, combines
them into a single data frame, and returns this summary data frame. Each
row in the resulting data frame represents a feature, and the columns
contain the properties of that feature for each group, including the
feature name, group name, proportion of NA values for that feature in
that group, randomness p-value, and maximum fold change over time.

## Usage

``` r
summarise_feature_property(se_obj_merged_list)
```

## Arguments

- se_obj_merged_list:

  A list of merged `SummarizedExperiment` objects, output of
  [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
  function

## Value

A data frame summarizing the feature properties across all groups, with
each row representing a feature and columns containing the properties
including the feature name, group name, and the proportion of NA values
for that feature in that group, randomness p-value, and maximum fold
change over time.

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
