# Impute missing values

Impute missing values for each group of samples in a list of
SummarizedExperiment objects. By default, replaces NA with the minimum
non-NA value per group / subject. A custom function can be supplied for
other strategies.

The input samples can contain replicates, or merged replicates (mean of
replicates).

## Usage

``` r
impute_groups(se_obj_list, fun = min, impute_by = c("group", "subject"))
```

## Arguments

- se_obj_list:

  A list of SummarizedExperiment objects, created by
  [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md)
  or
  [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md),
  each corresponds to one group of samples.

- fun:

  Function applied to the non-NA values to generate the replacement
  value. Default: `min`. Use `function(x) min(x) / 2` for half-minimum,
  `median` for median imputation, etc.

- impute_by:

  `"group"` (default): compute replacement value from all non-NA values
  in the group. `"subject"`: compute per-subject replacement value
  (requires Subject column). If a subject is all NA, fall back to
  group-level replacement.

## Value

A list of SummarizedExperiment objects with missing values imputed. Each
object corresponds to one group of samples.

## Examples

``` r
data(example_obj)
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)
example_obj_merged_imp_list <- impute_groups(example_obj_merged_list)
#> Group IFNbeta: no missing values.
#> Group IFNgamma: no missing values.
#> Group LPS: no missing values.
#> Group untreated: no missing values.
```
