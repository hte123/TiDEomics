# Prepare TiDEomics input

One-step data preparation wrapper for TiDEomics: creates the
SummarizedExperiment object, normalises to starting time point, filters
out group-specific features, runs variance decomposition, and filters
out noisy features with high residual variance. Returns both unmerged
and merged (mean of replicates) outputs, with all-feature and filtered
versions of each, ready for downstream analysis.

## Usage

``` r
prepare_tide(
  data,
  sample_ann,
  subject_col = NULL,
  replicate_col = NULL,
  batch_col = NULL,
  feature_threshold = 0,
  filter_ratio = 0.5,
  min_groups = 2,
  keep = c("below_quantile", "top_n", "threshold"),
  residual_threshold = NULL,
  n_keep = 5000,
  interaction = FALSE,
  n_cores = 1
)
```

## Arguments

- data:

  A data frame with rows as features (e.g., genes, proteins) and columns
  as samples. The first column should contain feature identifiers (e.g.,
  gene symbols) and be named 'Feature'. The data should have been
  normalised and log transformed. Missing values are allowed.

- sample_ann:

  A data frame containing sample annotations with required columns:
  'Sample', 'Group', and 'Time'. 'Replicate' and 'Batch' are optional
  columns and will be auto-generated if not provided. 'Time',
  'Replicate' and 'Batch' should be numeric.

- subject_col:

  Optional: name of a column in `sample_ann` identifying biological
  subjects measured repeatedly across time points (e.g., `"PatientID"`).
  If provided, the column is renamed to 'Subject' and used for
  repeated-measures analyses downstream. If NULL (default), all samples
  are treated as independent, e.g. for cell culture experiments.

- replicate_col:

  Optional: name of a column in `sample_ann` identifying replicate IDs.
  If provided, the column is renamed to 'Replicate'. If NULL (default),
  the function looks for a column named 'Replicate'; if absent,
  replicate IDs are auto-generated within each Group and Time (and
  Subject, if provided).

- batch_col:

  Optional: name of a column in `sample_ann` identifying batch
  information. If provided, the column is renamed to 'Batch'. If NULL
  (default), the function looks for a column named 'Batch'; if absent,
  all samples are assigned to batch 1.

- feature_threshold:

  Threshold for
  [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md).
  Feature values \> feature_threshold are considered expressed (default:
  0).

- filter_ratio:

  Minimum ratio of time points a feature must be expressed in, to be
  considered expressed in a group (default: 0.5, meaning that a feature
  must be \> feature_threshold in \>=50% time points in a group to be
  considered expressed). Passed to
  [`group_specific_features()`](https://hte123.github.io/TiDEomics/reference/group_specific_features.md).

- min_groups:

  Minimum number of groups a feature must be expressed in, with the
  above filter_ratio, to pass the cross group filter (default: 2).
  Passed to
  [`group_specific_features()`](https://hte123.github.io/TiDEomics/reference/group_specific_features.md).

- keep:

  How to select features after variance decomposition:
  `"below_quantile"` (default): keep features with residual variance
  below `residual_threshold` quantile (0-1 scale). `"top_n"`: keep the
  top `n_keep` lowest-residual features. `"threshold"`: keep features
  with residual percentage below `residual_threshold` (0-100 scale; use
  100 to skip filtering).

- residual_threshold:

  Numeric: when `keep = "below_quantile"`, the quantile of residual
  variance to use as cutoff (0-1 scale). When `keep = "threshold"`, the
  percentage cutoff (0-100 scale, 100 = keep all features). (default:
  NULL, auto-assigned as 0.75 for "below_quantile", 100 for "threshold")

- n_keep:

  Number of features to keep when `keep = "top_n"` (default: 5000).

- interaction:

  Passed to
  [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md):
  if TRUE, adds `(1|Group:Time)` (or `(1|Subject:Time)` when Subject
  present) interaction term (default: FALSE).

- n_cores:

  Number of cores for variance decomposition (default: 1).

## Value

A named list with elements. The elements suffixed with "\_filtered" have
passed both the cross-group and residual filters, and are split into
\$orig and \$norm sub-lists, with variance decomposition and residual
filtering done per assay.

1.  `se` (SE with replicates, all features, both assays);

2.  `se_filtered` (a list of two SEs with replicates, filtered features
    by `orig` / `norm`, use for WGCNA);

3.  `merged_list` (per-group merged SEs, all features, use for Trendy);

4.  `merged_list_filtered` (a list of two per-group merged SEs, filtered
    features by `orig` / `norm`);

5.  `merged_se` (single merged SE, all groups combined, all features,
    use for plot_modules_v/h);

6.  `merged_se_filtered` (a list of two single merged SEs, filtered
    features by `orig` / `norm`, use for plot_modules_v/h);

7.  `variance` (a list of two variance decomposition data.frames,
    computed from `orig` / `norm` assay);

8.  `feature_property` (feature property summary, all features);

9.  `filter_summary` (per-stage filtering statistics);

10. `group_specific_filter` (a vector of group-specific features,
    removed by the cross-group filter); 11-13. `DE`, `enrichment`,
    `WGCNA` (empty on initialization).

## Details

The pipeline:

1.  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md):
    validate and build SE

2.  [`normalise_to_start()`](https://hte123.github.io/TiDEomics/reference/normalise_to_start.md):
    add time-0-normalised assay

3.  [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md)
    -\>
    [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md)
    -\>
    [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
    -\>
    [`summarise_feature_property()`](https://hte123.github.io/TiDEomics/reference/summarise_feature_property.md):
    per-feature expression statistics

4.  [`group_specific_features()`](https://hte123.github.io/TiDEomics/reference/group_specific_features.md):
    filter out group-specific features

5.  [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md)
    : LMM variance decomposition on both `"orig"` and `"norm"` assays

6.  Residual filter: filter out features with high residual variance,
    calculated from
    [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md)
    on each assay

7.  Produce filtered-unmerged and filtered-merged SEs for each assay

## Examples

``` r
data(tutorial_data)
data(tutorial_sample_info)
tide <- prepare_tide(tutorial_data, tutorial_sample_info,
    keep = "threshold", residual_threshold = 100)
#> No Subject column specified. Samples treated as independent. For repeated-measures designs, set subject_col to the column identifying biological subjects.
#> Converting 'Group' column to factor. Default order is alphabetical.
#> Converting 'Replicate' column to factor. Default order is numerical.
#> Converting 'Batch' column to factor. Default order is numerical.
#> Normalising to group baseline at each feature's first non-NA time point.
#> Preparing TiDEomics input: 500 features, 40 samples, 4 groups
#> Filtering criteria: >=50% values >0 in >=2 of groups: IFNbeta, IFNgamma, LPS, untreated
#> Cross-group filter (Exp_ratio >= 0.5 in >= 2 groups): kept 340 of 500 features (68%)
#> --- assay: orig ---
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
#> Residual filter (threshold): kept 334 of 340 features (98.2%)
#> --- assay: norm ---
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
#> Residual filter (threshold): kept 333 of 340 features (97.9%)
tide$filter_summary
#>                stage assay n_features pct_kept                       threshold
#> 1              input  <NA>        500    100.0                            <NA>
#> 2 cross_group_filter  <NA>        340     68.0 Exp_ratio >= 0.5 in >= 2 groups
#> 3    residual_filter  orig        334     66.8       skipped (threshold = 100)
#> 4    residual_filter  norm        333     66.6       skipped (threshold = 100)
```
