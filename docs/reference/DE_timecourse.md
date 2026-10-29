# DE over time course

Differential expression analysis of temporal changes within and between
groups by limma, using time as a factor or a numeric variable for linear
or spline models.

## Usage

``` r
DE_timecourse(
  se_obj,
  group = NULL,
  reference_group = NULL,
  assay,
  group_mode = c("per_group", "joint_group"),
  time_mode = c("factor", "numeric", "spline"),
  spline_df = 3,
  trend = FALSE,
  pAdjustMethod = "BH"
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md).

- group:

  A character vector specifying which groups to analyse. If NULL, all
  groups in the `Group` column are used.

- reference_group:

  A single reference group for a `joint_group` model. If NULL, the first
  group in `group` or `colData(se_obj)$Group` is used and reported by a
  message.

- assay:

  The assay to use in the SummarizedExperiment object: a numeric index
  or character name, e.g. 1 or "orig" for original data, 2 or "norm" for
  time 0 normalised data. The selected assay should contain
  log-transformed, normalised values. A warning is issued if the data
  appears to be un-logged raw counts.

- group_mode:

  `"per_group"` or `"joint_group"`.

  - `"per_group"` tests each group separately.

  - `"joint_group"` fits one group-by-time model across the selected
    groups.

- time_mode:

  `"factor"`, `"numeric"`, or `"spline"`.

  - `"factor"` treats time as a categorical variable. The earliest
    available time point is the reference, and each later time point is
    compared with it. All temporal terms are tested jointly with an
    omnibus F-test. With one temporal term, this is equivalent to a
    two-sided t-test.

  - `"numeric"` treats time as continuous and tests one linear trend.

  - `"spline"` models time with a natural cubic spline basis, allowing a
    smooth non-linear trajectory. The spline coefficients are tested
    jointly with an omnibus F-test.

- spline_df:

  An integer of degrees of freedom when `time_mode = "spline"`,
  controlling the flexibility of the fitted trajectory. Larger values
  allow more complex trajectories but may increase overfitting. It must
  be at least 2 and smaller than the number of time points used by the
  model. Default: 3.

- trend:

  Logical, passed to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html). Use
  `TRUE` only when the mean-variance relationship is expected to vary
  with abundance. The default `FALSE` is usually appropriate for
  log-transformed, normalised data unless diagnostics support modelling
  a mean-variance trend.

- pAdjustMethod:

  P-value adjustment method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for the
  returned result tables. Must be one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html).
  Default: `"BH"`.

## Value

A list with `results`, `fits`, `design`, and `metadata`.

In a `per_group` model, `results` contains a table named `time_effect`.

In a `joint_group` model, `results` contains tables named
`reference_time_effect`, `global_time_effect`, `group_effect`,
`global_group_effect`, and `group_time_interaction`.

Result tables report an omnibus `F` statistic, `P.Value`, and
`adj.P.Val`. For a table testing a single model term, `F = t^2`.

Columns in `results` tables have effect-specific names. For time effect
tables,

- `time_mode = "factor"` uses `logFC_TimeFactor<time>` (log-scale
  changes relative to the reference time)

- `time_mode = "numeric"` uses `logFC_TimeSlope` (log-scale changes per
  time unit)

- `time_mode = "spline"` uses `SplineBasis_1`, `SplineBasis_2`, etc.

In `joint_group` results, group-by-time differences use names such as

- `time_mode = "factor"`: `delta_logFC_GroupX_TimeFactor<time>`

- `time_mode = "numeric"`: `delta_logFC_GroupX_TimeSlope`

- `time_mode = "spline"`: `delta_SplineBasis_GroupX_1`, etc.

The `global_group_effect` table reports group main-effect model
coefficients as `GroupX` (`X` is a non-reference group), alongside
group-by-time terms; `group_effect` table reports each non-reference
group's estimated difference from the reference group at the reference
time (`ReferenceTime`) by columns of `logFC_GroupX_at_reference_time`.

These names apply only to `results` table columns; `fits` and `design`
retain the original model coefficient names.

Features with missing values are fitted using their available
observations when the model remains estimable. Features with
insufficient observations or non-estimable designs return in the results
with `NA` for model coefficients, test statistics, and p-values.

## Details

If `colData(se_obj)` contains `Subject`, repeated measurements are
modelled with
[`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html).

Missingness patterns are defined by whether each Group-Time combination
contains at least one finite value.

The joint model is fit on the shared time points and shared features of
all selected groups. Excluded or non-estimable features are retained in
the results with `NA`.

Multiple testing correction is applied within each group for `per_group`
results and separately within each `joint_group` result table, after
combining all missingness patterns.

`ReferenceTime` is the earliest time point available in each fitting
pattern and is used for `group_effect` in all time modes and for
temporal comparisons in factor mode. It is `NA` for numeric and spline
temporal tests.

In a joint model,

- `reference_time_effect` tests whether the reference group changes over
  time

- `group_time_interaction` tests whether the other groups have different
  temporal patterns

- `global_time_effect` tests all time and group-by-time terms jointly,
  asking whether any selected group has a temporal effect

- `group_effect` tests selected group differences at the earliest
  (reference) time point used by each missingness pattern

- `global_group_effect` tests all group and group-by-time terms jointly,
  asking whether group trajectories differ anywhere in the fitted time
  course.

## Examples

``` r
data(example_obj)

per_group_result <- DE_timecourse(
    example_obj, assay = 2, group_mode = "per_group"
)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Warning: Zero sample variances detected, have been offset away from zero
#> Some Group-Time combinations contain only one sample; replicate-level variance cannot be estimated for them.
#> Warning: Zero sample variances detected, have been offset away from zero
#> Some Group-Time combinations contain only one sample; replicate-level variance cannot be estimated for them.
#> Warning: Zero sample variances detected, have been offset away from zero
joint_result <- DE_timecourse(
    example_obj, assay = 2,
    group = c("IFNbeta", "IFNgamma"),
    group_mode = "joint_group", reference_group = "IFNgamma"
)
#> Warning: Zero sample variances detected, have been offset away from zero
#> Warning: Zero sample variances detected, have been offset away from zero
utils::head(joint_result$results$global_time_effect)
#>    Feature Group               Test TimePointsTested ReferenceTime
#> 1   Scpep1  <NA> global_time_effect     0,2,4,6,8,24             0
#> 2      Dbt  <NA> global_time_effect     0,2,4,6,8,24             0
#> 3   Tmigd3  <NA> global_time_effect     0,2,4,6,8,24             0
#> 4    Clcn4  <NA> global_time_effect     0,2,4,6,8,24             0
#> 5 Serpinf1  <NA> global_time_effect     0,2,4,6,8,24             0
#> 6   Sema6b  <NA> global_time_effect     0,2,4,6,8,24             0
#>   logFC_TimeFactor2 logFC_TimeFactor4 logFC_TimeFactor6 logFC_TimeFactor8
#> 1        0.05264859         0.2547546         0.1546295        0.26477700
#> 2        0.17186536        -0.3536703        -0.2395246        0.05240711
#> 3        7.08175472         4.8269434         3.3569599        3.29526204
#> 4       -0.26992280        -0.5368547        -0.5466864       -0.57455595
#> 5        1.63730552        -2.9552536        -2.9552536        0.96910948
#> 6        0.32081854        -0.9809103        -1.5462760       -1.20164504
#>   logFC_TimeFactor24 delta_logFC_GroupIFNbeta_TimeFactor2
#> 1        0.267953712                           -0.1501064
#> 2        0.342456114                           -0.3522534
#> 3        3.123314469                           -0.6431340
#> 4       -0.002654067                           -0.1514822
#> 5       -2.955253581                           -2.5650979
#> 6       -1.432630206                           -0.9719256
#>   delta_logFC_GroupIFNbeta_TimeFactor4 delta_logFC_GroupIFNbeta_TimeFactor6
#> 1                          -0.74026758                           -0.6114845
#> 2                          -0.06397471                           -1.0366844
#> 3                           0.62430724                            2.3748026
#> 4                          -0.39087707                           -0.9386199
#> 5                           1.95406365                            1.9374434
#> 6                          -1.20320998                           -0.6729592
#>   delta_logFC_GroupIFNbeta_TimeFactor8 delta_logFC_GroupIFNbeta_TimeFactor24
#> 1                           -0.6359364                          1.437769e-01
#> 2                           -1.8181863                         -9.522763e-01
#> 3                            1.3019406                         -2.788658e+00
#> 4                           -1.0699083                         -3.637490e-01
#> 5                           -2.0669718                          2.512148e-15
#> 6                           -1.5413149                         -7.246965e-01
#>         AveExpr         F      P.Value    adj.P.Val  B
#> 1 -0.0003742715  2.569648 5.711542e-02 0.0887222039 NA
#> 2 -0.3563589895  5.037223 4.166249e-03 0.0090570640 NA
#> 3  3.6864772562  3.770741 1.420769e-02 0.0258321681 NA
#> 4 -0.5646653651 13.393412 2.844185e-05 0.0001673050 NA
#> 5 -1.1049378449  1.097039 4.292959e-01 0.5153483987 NA
#> 6 -1.2329493518 10.888941 8.905665e-05 0.0003710694 NA
```
