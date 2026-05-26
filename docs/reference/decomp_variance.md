# Variance decomposition

Variance decomposition by linear mixed models (LMM), to estimate the
contribution of Group and Time to the variance of each feature. Modified
from `PALMO::lmeVariance()` function.

Group and Time are always included as random effects. The formula is
extended automatically based on colData:

- Subject column present: adds `(1|Subject)` to model between-subject
  differences

- `interaction = TRUE`: adds `(1|Group:Time)` (or `(1|Subject:Time)`
  when Subject present) to capture interaction variance. Requires \>= 2
  replicates per combination. Default: FALSE.

Output always includes: Group, Time, Residual. Subject is included when
present. All values are percentages of total variance.

Allow using original data or data normalised to time point 0.

## Usage

``` r
decomp_variance(
  se_obj,
  features = NULL,
  fixed_effect_var = NULL,
  interaction = FALSE,
  assay = c(1, 2),
  core = 1
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)

- features:

  A vector of features to include. If NULL, all features.

- fixed_effect_var:

  Column names to include as fixed effects (regressed out, not appearing
  as variance components). Default is `NULL` (no fixed effects). Set to
  `"Batch"` to regress out batch effects, or provide other columns
  (e.g., `c("Sex", "Age")`).

- interaction:

  Logical. If TRUE, adds `(1|Group:Time)` (or `(1|Subject:Time)` when
  Subject present) to the model to capture interaction variance.
  Default: FALSE.

- assay:

  1 for original data, or 2 for data normalised to time 0

- core:

  Number of cores for parallel processing (default: 2)

## Value

A data frame with variance decomposition results (percentages).

## References

https://github.com/aifimmunology/PALMO/blob/main/R/lmeVariance.R

## Examples

``` r
data("example")
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.

var_decomp <- decomp_variance(example_obj, assay = 1)
#> LMM: exp ~ (1|Group) + (1|Time)  |  Output: Group, Time, Residual
plot_variance(var_decomp, rank = "Time", top_n = 20)
#> Features not specified. Plotting top 20 features ranked by Time.
```
