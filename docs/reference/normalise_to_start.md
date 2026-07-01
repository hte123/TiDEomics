# Normalise to starting time point

Normalise to starting time point, to make the mean of starting time
point samples 0.

Two modes are available:

- `by_subject = FALSE` (default): computes the group-level mean at the
  first non-NA time point per feature and subtracts it. All samples
  within a group share the same baseline.

- `by_subject = TRUE`: computes each subject's value at the first non-NA
  time point per feature and subtracts that subject-specific baseline.
  Removes between-subject baseline differences. Use when each subject
  has their own initial condition. Requires a Subject column in colData.

## Usage

``` r
normalise_to_start(se_obj, by_subject = FALSE)
```

## Arguments

- se_obj:

  A SummarizedExperiment object created by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)

- by_subject:

  If FALSE (default), use group-level baseline. If TRUE, use
  subject-level baseline (requires Subject column in colData).

## Value

A SummarizedExperiment object with normalised data in the second assay
slot.

## Examples

``` r
data(example_obj)
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
```
