# Create object

Check format of input data and sample annotation, create a
SummarizedExperiment object

## Usage

``` r
create_input(data, sample_ann)
```

## Arguments

- data:

  A data frame with rows as features (e.g., genes, proteins) and columns
  as samples. The first column should contain feature identifiers (e.g.,
  gene symbols) and named as 'Feature'. The data should have been log
  transformed (and normalised if needed). Missing values are allowed in
  some of the downstream analyses.

- sample_ann:

  A data frame containing sample annotations with required columns:
  'Sample', 'Group', and 'Time'. 'Replicate' and 'Batch' are optional
  columns and will be set to 1 for all samples if not provided. 'Time',
  'Replicate' and 'Batch' should be numeric.

## Value

A SummarizedExperiment object containing the input data and sample
annotations.

## Examples

``` r
data <- data.frame(Feature = c("Gene1", "Gene2"),
                   Sample1 = c(1, 2), Sample2 = c(3, 4))
sample_ann <- data.frame(
    Sample = c("Sample1", "Sample2"),
    Group = c("A", "A"),
    Time = c(0, 1),
    Replicate = c(1, 1),
    Batch = c(1, 1)
)
se_obj <- create_input(data, sample_ann)
#> Converting 'Group' column to factor. Default order is alphabetical.
#> Converting 'Replicate' column to factor. Default order is numerical.
#> Converting 'Batch' column to factor. Default order is numerical.
```
