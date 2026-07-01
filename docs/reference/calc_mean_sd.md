# Calculate mean and SD

Calculate mean and standard deviation for each group and time point

## Usage

``` r
calc_mean_sd(se_obj)
```

## Arguments

- se_obj:

  A SummarizedExperiment object with assays containing expression data.
  The first assay should be the original data, and the second (if
  present) should be normalized data.

## Value

A list containing two data frames: one for original values and one for
normalized values (if available). Each data frame includes columns for
Mean, SD, Group, Time, and Feature.

## Examples

``` r
data(example_obj)
table_mean_sd_list <- calc_mean_sd(example_obj)
table_mean_sd <- table_mean_sd_list$norm
```
