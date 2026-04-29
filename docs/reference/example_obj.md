# SummarizedExperiment object for runnable examples

A subset of
`data_obj <- create_input(data = tutorial_data, sample_ann = tutorial_sample_info)`
for use in runnable examples in function documentation.

## Usage

``` r
data(example)
```

## Format

A SummarizedExperiment object with assays of 100 rows and 40 columns,
colData of 40 rows and 5 columns:

- colData:

  tutorial_sample_info

- assays:

  tutorial_data first 100 rows, "Feature" column as rownames

## Source

GSE263759

## Details

Code for preparing the data is available in `data-raw/tutorial_input.R`
