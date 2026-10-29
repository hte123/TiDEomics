# SummarizedExperiment object for the simulated second dataset (Exp2)

Exp2 is a noise-perturbed copy of `example_obj` (Exp1): each gene's
expression receives independent Gaussian noise scaled to 0.5x its
standard deviation. The first 5 feature names have the suffix `_Exp2` to
introduce difference in features; the remaining 95 features are shared
with Exp1.

## Usage

``` r
data(example_obj_2)
```

## Format

A SummarizedExperiment object matching the structure of `example_obj`
(100 rows, 40 columns)

## Source

GSE263759

## Details

Code for producing the data is in `inst/script/generate_example_multi.R`
