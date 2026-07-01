# Warn if assay data looks like raw counts rather than log-transformed values

limma expects log-transformed input. Raw counts (integers, large max)
produce unreliable results. This check warns if the data has
characteristics of un-logged counts: all non-negative, \>90% of values
are integer-like, and max value \> 100.

## Usage

``` r
.check_limma_input(mat, assay_name)
```

## Arguments

- mat:

  A numeric matrix of expression values.

- assay_name:

  Name of the assay (for the warning message).
