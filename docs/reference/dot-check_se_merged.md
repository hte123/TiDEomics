# Check that an SE has merged data (one value per Group x Time)

After
[`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md),
each SE should have exactly one sample per combination of Group and
Time.

## Usage

``` r
.check_se_merged(se_obj, arg = "se_obj")
```

## Arguments

- se_obj:

  A SummarizedExperiment object.

- arg:

  Name of the argument (for error message).
