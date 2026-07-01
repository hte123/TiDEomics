# Check that an SE list has merged data (one value per Group x Time)

After
[`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md),
each SE should have exactly one sample per combination of Group and
Time.

## Usage

``` r
.check_se_list_merged(se_obj_list, arg = "se_obj_list")
```

## Arguments

- se_obj_list:

  A list of SummarizedExperiment objects.

- arg:

  Name of the argument (for error message).
