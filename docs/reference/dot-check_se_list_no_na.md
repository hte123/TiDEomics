# Check that an SE list has no missing values (imputed)

After
[`impute_groups()`](https://hte123.github.io/TiDEomics/reference/impute_groups.md),
no assay should contain NA values.

## Usage

``` r
.check_se_list_no_na(se_obj_list, arg = "se_obj_list")
```

## Arguments

- se_obj_list:

  A list of SummarizedExperiment objects.

- arg:

  Name of the argument (for error message).
