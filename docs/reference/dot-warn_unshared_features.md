# Restrict module features to those shared by all datasets

Restrict the module table to the features present in every
SummarizedExperiment in `se_obj_list`. When non-shared features are
present, it warns with the count per dataset and errors when no feature
is shared.

## Usage

``` r
.warn_unshared_features(module, se_obj_list, assay, names = NULL)
```

## Arguments

- module:

  A data frame with columns "Feature" and "Module"

- se_obj_list:

  A list of SummarizedExperiment objects

- assay:

  The assay to use in the SummarizedExperiment objects: a numeric index
  or character name

- names:

  Optional character vector of dataset names used in warning, default:
  NULL

## Value

The subset of `module` with features present in all datasets
