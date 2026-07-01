# Resolve assay argument

Accepts assay as a name or index, validates it against the object, and
returns the resolved value. Open to any assay present in the object.

## Usage

``` r
.match_assay(assay, se_obj)
```

## Arguments

- assay:

  A numeric index or character name.

- se_obj:

  A SummarizedExperiment object.

## Value

The resolved assay identifier (name if available, otherwise index).
