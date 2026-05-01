# Plot coefficient of variation (CV)

Plot the distribution of coefficient of variation (CV) for each feature
with replicates. The CV is calculated as the standard deviation divided
by the mean of the abundance values.

## Usage

``` r
plot_cv(se_obj, fontsize = 8)
```

## Arguments

- se_obj:

  A SummarizedExperiment object, produced by
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)
  function, containing the abundance data and associated sample
  information.

- fontsize:

  (Optional) An integer specifying the font size for the plot (default
  is 8).

## Value

A plot showing the distribution of coefficient of variation (CV) for
each feature, grouped by Time and Group. If there are no replicates, a
message will be printed indicating that CV cannot be calculated.

## Examples

``` r
data("example")
plot_cv(example_obj)
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
```
