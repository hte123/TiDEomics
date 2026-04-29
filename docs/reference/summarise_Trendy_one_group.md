# Summarise Trendy results (one group)

Summarise Trendy results into data frame (one group)

## Usage

``` r
summarise_Trendy_one_group(res, ...)
```

## Arguments

- res:

  Result of the Trendy analysis, output of
  [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)

- ...:

  Additional arguments to be passed to the
  [`Trendy::topTrendy`](https://rdrr.io/pkg/Trendy/man/topTrendy.html)

## Value

A data frame of summary results of Trendy, including breakpoints,
segment slopes and p-values for each fitted feature.
