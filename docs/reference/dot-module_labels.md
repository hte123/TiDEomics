# Normalise and order module labels

Converts numeric or digit-string module labels to a factor with levels
in natural numeric order. When `prefix = TRUE` (default), labels are
formatted as `"M0", "M1", "M2", ...`. When `prefix = FALSE`, bare digit
strings are returned (e.g. `"0", "1", "2"`). WGCNA colour names (e.g.
`"turquoise"`, `"blue"`) are unaffected by `prefix` and returned as a
factor preserving their original order of first appearance.

## Usage

``` r
.module_labels(x, prefix = TRUE)
```

## Arguments

- x:

  A vector of module labels (numeric, character digits, or WGCNA colour
  names).

- prefix:

  Logical; if `TRUE` (default), numeric labels are prefixed with `"M"`.
  Ignored for WGCNA colour names.

## Value

A factor of module labels. Numeric labels are ordered by their numeric
value; non-numeric labels preserve their original order of first
appearance.

## Details

Callers that need a different ordering (e.g. by module size) should
re-factor the result after calling this function, as
[`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)
does.
