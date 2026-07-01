# Build a single-column textbox grob

Stacks term text strings vertically using the local `.textbox_grob`.
Returns an invisible placeholder grob if `df` has no rows.

## Usage

``` r
.multicol_textbox_col_grob(df, ...)
```

## Arguments

- df:

  A data.frame with columns `text`, `col`, `fontsize`, and optionally
  `fontfamily` and `fontface`.

- ...:

  Passed to `.textbox_grob`.

## Value

A textbox grob.
