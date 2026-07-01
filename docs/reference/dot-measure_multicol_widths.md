# Measure column widths for multi-column textbox layout

Computes the maximum grob width per category in millimetres, used to
align columns across modules.

## Usage

``` r
.measure_multicol_widths(text, textbox_args = list())
```

## Arguments

- text:

  A per-module list of category term data.frames (output from
  `.reformat_cat_terms`).

- textbox_args:

  Additional arguments for `.multicol_textbox_col_grob`.

## Value

A named numeric vector of column widths in mm.
