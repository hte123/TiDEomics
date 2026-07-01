# Multi-column textbox grob

Arranges per-category textbox grobs side-by-side in absolute mm units.
Returns a `gTree` of class `"multicol_textbox"` with attributes
`column_labels`, `column_widths_mm`, `column_gap_mm`, and `padding_mm`
for downstream header decoration.

## Usage

``` r
.multicol_textbox_grob(
  x,
  col_gap = unit(2, "mm"),
  background_gp = grid::gpar(fill = "#F7F7F7", col = "#CCCCCC"),
  padding = unit(c(0, 0, 0, 0), "mm"),
  column_widths_mm = NULL,
  ...
)
```

## Arguments

- x:

  A named list of category term data.frames (one per column).

- col_gap:

  Gap between columns as a
  [`grid::unit`](https://rdrr.io/r/grid/unit.html) (default: 2 mm).

- background_gp:

  Background graphical parameters.

- padding:

  Padding around content as a
  [`grid::unit`](https://rdrr.io/r/grid/unit.html) vector (top, left,
  bottom, right). Default: 0 mm all sides.

- column_widths_mm:

  Optional numeric vector of pre-computed column widths in mm. If NULL,
  auto-computed from grob widths.

- ...:

  Passed to `.multicol_textbox_col_grob`.

## Value

A `gTree` of class `"multicol_textbox"`.
