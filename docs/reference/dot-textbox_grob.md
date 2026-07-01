# Local textbox grob

Stacks text strings vertically with per-string graphical parameters.
Adapted from `ComplexHeatmap:::textbox_grob` to avoid dependency on
unexported ComplexHeatmap internals.

## Usage

``` r
.textbox_grob(
  text,
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  just = "centre",
  gp = grid::gpar(),
  background_gp = grid::gpar(col = "transparent", fill = "transparent"),
  round_corners = FALSE,
  r = unit(0.1, "snpc"),
  line_space = unit(4, "pt"),
  text_space = unit(4, "pt"),
  max_width = unit(100, "mm"),
  padding = unit(4, "pt"),
  first_text_from = "top",
  add_new_line = FALSE,
  word_wrap = FALSE
)
```

## Arguments

- text:

  Character vector of text strings to display.

- x:

  X position as a [`grid::unit`](https://rdrr.io/r/grid/unit.html)
  (default: 0.5 npc).

- y:

  Y position as a [`grid::unit`](https://rdrr.io/r/grid/unit.html)
  (default: 0.5 npc).

- just:

  Justification of the viewport (default: `"centre"`).

- gp:

  A [`grid::gpar`](https://rdrr.io/r/grid/gpar.html) object with
  graphical parameters for text. Elements `col`, `fontsize`,
  `fontfamily`, and `fontface` are recycled to match `length(text)`.

- background_gp:

  A [`grid::gpar`](https://rdrr.io/r/grid/gpar.html) object for the
  background rectangle (default: transparent fill and border).

- round_corners:

  Logical; whether to draw rounded corners on the background rectangle
  (default: `FALSE`).

- r:

  Corner radius as a [`grid::unit`](https://rdrr.io/r/grid/unit.html)
  when `round_corners = TRUE` (default: 0.1 snpc).

- line_space:

  Spacing between lines as a
  [`grid::unit`](https://rdrr.io/r/grid/unit.html) (default: 4 pt).

- text_space:

  Spacing between inline text segments as a
  [`grid::unit`](https://rdrr.io/r/grid/unit.html) (default: 4 pt).

- max_width:

  Maximum width as a [`grid::unit`](https://rdrr.io/r/grid/unit.html)
  (default: 100 mm).

- padding:

  Padding around content as a
  [`grid::unit`](https://rdrr.io/r/grid/unit.html) vector of length 1,
  2, or 4 (default: 4 pt all sides).

- first_text_from:

  `"top"` (default) or `"bottom"`.

- add_new_line:

  Logical; force each text string to a new line (default: `FALSE`).

- word_wrap:

  Logical; enable word wrapping (default: `FALSE`).

## Value

A [`grid::gTree`](https://rdrr.io/r/grid/grid.grob.html) of class
`"textbox"`.
