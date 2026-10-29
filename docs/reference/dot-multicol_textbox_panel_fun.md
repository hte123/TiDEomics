# Create a panel function for a reusable textbox layout

This draws only the textbox panel. Heatmap connectors are drawn by the
outer annotation layout so profile and textbox share coordinates.

## Usage

``` r
.multicol_textbox_panel_fun(
  layout,
  background_gp,
  show_headers = FALSE,
  header_gp = grid::gpar(fontface = "bold"),
  header_offset = unit(1.5, "mm")
)
```

## Arguments

- layout:

  Output from `.multicol_textbox_layout`.

- background_gp:

  Background graphical parameters.

- show_headers:

  Logical; draw category headers for the first module.

- header_gp:

  Graphical parameters for category headers.

- header_offset:

  Offset above the textbox in grid units.

## Value

A panel function accepting `index` and `nm`.
