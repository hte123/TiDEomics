# Build the selected-module profile and textbox annotation

Combines aligned profile and enrichment panels inside a row `anno_zoom`.

## Usage

``` r
.make_selected_anno_zoom(
  profile_row_groups,
  profile_width_mm,
  link_width_mm,
  heatmap_height_mm = NULL,
  profile_height_mm = NULL,
  profile_panel_fun,
  text_panel_fun = NULL,
  text_row_groups = NULL,
  text_heights_mm = numeric(0),
  textbox_width_mm = 0,
  textbox_gap_mm = 2,
  link_gp = grid::gpar(fill = "grey98", col = "grey50", lwd = 1)
)
```

## Value

A ComplexHeatmap row annotation.
