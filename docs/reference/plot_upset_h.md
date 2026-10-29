# UpSet plot of feature sets with horizontally aligned enrichment terms

Combines an UpSet plot with row annotations for intersection size and
optional textboxes of top enriched functional terms.

## Usage

``` r
plot_upset_h(
  set_list,
  enrich_list = NULL,
  enrich_category = NULL,
  enrich_top_n = 3,
  enrich_rank_by = "p.adjust",
  enrich_threshold = 0.05,
  enrich_wrap_width = NULL,
  set_order = NULL,
  intersection_order = NULL,
  intersection_labels = NULL,
  intersection_colours = NULL,
  fontsize = 8,
  ...
)
```

## Arguments

- set_list:

  A named list of character vectors of feature sets (e.g.,
  `list(Exp1_Up = ..., Exp1_Down = ..., Exp2_Up = ..., Exp2_Down = ...)`).
  Empty input sets and zero-size intersections are excluded from the
  plot.

- enrich_list:

  A named enrichment list for intersections, such as the output of
  `enrichGO_list()$all` or
  [`enrichR_list()`](https://hte123.github.io/TiDEomics/reference/enrichR_list.md).
  Each named element should be a data.frame with `Cluster`,
  `Description`, and the column named by `enrich_rank_by`. If NULL
  (default), no enrichment terms are displayed. The `Cluster` column may
  contain binary combination names or displayed UpSet labels. Displayed
  labels are auto-generated from `names(set_list)` by joining
  participating set names with `_`, and may be overridden with
  `intersection_labels`. See
  [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)
  for the same per-category enrichment parameter convention. When
  `enrich_category` is a vector, each category is shown as a separate
  column in the annotation textbox.

- enrich_category:

  Category name(s) of enrichment to display: a single string (e.g.
  `"BP"`), or a vector for multiple categories (e.g. `c("BP", "CC")`).
  Categories must be present in `enrich_list`.

- enrich_top_n:

  Top N enrichment terms per intersection. A single integer (default: 3)
  is recycled for all categories in `enrich_category`. Supply an integer
  vector of matching length for per-category control.

- enrich_rank_by:

  Numeric column used to rank enrichment terms within each intersection.
  A single string (default: `"p.adjust"`) is recycled for all categories
  in `enrich_category`. Supply a character vector of matching length for
  per-category control. Also used for colour-coding: terms with values
  below `enrich_threshold` are shown in black, terms above in grey.

- enrich_threshold:

  Threshold for colour-coding enrichment terms (default: 0.05). Terms
  where `enrich_rank_by` column value falls below this threshold are
  drawn in black, terms above in grey. Supply a vector matching
  `enrich_category` for per-category control, or set to `NULL` to skip
  colour-coding (all terms black).

- enrich_wrap_width:

  Optional maximum enrichment term width in cm. The default `NULL`
  disables wrapping.

- set_order:

  Optional order of the sets. When `NULL` (default), sets are ordered by
  decreasing size. Pass a character vector of set names to control
  display order.

- intersection_order:

  Optional order of intersection rows. When `NULL` (default),
  intersections are ordered by decreasing combination size. Pass a
  character vector of combination names (as returned by
  [`ComplexHeatmap::comb_name()`](https://rdrr.io/pkg/ComplexHeatmap/man/comb_name.html)),
  or a numeric permutation of their indices, to specify a custom order.

- intersection_labels:

  Optional named character vector mapping `comb_name` binary patterns to
  custom intersection names, overriding the auto-generated defaults.
  E.g., `c("1010" = "Concordant_Up", "0101" = "Concordant_Down")`.
  Entries for patterns not present in the combination matrix cause an
  error.

- intersection_colours:

  Optional named colour vector for intersections. Names must match the
  binary combination patterns used by `intersection_labels` (e.g.
  `"1010"`). Colours are used for the intersection bars and combination
  marks.

- fontsize:

  Base font size (default: 8).

- ...:

  Additional arguments passed to
  [`ComplexHeatmap::UpSet()`](https://rdrr.io/pkg/ComplexHeatmap/man/UpSet.html).

## Value

A ComplexHeatmap UpSet plot, with optional enrichment term annotations.

## Details

The binary combination pattern and displayed intersection label are
printed when the function is called. These mappings can be used as keys
for `intersection_labels` and `intersection_colours`.

## Examples

``` r
# Mock feature sets with reproducible intersections
genes <- paste0("G", 1:200)
set_list <- list(Exp1_Up = genes[1:60], Exp1_Down = genes[61:100],
    Exp2_Up = genes[c(1:40, 101:115)], Exp2_Down = genes[c(61:80, 116:130)])
enrich_list <- list(BP = data.frame(
    Cluster = c("Exp1_Up_Exp2_Up", "Exp1_Down_Exp2_Down"),
    Description = c("term A", "term B"), p.adjust = c(0.001, 0.04)))
intersection_colours <- c("1010" = "#0072B5", "0101" = "#D55E00")
plot_upset_h(set_list, enrich_list = enrich_list,
    enrich_category = "BP", intersection_colours = intersection_colours)
#> Intersection pattern to label mapping:
#>   1010: Exp1_Up_Exp2_Up
#>   0101: Exp1_Down_Exp2_Down
#>   1000: Exp1_Up
#>   0100: Exp1_Down
#>   0010: Exp2_Up
#>   0001: Exp2_Down

# Custom intersection order uses the binary combination patterns.
plot_upset_h(set_list,
    intersection_order = c("1010", "1000", "0101", "0100", "0010", "0001"))
#> Intersection pattern to label mapping:
#>   1010: Exp1_Up_Exp2_Up
#>   0101: Exp1_Down_Exp2_Down
#>   1000: Exp1_Up
#>   0100: Exp1_Down
#>   0010: Exp2_Up
#>   0001: Exp2_Down
```
