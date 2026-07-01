# Build marked-feature term_list for multi-column textbox

Build a per-module term_list from `mark_features` for use as a column in
`.anno_multicol_textbox`. Accepts a named list (per-category colours
from [`ggsci::pal_jco()`](https://nanx.me/ggsci/reference/pal_jco.html))
or a character vector (all black). Colours are taken from
`mark_state$feat_col`.

## Usage

``` r
.build_marked_features_anno(
  mark_features,
  module_df,
  mark_state,
  fontsize,
  mod_levels
)
```

## Arguments

- mark_features:

  A named list of character vectors, mapping category names to feature
  IDs, or a plain character vector of feature IDs (all shown in black).

- module_df:

  A data frame with columns `Feature` and `Module`.

- mark_state:

  A list with element `feat_col`, a named character vector mapping
  feature IDs to colours.

- fontsize:

  Numeric font size for term text.

- mod_levels:

  Character vector of module levels to include.

## Value

A named list of per-module data.frames (columns `text`, `col`,
`fontsize`), or `NULL` if no marked features are found in the module
data.
