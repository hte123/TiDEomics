# Reformat enrichment category terms into per-module text lists

Transposes `cat_terms` from category-first (`cat_terms[[cate]][[mod]]`)
to module-first (`text[[mod]][[cate]]`) layout expected by
`.multicol_textbox_grob`. Missing categories get empty data.frames.

## Usage

``` r
.reformat_cat_terms(cat_terms)
```

## Arguments

- cat_terms:

  A named list of categories, each containing a named list of modules
  with term annotation data.frames.

## Value

A named list of modules, each containing a named list of categories with
term data.frames.
