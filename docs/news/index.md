# Changelog

## TiDEomics 0.99.0

NEW FEATURES

- Added a `NEWS.md` file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- New package.
- Initial submission to Bioconductor.

BUG FIXES

- Your bug fixes. See more details at
  <http://bioconductor.org/developers/package-guidelines/#news>.

## TiDEomics 0.99.1

- Added input validation to all user-facing function arguments with
  descriptive errors.
- Added [`match.arg()`](https://rdrr.io/r/base/match.arg.html) for
  fixed-set arguments (`theme_custom$legend_position`,
  `prepare_WGCNA$networkType`).
- Added `assay` parameter +
  [`.match_assay()`](https://hte123.github.io/TiDEomics/reference/dot-match_assay.md)
  to
  [`plot_cor_matrix()`](https://hte123.github.io/TiDEomics/reference/plot_cor_matrix.md).
- Replaced `WGCNA::unsignedKME()` (unexported) with
  `abs(WGCNA::signedKME())`.
- Removed `S4Vectors::` calls; added `S4Vectors`, `SingleCellExperiment`
  to Suggests.
- Fixed roxygen inconsistencies (defaults, return types).
- Added test coverage for uncovered branches across multiple functions.
- Reordered vignette; fixed prepare_tide mapping table and assay naming
  notes.

## TiDEomics 0.99.2

NEW FEATURES

- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)
  now supports multi-category enrichment annotation via the
  `enrich_category` parameter, which accepts a vector of category names
  (e.g. `c("BP", "CC", "Hub features")`). Multiple enrichment columns
  are displayed side-by-side in a single annotation block using the new
  internal `.anno_multicol_textbox` system. The new `enrich_p_threshold`
  parameter controls per-category p-value colour-coding of enrichment
  terms.
- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)
  can now display hub features as a dedicated column in the enrichment
  textbox via `enrich_category = "Hub features"`, using the same
  `mark_features` colours as `anno_mark`. When hub features are shown as
  a textbox, the left-side `anno_mark` is automatically suppressed.
- New internal function
  [`.module_labels()`](https://hte123.github.io/TiDEomics/reference/dot-module_labels.md)
  normalises and orders module labels: numeric labels are ordered by
  value, WGCNA colour names preserve their original order. A `prefix`
  parameter controls `"M"` prefixing.
- [`extract_hubs()`](https://hte123.github.io/TiDEomics/reference/extract_hubs.md)
  extracts the top N hub features per WGCNA module using intra-modular
  connectivity (kME). Returns a character vector suitable for
  `mark_features` in
  [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md).
- `enrich_rank_by` and `enrich_top_n` in
  [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)
  now accept per-category vectors matching the length of
  `enrich_category`, enabling different ranking columns and top-N
  cutoffs for each enrichment column.
- [`as_DeeDeeExperiment()`](https://hte123.github.io/TiDEomics/reference/as_DeeDeeExperiment.md)
  converts
  [`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md)
  output to a `DeeDeeExperiment` object. DE results are flattened from
  the nested contrast-by-time format and columns are renamed to the
  `DeeDeeExperiment` convention (`logFC` → `log2FoldChange`, etc.).
  Enrichment results from
  [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)
  are automatically unwrapped.
- New
  [`.check_numeric()`](https://hte123.github.io/TiDEomics/reference/dot-check_numeric.md)
  validator for numeric arguments; `.check_*` validation added across
  all remaining exported functions with previously unchecked arguments
  (`extract_hubs`, `WGCNA_module`, `decomp_variance`, `create_input`,
  `calc_feature_property`, `split_groups`, `plot_trend`, `plot_volcano`,
  `plot_breakpoints`, `plot_DE_between_group`, `prepare_WGCNA`,
  `impute_groups`, `enrichGO_rank`, `enrichR_list`, `prepare_tide`,
  `plot_pca_3D`).
- [`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md)
  provides end-to-end data preparation: creates the
  SummarizedExperiment, normalises to starting time point, filters
  features by detection completeness (via
  [`group_specific_features()`](https://hte123.github.io/TiDEomics/reference/group_specific_features.md))
  and residual variance, and returns both unmerged (replicate-level) and
  merged (mean per Group-Time) outputs ready for downstream analysis.
- [`enrich_msigdb()`](https://hte123.github.io/TiDEomics/reference/enrich_msigdb.md):
  hypergeometric enrichment against MSigDB via `msigdbr`, with
  `species`/`db_species` support and
  [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)-compatible
  output.
- [`summarise_module_metrics()`](https://hte123.github.io/TiDEomics/reference/summarise_module_metrics.md):
  per-module kME statistics (`MeanKME`, `MeanKME2`, `MedianKME`,
  `SDKME`, `MinKME`, `MaxKME`). Auto-detects signed/unsigned networks;
  preserves kME sign.

SIGNIFICANT USER-VISIBLE CHANGES

- [`WGCNA_module()`](https://hte123.github.io/TiDEomics/reference/WGCNA_module.md)
  now returns modules ordered by decreasing module size instead of
  alphabetically, matching WGCNA convention. The Module column retains
  bare numeric labels (e.g. `"0"`, `"1"`, `"2"`) for backwards
  compatibility.
- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)
  now normalises module labels internally with
  [`.module_labels()`](https://hte123.github.io/TiDEomics/reference/dot-module_labels.md)
  and automatically excludes the grey module (`"M0"`). The
  multi-category enrichment system supersedes the previous
  single-category path.
- New internal annotation infrastructure in `anno_multicol_textbox.R`
  provides a local `.textbox_grob` (avoids `ComplexHeatmap:::`
  internals), `.reformat_cat_terms`, and `.anno_multicol_textbox` for
  building multi-column enrichment annotations.
- [`.module_labels()`](https://hte123.github.io/TiDEomics/reference/dot-module_labels.md)
  gained a `prefix` parameter to control `"M"`-prefixing; WGCNA colour
  names now preserve first-appearance order instead of being silently
  alphabetised by [`factor()`](https://rdrr.io/r/base/factor.html).
- [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md)
  now returns `ref_groups` and `all_groups` elements for cleaner
  downstream use in
  [`plot_DE_between_group()`](https://hte123.github.io/TiDEomics/reference/plot_DE_between_group.md).
- [`plot_pca()`](https://hte123.github.io/TiDEomics/reference/plot_pca.md)
  and
  [`plot_umap()`](https://hte123.github.io/TiDEomics/reference/plot_umap.md)
  no longer have a `plot` parameter; they always return a list with
  `p_list` (named list of ggplot objects) and `pca` / `umap_layout` (the
  underlying data). Use `$pca` or `$umap_layout` to access the data
  previously returned by `plot = FALSE`.

BUG FIXES

- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md):
  when `mark_features` is a plain character vector and `enrich_category`
  includes `"Hub features"`, the features now correctly appear in the
  textbox column (previously they were silently dropped).
- [`extract_hubs()`](https://hte123.github.io/TiDEomics/reference/extract_hubs.md):
  fixed a bug where a missing kME column would cause an error instead of
  skipping the module with a message.
- [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md):
  fixed `tb_name` → `t` typo causing “object not found” error; return
  element renamed to `ref_groups` (plural) for consistency.
- [`WGCNA_module()`](https://hte123.github.io/TiDEomics/reference/WGCNA_module.md):
  fixed loss of feature names (row.names) when using
  [`.module_labels()`](https://hte123.github.io/TiDEomics/reference/dot-module_labels.md),
  which caused empty downstream data.
- [`summarise_module_pattern()`](https://hte123.github.io/TiDEomics/reference/summarise_module_pattern.md):
  fixed `apply`/`table`/`data.frame` edge-case bug (error with
  single-column input) by replacing with
  [`tidyr::unite()`](https://tidyr.tidyverse.org/reference/unite.html) +
  [`dplyr::count()`](https://dplyr.tidyverse.org/reference/count.html).
  Removed `print_top_n`/`top_n` arguments.
- [`plot_volcano()`](https://hte123.github.io/TiDEomics/reference/plot_volcano.md):
  group/time arguments are now conditionally validated with
  [`missing()`](https://rdrr.io/r/base/missing.html), and
  `.check_character` applied to the active parameter set.
- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md):
  enrichment annotation blocks are now filtered to exclude modules with
  no terms across all categories (previously showed empty textboxes).
- [`.check_positive()`](https://hte123.github.io/TiDEomics/reference/dot-check_positive.md)
  now accepts vectors (previously required scalar length 1), fixing
  `prepare_WGCNA(powers = seq(1, 20))`.
- [`.match_assay()`](https://hte123.github.io/TiDEomics/reference/dot-match_assay.md)
  fixed: numeric assay indices on unnamed assays are now kept as-is
  instead of resolving to `NA` or empty string.
- `example_obj.rda` regenerated with properly named assays
  (`c("orig", "norm")`) matching the current
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)
  output.
- Data generation script (`inst/script/tutorial_input.R`) updated for
  GEOquery v2 (returns `SummarizedExperiment`; uses
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  instead of
  [`pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html)).
- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md):
  known issue — a blank plot may appear before the heatmap in RStudio
  when `enrich_list` is provided. This is an upstream ComplexHeatmap
  behaviour (physical-device measurement for `anno_link` in RStudio,
  introduced in v2.5.x). Saving directly to file or using knitr/R
  Markdown avoids this.
- [`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md)
  uses
  [`group_specific_features()`](https://hte123.github.io/TiDEomics/reference/group_specific_features.md)
  for the ID filtering step with a `min_groups` parameter (default 2),
  replacing inline dplyr logic. `min_groups` is validated to not exceed
  the number of groups.
- Non-ASCII characters removed from all R source files: right arrow
  (`->` replaces U+2192), multiplication sign (`x` replaces U+00D7), and
  approximately-equal (`~` replaces U+2248). Em dashes were already
  removed in a previous update.
- Tutorial data (`tutorial_data`) now normalised to `log2(CPM + 1)`
  (counts-per-million with pseudocount, then log2-transformed) instead
  of raw counts, providing a more realistic starting point for
  demonstration.
- Assays are now named `"orig"` (original data) and `"norm"`
  (normalised-to-start data).
  [`DE_between_time()`](https://hte123.github.io/TiDEomics/reference/DE_between_time.md),
  [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md),
  and
  [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md)
  no longer has a default assay. Both names and numeric indices are
  accepted across all functions.
- Computation separated from visualization:
  [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md)
  no longer produces an inline plot; use the new
  [`plot_DE_between_group()`](https://hte123.github.io/TiDEomics/reference/plot_DE_between_group.md)
  function instead. The `plot` and `fontsize` arguments have been
  removed.
- [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)
  now has `replicate_col` and `batch_col` arguments, allowing custom
  column names for replicate and batch information (matching the
  existing `subject_col` pattern).
- [`match.arg()`](https://rdrr.io/r/base/match.arg.html) is now used for
  fixed-choice arguments: `method` in
  [`plot_cor_matrix()`](https://hte123.github.io/TiDEomics/reference/plot_cor_matrix.md),
  `facet_by` in
  [`plot_distribution()`](https://hte123.github.io/TiDEomics/reference/plot_distribution.md),
  `confidence` in `decouple_score()`, and `category` in
  [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)
  and
  [`enrichGO_rank()`](https://hte123.github.io/TiDEomics/reference/enrichGO_rank.md).
- Input validation added across all exported functions:
  `SummarizedExperiment` type checks, p-value range checks (0–1),
  positive-integer checks, logical-flag checks, descriptive errors for
  incorrect argument types, and provenance checks (e.g.,
  [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
  validates its input was produced by
  [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md),
  [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)
  validates imputed data).
- [`DE_between_time()`](https://hte123.github.io/TiDEomics/reference/DE_between_time.md)
  and
  [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md)
  now warn if the selected assay appears to contain un-logged raw counts
  (limma requires log-transformed, normalised input).
- `magrittr` dependency fully removed. All `%>%` replaced with native
  `|>` pipe throughout the package and vignettes. Magrittr helper
  functions replaced with base R equivalents (`set_names`/`set_colnames`
  -\> [`stats::setNames`](https://rdrr.io/r/stats/setNames.html),
  `set_rownames` -\> base, `extract2` -\> `nrow`).
- Single-plot functions now
  [`return()`](https://rdrr.io/r/base/function.html) their plot objects
  instead of [`print()`](https://rdrr.io/r/base/print.html) side effects
  (`plot_trend`, `plot_cor_matrix`, `plot_cv`, `plot_distribution`,
  `plot_variance`, `plot_volcano`, `plot_breakpoints`). This enables
  programmatic composition and testing of returned plots.

Documentation

- Bioconductor installation instructions added to README and vignette.
- Data generation scripts moved from `data-raw/` to `inst/script/`, with
  new scripts for generating `example_net` and `example_res_list` used
  in `@examples`.
- `@examples` no longer include downstream function calls.
- Removed `sessioninfo` from Suggests.
- Vignette restructured with two paths: one-step
  [`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md)
  wrapper and manual step-by-step, with cross-references between them.
- `@returns` formatting standardised to plain prose across all functions
  (replacing `\describe{\item{...}}` blocks).
- Assays are now stored as matrices instead of data.frames throughout
  the pipeline, fixing iSEE compatibility.
- [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)
  now returns `unmerged_all` and `unmerged_simplified` (lists of
  `enrichResult` objects) for `DeeDeeExperiment` compatibility.
- [`.prepare_gene_list()`](https://hte123.github.io/TiDEomics/reference/dot-prepare_gene_list.md)
  drops empty factor levels after filtering, eliminating spurious “Gene
  list 0 is empty” messages.
- [`enrich_msigdb()`](https://hte123.github.io/TiDEomics/reference/enrich_msigdb.md)
  removed redundant `qvalue` column; `qvalue` dropped from Suggests.

BUG FIXES

- [`as_DeeDeeExperiment()`](https://hte123.github.io/TiDEomics/reference/as_DeeDeeExperiment.md):
  fixed `reducedDims()` not found (DDE 1.3.0 NAMESPACE bug),
  `as(se, "SCE")` coercion failure, empty de_results rejection, and
  enrichR column mapping for DDE compatibility.
- [`plot_distribution()`](https://hte123.github.io/TiDEomics/reference/plot_distribution.md),
  [`plot_missing()`](https://hte123.github.io/TiDEomics/reference/plot_missing.md),
  [`plot_ID()`](https://hte123.github.io/TiDEomics/reference/plot_ID.md),
  [`plot_cv()`](https://hte123.github.io/TiDEomics/reference/plot_cv.md):
  added [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  for matrix assay → dplyr compatibility.
- [`plot_cv()`](https://hte123.github.io/TiDEomics/reference/plot_cv.md):
  fixed `.` pronoun in native pipe (replaced with
  `rownames_to_column()`).
- [`plot_umap()`](https://hte123.github.io/TiDEomics/reference/plot_umap.md):
  fixed `p1` variable reference in `circle=TRUE` path.
- [`plot_breakpoints()`](https://hte123.github.io/TiDEomics/reference/plot_breakpoints.md):
  empty-list input now errors meaningfully.
- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md):
  Hub features Legend lacked `labels` argument.

Documentation

- Vignette: added
  [`enrich_msigdb()`](https://hte123.github.io/TiDEomics/reference/enrich_msigdb.md)
  section.
