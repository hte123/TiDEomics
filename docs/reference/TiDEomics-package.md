# TiDEomics: Time-course Differential Expression analysis of omics data

TiDEomics provides a workflow for **multi-group time-course** omics data
analysis, analysing **time-dominant**, **group-dominant**, and
**group-specific temporal** effects through pairwise differential
expression, variance decomposition, and co-expression module analysis
(WGCNA). A **residual variance** filtering strategy prioritises features
with structured differential expression. It supports datasets with
missing values (e.g. mass spectrometry-based proteomics) and operates on
`SummarizedExperiment` objects for compatibility with the Bioconductor
ecosystem.

## Details

Key functionality:

- **Data preparation & normalisation:**
  [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md),
  [`prepare_tide()`](https://hte123.github.io/TiDEomics/reference/prepare_tide.md),
  [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md),
  [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md),
  [`merge_groups()`](https://hte123.github.io/TiDEomics/reference/merge_groups.md),
  [`normalise_to_start()`](https://hte123.github.io/TiDEomics/reference/normalise_to_start.md),
  [`impute_groups()`](https://hte123.github.io/TiDEomics/reference/impute_groups.md)

- **Quality control & exploration:**
  [`plot_missing()`](https://hte123.github.io/TiDEomics/reference/plot_missing.md),
  [`plot_distribution()`](https://hte123.github.io/TiDEomics/reference/plot_distribution.md),
  [`plot_ID()`](https://hte123.github.io/TiDEomics/reference/plot_ID.md),
  [`plot_cv()`](https://hte123.github.io/TiDEomics/reference/plot_cv.md),
  [`plot_cor_matrix()`](https://hte123.github.io/TiDEomics/reference/plot_cor_matrix.md),
  [`plot_pca()`](https://hte123.github.io/TiDEomics/reference/plot_pca.md),
  [`plot_pca_3D()`](https://hte123.github.io/TiDEomics/reference/plot_pca_3D.md),
  [`plot_pca_arrows()`](https://hte123.github.io/TiDEomics/reference/plot_pca_arrows.md),
  [`plot_pca_by_group()`](https://hte123.github.io/TiDEomics/reference/plot_pca_by_group.md),
  [`plot_umap()`](https://hte123.github.io/TiDEomics/reference/plot_umap.md),
  [`plot_umap_by_group()`](https://hte123.github.io/TiDEomics/reference/plot_umap_by_group.md)

- **Feature properties & variance decomposition:**
  [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md),
  [`summarise_feature_property()`](https://hte123.github.io/TiDEomics/reference/summarise_feature_property.md),
  [`group_specific_features()`](https://hte123.github.io/TiDEomics/reference/group_specific_features.md),
  [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md),
  [`plot_variance()`](https://hte123.github.io/TiDEomics/reference/plot_variance.md),
  [`plot_trend()`](https://hte123.github.io/TiDEomics/reference/plot_trend.md)

- **Pairwise differential expression:**
  [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md),
  [`DE_between_time()`](https://hte123.github.io/TiDEomics/reference/DE_between_time.md),
  [`plot_DE_between_group()`](https://hte123.github.io/TiDEomics/reference/plot_DE_between_group.md),
  [`plot_DE_between_time()`](https://hte123.github.io/TiDEomics/reference/plot_DE_between_time.md),
  [`plot_volcano()`](https://hte123.github.io/TiDEomics/reference/plot_volcano.md)

- **Segmentation regression with Trendy:**
  [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md),
  [`summarise_Trendy()`](https://hte123.github.io/TiDEomics/reference/summarise_Trendy.md),
  [`plot_segments()`](https://hte123.github.io/TiDEomics/reference/plot_segments.md),
  [`plot_breakpoints()`](https://hte123.github.io/TiDEomics/reference/plot_breakpoints.md)

- **Co-expression module identification with WGCNA:**
  [`prepare_WGCNA()`](https://hte123.github.io/TiDEomics/reference/prepare_WGCNA.md),
  [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md),
  [`WGCNA_module()`](https://hte123.github.io/TiDEomics/reference/WGCNA_module.md),
  [`extract_hubs()`](https://hte123.github.io/TiDEomics/reference/extract_hubs.md),
  [`plot_WGCNA()`](https://hte123.github.io/TiDEomics/reference/plot_WGCNA.md),
  [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md),
  [`plot_modules_v()`](https://hte123.github.io/TiDEomics/reference/plot_modules_v.md),
  [`summarise_module_pattern()`](https://hte123.github.io/TiDEomics/reference/summarise_module_pattern.md),
  [`summarise_module_metrics()`](https://hte123.github.io/TiDEomics/reference/summarise_module_metrics.md)

- **Functional enrichment:**
  [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md),
  [`enrichGO_rank()`](https://hte123.github.io/TiDEomics/reference/enrichGO_rank.md),
  [`enrichR_list()`](https://hte123.github.io/TiDEomics/reference/enrichR_list.md),
  [`enrich_msigdb()`](https://hte123.github.io/TiDEomics/reference/enrich_msigdb.md),
  [`plot_GO()`](https://hte123.github.io/TiDEomics/reference/plot_GO.md)

- **Interoperability & settings:**
  [`flatten_DE()`](https://hte123.github.io/TiDEomics/reference/flatten_DE.md),
  [`flatten_enrich()`](https://hte123.github.io/TiDEomics/reference/flatten_enrich.md),
  [`set_custom_palette()`](https://hte123.github.io/TiDEomics/reference/set_custom_palette.md),
  [`get_custom_palette()`](https://hte123.github.io/TiDEomics/reference/get_custom_palette.md),
  [`theme_custom()`](https://hte123.github.io/TiDEomics/reference/theme_custom.md)

Two experimental designs are auto-detected from the sample annotation:
independent samples and repeated measures.

See
[`vignette("TiDEomics")`](https://hte123.github.io/TiDEomics/articles/TiDEomics.md)
for a step-by-step tutorial.

## See also

Useful links:

- <https://hte123.github.io/TiDEomics>

- <https://github.com/hte123/TiDEomics>

- Report bugs at <https://github.com/hte123/TiDEomics/issues>

## Author

**Maintainer**: Tianen He <tianen.he@ndm.ox.ac.uk>
([ORCID](https://orcid.org/0000-0001-6864-0723))

Authors:

- Tianen He <tianen.he@ndm.ox.ac.uk>
  ([ORCID](https://orcid.org/0000-0001-6864-0723))
