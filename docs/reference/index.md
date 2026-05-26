# Package index

## Data preparation

- [`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)
  : Create object
- [`normalise_to_start()`](https://hte123.github.io/TiDEomics/reference/normalise_to_start.md)
  : Normalise to time 0
- [`split_groups()`](https://hte123.github.io/TiDEomics/reference/split_groups.md)
  : Split groups
- [`merge_groups()`](https://hte123.github.io/TiDEomics/reference/merge_groups.md)
  : Merge groups into one object
- [`merge_replicates()`](https://hte123.github.io/TiDEomics/reference/merge_replicates.md)
  : Merge replicates
- [`impute_groups()`](https://hte123.github.io/TiDEomics/reference/impute_groups.md)
  : Impute missing values

## Quality control

- [`plot_distribution()`](https://hte123.github.io/TiDEomics/reference/plot_distribution.md)
  : Abundance distribution plot
- [`plot_missing()`](https://hte123.github.io/TiDEomics/reference/plot_missing.md)
  : Plot missing rate
- [`plot_ID()`](https://hte123.github.io/TiDEomics/reference/plot_ID.md)
  : Plot number of identified features
- [`plot_cv()`](https://hte123.github.io/TiDEomics/reference/plot_cv.md)
  : Plot coefficient of variation (CV)

## Feature properties

- [`calc_mean_sd()`](https://hte123.github.io/TiDEomics/reference/calc_mean_sd.md)
  : Calculate mean and SD
- [`plot_trend()`](https://hte123.github.io/TiDEomics/reference/plot_trend.md)
  : Plot feature abundance over time
- [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
  : Calculate feature property
- [`summarise_feature_property()`](https://hte123.github.io/TiDEomics/reference/summarise_feature_property.md)
  : Summarise feature properties
- [`group_specific_features()`](https://hte123.github.io/TiDEomics/reference/group_specific_features.md)
  : Group specific features
- [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md)
  : Variance decomposition
- [`plot_variance()`](https://hte123.github.io/TiDEomics/reference/plot_variance.md)
  : Plot variance decomposition

## Sample relationship

- [`plot_cor_matrix()`](https://hte123.github.io/TiDEomics/reference/plot_cor_matrix.md)
  : Plot correlation matrix
- [`plot_pca()`](https://hte123.github.io/TiDEomics/reference/plot_pca.md)
  : Plot PCA
- [`plot_pca_3D()`](https://hte123.github.io/TiDEomics/reference/plot_pca_3D.md)
  : Plot PCA in 3D
- [`plot_pca_arrows()`](https://hte123.github.io/TiDEomics/reference/plot_pca_arrows.md)
  : Plot PCA with arrows
- [`plot_pca_by_group()`](https://hte123.github.io/TiDEomics/reference/plot_pca_by_group.md)
  : Plot PCA by group
- [`plot_umap()`](https://hte123.github.io/TiDEomics/reference/plot_umap.md)
  : Plot UMAP
- [`plot_umap_by_group()`](https://hte123.github.io/TiDEomics/reference/plot_umap_by_group.md)
  : Plot UMAP by group

## Pairwise differential expression

- [`DE_between_group()`](https://hte123.github.io/TiDEomics/reference/DE_between_group.md)
  : DE between groups
- [`DE_between_time()`](https://hte123.github.io/TiDEomics/reference/DE_between_time.md)
  : DE between time points
- [`plot_DE_between_time()`](https://hte123.github.io/TiDEomics/reference/plot_DE_between_time.md)
  : DE number between time points
- [`plot_volcano()`](https://hte123.github.io/TiDEomics/reference/plot_volcano.md)
  : Volcano plot of DE results

## Segmented regression with Trendy

- [`run_Trendy()`](https://hte123.github.io/TiDEomics/reference/run_Trendy.md)
  : Segmented regression analysis
- [`summarise_Trendy()`](https://hte123.github.io/TiDEomics/reference/summarise_Trendy.md)
  : Summarise Trendy results
- [`extract_segment_trends()`](https://hte123.github.io/TiDEomics/reference/extract_segment_trends.md)
  : Extract feature trends
- [`plot_segments()`](https://hte123.github.io/TiDEomics/reference/plot_segments.md)
  : Plot segmented regression
- [`plot_breakpoints()`](https://hte123.github.io/TiDEomics/reference/plot_breakpoints.md)
  : Plot breakpoint distribution

## Temporal module identification with WGCNA

- [`prepare_WGCNA()`](https://hte123.github.io/TiDEomics/reference/prepare_WGCNA.md)
  : Prepare data and choose power for WGCNA
- [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)
  : Weighted gene co-expression network analysis
- [`plot_WGCNA()`](https://hte123.github.io/TiDEomics/reference/plot_WGCNA.md)
  : Plot WGCNA results
- [`WGCNA_module()`](https://hte123.github.io/TiDEomics/reference/WGCNA_module.md)
  : Convert WGCNA output to feature-module data frame
- [`plot_modules_h()`](https://hte123.github.io/TiDEomics/reference/plot_modules_h.md)
  : Plot modules (horizontal layout)
- [`plot_modules_v()`](https://hte123.github.io/TiDEomics/reference/plot_modules_v.md)
  : Plot modules (vertical layout)
- [`summarise_module_pattern()`](https://hte123.github.io/TiDEomics/reference/summarise_module_pattern.md)
  : Summarise module patterns

## Functional enrichment

- [`enrichGO_rank()`](https://hte123.github.io/TiDEomics/reference/enrichGO_rank.md)
  : GO enrichment with ranked gene list
- [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)
  : GO enrichment with gene sets
- [`enrichR_list()`](https://hte123.github.io/TiDEomics/reference/enrichR_list.md)
  : Gene-set enrichment via enrichR
- [`plot_GO()`](https://hte123.github.io/TiDEomics/reference/plot_GO.md)
  : Plot GO enrichment

## Visualization settings

- [`set_custom_palette()`](https://hte123.github.io/TiDEomics/reference/set_custom_palette.md)
  : Set custom color palette
- [`get_custom_palette()`](https://hte123.github.io/TiDEomics/reference/get_custom_palette.md)
  : Get custom color palette
- [`theme_custom()`](https://hte123.github.io/TiDEomics/reference/theme_custom.md)
  : Custom ggplot2 theme

## Data sets in tutorial and examples

- [`tutorial_data`](https://hte123.github.io/TiDEomics/reference/tutorial_data.md)
  : Dataset for TiDEomics tutorial, expression matrix

- [`tutorial_sample_info`](https://hte123.github.io/TiDEomics/reference/tutorial_sample_info.md)
  : Dataset for TiDEomics tutorial, sample information

- [`example`](https://hte123.github.io/TiDEomics/reference/example_obj.md)
  : SummarizedExperiment object for runnable examples

- [`example_res_list`](https://hte123.github.io/TiDEomics/reference/example_res_list.md)
  :

  `run_Trendy` output object for runnable examples

- [`example_net`](https://hte123.github.io/TiDEomics/reference/example_net.md)
  :

  [`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)
  output object for runnable examples
