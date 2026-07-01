# Group specific features

Identify features that are unique to selected groups (present in those
groups but not in any others) and annotate them with gene names using
the
[`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html)
function. The function also provides an option to perform Gene Ontology
(GO) enrichment analysis on the identified unique features using the
[`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)
function and visualize the results with a dot plot.

## Usage

``` r
group_specific_features(
  property_random_fc,
  groups = NULL,
  filter_ratio = 0.5,
  group_pct = 1,
  genename = TRUE,
  GO = TRUE,
  OrgDb = NULL,
  keytype = NULL,
  ...
)
```

## Arguments

- property_random_fc:

  A data frame containing the results of the
  [`calc_feature_property()`](https://hte123.github.io/TiDEomics/reference/calc_feature_property.md)
  function, which includes the feature names, group names, and the
  proportion of expressed values for each feature in each group. The
  data frame should have at least the following columns: "Feature",
  "Group", "Exp_ratio" and "Exp_threshold".

- groups:

  A character vector of group names to be compared. If NULL, all groups
  in the input data will be used (default is NULL).

- filter_ratio:

  A numeric value between 0 and 1 specifying the minimum proportion of
  expressed time points (minimum Exp_ratio) for a feature to be
  considered present (default is 0.5, meaning that a feature must have
  \>=50% values \> threshold in a group (\>= 0.5 Exp_ratio) to be
  considered present in that group).

- group_pct:

  A numeric value between 0 and 1 specifying the percentage of groups in
  which a feature must be present. (default is 1, meaning that a feature
  must be present in all specified groups).

- genename:

  A logical value indicating whether to output a table of gene names. If
  TRUE, the function will use the
  [`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html)
  function to annotate the features with gene names based on the
  specified `OrgDb` and `keytype`. (default is TRUE).

- GO:

  A logical value indicating whether to perform Gene Ontology (GO)
  enrichment analysis on the identified unique features. If TRUE, the
  function will use the
  [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)
  function to perform GO enrichment analysis and visualize the results
  with a dot plot. (default is TRUE).

- OrgDb:

  An OrgDb object from the `AnnotationDbi` package corresponding to the
  organism of interest (e.g., `org.Hs.eg.db` for human, `org.Mm.eg.db`
  for mouse). This will be used for gene annotation with the
  [`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html)
  function.

- keytype:

  A character string specifying the type of gene identifiers used in the
  row names of the assay data (e.g., "SYMBOL", "ENTREZID", "ENSEMBL").
  This will be used for gene annotation with the
  [`clusterProfiler::bitr()`](https://rdrr.io/pkg/clusterProfiler/man/bitr.html)
  function. Available key types depend on the `OrgDb` database and can
  be checked with the
  [`AnnotationDbi::keytypes`](https://rdrr.io/pkg/AnnotationDbi/man/AnnotationDb-class.html)
  function.

- ...:

  Additional arguments to be passed to the
  [`enrichGO_list()`](https://hte123.github.io/TiDEomics/reference/enrichGO_list.md)
  function for GO enrichment analysis (e.g., `pvalueCutoff`,
  `qvalueCutoff`, etc.).

## Value

A named list with element `features` (character vector of features
identified as unique to the specified groups based on the filtering
criteria). If `genename` is TRUE, element `genename` contains a
data.frame of gene annotations. If `GO` is TRUE and enrichment succeeds,
element `GO` contains a dot plot of GO enrichment results. Returns
`NULL` if no unique features pass the filter.

## Examples

``` r
data(example_obj)
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.
example_obj_list <- split_groups(example_obj)
example_obj_merged_list <- merge_replicates(example_obj_list)

example_obj_merged_list <-
    calc_feature_property(example_obj_merged_list, threshold = 0)
property_random_fc <- summarise_feature_property(example_obj_merged_list)

group_specific_features(property_random_fc, groups = c("untreated"),
    genename = FALSE, GO = FALSE)
#> Filtering criteria: >=50% values >0 in >=1 of groups: untreated
#> No unique features with the specified filter and groups.
#> NULL
```
