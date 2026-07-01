# Gene set enrichment via MSigDB

Enrichment analysis against MSigDB gene sets via hypergeometric test.
Supports MSigDB categories, subcategories, and specified gene sets.
Requires the `msigdbr` package.

Use cases:

- Drug perturbations: `category = "C2", subcategory = "CGP"` (chemical
  and genetic perturbations)

- TF targets: `category = "C3", subcategory = "TFT:GTRD"`

- Hallmarks: `category = "H", subcategory = NULL` (no subcategory)

- GO: `category = "C5", subcategory = "BP"` (GO biological process)

Multiple testing: P-values are adjusted (Benjamini-Hochberg) per module
across the gene-set tests within that module.

## Usage

``` r
enrich_msigdb(
  gene_list,
  universe,
  category = NULL,
  subcategory = NULL,
  gene_sets = NULL,
  species = "Homo sapiens",
  db_species = "HS",
  name = NULL,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)
```

## Arguments

- gene_list:

  A named list of gene vectors, or a data.frame with `Feature` and
  `Module` columns from
  [`WGCNA_module()`](https://hte123.github.io/TiDEomics/reference/WGCNA_module.md).

- universe:

  Background genes (required).

- category:

  MSigDB category (default: `NULL`). See
  `msigdbr::msigdbr_collections(db_species = "HS"/"MM")` for available
  options. Not required when `gene_sets` is provided; in that case all
  collections are searched to locate the named gene sets.

- subcategory:

  MSigDB subcategory (default: `NULL`). Set to NULL to include all
  subcategories of `category`.

- gene_sets:

  Optional character vector of specific MSigDB gene set names to include
  (e.g., `c("HALLMARK_APOPTOSIS", "HALLMARK_P53_PATHWAY")`). When
  provided, `msigdbr` is queried across all collections to find these
  gene sets, so `category` and `subcategory` are not required. If NULL
  (default), all gene sets in the category/subcategory are used.

- species:

  Species of input genes (default: `"Homo sapiens"`). See
  [`msigdbr::msigdbr_species()`](https://igordot.github.io/msigdbr/reference/msigdbr_species.html)
  for available options.

- db_species:

  Species of MSigDB genesets (default: `"HS"`). Use "MM" for
  mouse-native gene sets. If species and db_species mismatch,
  `db_species` MSigDB will be used with ortholog mapping to `species`.

- name:

  Name for the output list element. If NULL, auto-derived from
  `category` (e.g., `"C2"` names the output element `"C2"`). If
  `category` is also NULL, it uses `"msigdb"`.

- pvalueCutoff:

  Adjusted p-value cutoff (default: 0.05).

- pAdjustMethod:

  P-value adjustment method (default: `"BH"`).

- minGSSize:

  Minimum gene set size (default: 10).

- maxGSSize:

  Maximum gene set size (default: 500).

## Value

A named list with one element: a data.frame with columns `Cluster`,
`ID`, `Description`, `GeneRatio`, `BgRatio`, `pvalue`, `p.adjust`,
`Count`, `geneID`.

## Examples

``` r
if (requireNamespace("msigdbr", quietly = TRUE)) {
    data(example_net)
    example_module <- WGCNA_module(example_net) |>
        dplyr::filter(Module %in% c("1", "2"))
    hallmark_msigdb <- enrich_msigdb(example_module, category = "MH",
        species = "Mus musculus", db_species = "MM",
        minGSSize = 1,
        pvalueCutoff = 0.9,
        universe = example_module$Feature)
}
#> Processing gene list: 1
#> Processing gene list: 2
```
