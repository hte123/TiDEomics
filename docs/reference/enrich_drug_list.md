# Enrich for drug targets

Enrichment analysis for drug targets of multiple gene sets using
enrichR, allowing for using different or same background genes for each
gene set

## Usage

``` r
enrich_drug_list(
  gene_list,
  drug_dbs = c("DSigDB", "DrugMatrix"),
  universe = NULL,
  universe_list = NULL,
  pvalueCutoff = 0.05
)
```

## Arguments

- gene_list:

  A list of gene vectors, need to be gene symbols

- drug_dbs:

  A character vector of drug-related databases in enrichR to use for
  enrichment analysis. Available databases can be checked with
  [`enrichR::listEnrichrDbs()`](https://rdrr.io/pkg/enrichR/man/listEnrichrDbs.html).
  (default is c("DSigDB", "DrugMatrix"))

- universe:

  Background genes for all input gene sets, used if `universe_list` is
  not provided

- universe_list:

  Background genes for each input gene set, a list of gene vectors with
  the same names as `gene_list`, need to be gene symbols.

- pvalueCutoff:

  (Optional) Adjusted p-value cutoff for filtering enriched terms
  (default is 0.05)

## Value

A data frame of enriched drug targets for each gene set in the input
list.

## Examples

``` r
library(dplyr)
data(example_net)
example_module <- data.frame(Module = as.factor(example_net$colors)) %>%
    tibble::rownames_to_column("Feature") %>% arrange(Module)
example_module_list <- example_module %>% filter(Module != 0) %>%
    split(as.character(.$Module)) %>%
    lapply(`[[`, "Feature")
# drugs_tb = enrich_drug_list(example_module_list, drug_dbs = c("DSigDB"),
#     universe = example_module$Feature, pvalueCutoff = 0.1)
```
