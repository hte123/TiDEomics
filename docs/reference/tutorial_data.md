# Dataset for TiDEomics tutorial, expression matrix

A subset of GSE263759 data set published in [Integrated time-series
analysis and high-content CRISPR screening delineate the dynamics of
macrophage immune
regulation](https://doi.org/10.1016/j.cels.2025.101346)

## Usage

``` r
data(tutorial_data)
```

## Format

A data.frame with 500 rows and 41 variables:

- Feature:

  Feature ID, e.g. gene symbols

- Sample1, Sample2, ...:

  log2(CPM + 1) normalised expression values

## Source

GSE263759

## Details

Code for preparing the data is available in
`inst/script/tutorial_input.R`

- Ensembl IDs were mapped to symbols, genes with all zero counts were
  excluded.

- Use time 0 untreated samples for other groups' time 0.

- Include "untreated", "IFNbeta", "IFNgamma", "LPS" groups.

- Sample 500 random genes.

- Normalised to log2(CPM + 1).
