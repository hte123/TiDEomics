# Dataset for TiDEomics tutorial, sample information

A subset of GSE263759 data set published in [Integrated time-series
analysis and high-content CRISPR screening delineate the dynamics of
macrophage immune
regulation](https://doi.org/10.1016/j.cels.2025.101346)

## Usage

``` r
data(tutorial_sample_info)
```

## Format

A data.frame with 40 rows and 5 variables:

- Sample:

  Sample ID

- Group:

  Experimental group (untreated, different treatments)

- Time:

  Time point

- Replicate:

  Replicate ID for each group and time point

- Batch:

  Batch information

## Source

GSE263759

## Details

Code for preparing the data is available in `data-raw/tutorial_input.R`

- Ensembl IDs were mapped to symbols, genes with all zero counts were
  excluded.

- Use time 0 untreated samples for other groups' time 0.

- Include "untreated", "IFNbeta", "IFNgamma", "LPS" groups.

- Sample 500 random genes.
