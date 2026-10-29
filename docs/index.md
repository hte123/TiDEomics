# TiDEomics: Time-course Differential Expression analysis of omics data

`TiDEomics` provides an integrated framework for **Ti**me-course
**D**ifferential **E**xpression analysis of **omics** data with
**multiple experimental groups / conditions**, including integration
across datasets or omics layers.

`TiDEomics` enables:

- Analysis of feature x sample format data from different omics
  technologies.
- Comparison of ≥2 time courses systematically.
- Filtering of features by residual variance, accounting for group and
  temporal structure.
- Identification of features (e.g. genes, proteins) differentially
  expressed by **time**, **group**, and **group-by-time** effects.
- Integration of quality control, data processing, sample-level and
  feature-level analysis, tailored for time-course multi-group data.
- Functional enrichment with integrated visualisation.
- **Multi-dataset/multi-omics** comparison and consensus module
  identification using shared features.

Existing tools for time-course analysis often focus on specific analysis
tasks or experimental designs, while TiDEomics provides a comprehensive
workflow with multi-group and multi-omics support.

![](reference/figures/TiDEomics_workflow_v0.99.6.png)

Example applications include:

- Transcriptomics analysis of a cell line’s response to different
  stimuli over time.
- Transcriptomics + proteomics analysis of multiple cell lines sampled
  at several time points.

`TiDEomics` supports datasets with missing values (e.g. mass
spectrometry-based proteomics), and operates on `SummarizedExperiment`
objects to ensure compatibility with the Bioconductor ecosystem.

`TiDEomics` supports both independent samples (e.g. cell culture) and
repeated measurements (e.g. patient longitudinal studies). Subject ID
can be supplied in
[`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)
and used in subject-aware functions,
e.g. [`DE_timecourse()`](https://hte123.github.io/TiDEomics/reference/DE_timecourse.md),
[`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md).

## Installation

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `TiDEomics` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("TiDEomics")
```

The development version can be installed from
[GitHub](https://github.com/hte123/TiDEomics) with:

``` r

if (!require("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("hte123/TiDEomics")
```

## Example

For detailed examples and explanations, please refer to the [tutorials,
applications](https://hte123.github.io/TiDEomics/articles/), and other
package documentation.

## Citation

Below is the citation output from using `citation('TiDEomics')` in R.
Please run this yourself to check for any updates on how to cite
**TiDEomics**.

``` r

print(citation("TiDEomics"), bibtex = TRUE)
#> To cite package 'TiDEomics' in publications use:
#> 
#>   He T (2026). _TiDEomics: Time-course Differential Expression analysis
#>   of omics data_. R package version 0.99.6,
#>   <https://github.com/hte123/TiDEomics>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {TiDEomics: Time-course Differential Expression analysis of omics data},
#>     author = {Tianen He},
#>     year = {2026},
#>     note = {R package version 0.99.6},
#>     url = {https://github.com/hte123/TiDEomics},
#>   }
```

Please note that the `TiDEomics` was only made possible thanks to many
other R and bioinformatics software authors, which are cited in the
vignettes.

## Code of Conduct

Please note that the `TiDEomics` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Code testing is possible thanks to
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* and
  *[BiocCheck](https://bioconductor.org/packages/3.24/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://hte123.github.io/TiDEomics) is
  updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.24/biocthis)*.
