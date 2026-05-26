# TiDEomics: Time-course Differential Expression analysis of omics data

`TiDEomics` is designed to streamline **Ti**me-course **D**ifferential
**E**xpression analysis of **omics** data with **multiple experimental
groups / conditions**, for example, different cell lines or different
treatments sampled at several time points.

The package’s main goals are:

- Compare multiple time courses systematically.
- Identify features (e.g. genes, proteins) and pathways differentially
  expressed by **time**, **condition**, and **both** factors (time x
  condition interactions).
- Provide utilities for quality control, data processing, sample-level
  and feature-level analysis, tailored for time-course multi-condition
  data.
- Output high-quality tables and figures to facilitate interpretation
  and reporting.

![](reference/figures/TiDEomics_workflow_v1.png)

The package supports datasets with missing values, and operates on
SummarizedExperiment objects to ensure compatibility with the
Bioconductor ecosystem.

`TiDEomics` supports both independent-sample designs (e.g. cell culture)
and repeated-measures designs (e.g. patient longitudinal studies, via
`subject_col` in
[`create_input()`](https://hte123.github.io/TiDEomics/reference/create_input.md)).
Downstream functions auto-detect the data structure.

## Installation

Install the development version from
[GitHub](https://github.com/hte123/TiDEomics) with:

``` r
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("hte123/TiDEomics")
```

## Example

For detailed examples and explanations, please refer to the
[tutorial](https://hte123.github.io/TiDEomics/articles/TiDEomics.html),
applications and other package documentation.

## Citation

Below is the citation output from using `citation('TiDEomics')` in R.
Please run this yourself to check for any updates on how to cite
**TiDEomics**.

``` r
print(citation("TiDEomics"), bibtex = TRUE)
#> To cite package 'TiDEomics' in publications use:
#> 
#>   He T (2026). _TiDEomics: Time-course Differential Expression analysis
#>   of omics data_. R package version 0.99.0,
#>   <https://github.com/hte123/TiDEomics>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {TiDEomics: Time-course Differential Expression analysis of omics data},
#>     author = {Tianen He},
#>     year = {2026},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/hte123/TiDEomics},
#>   }
```

Please note that the `TiDEomics` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `TiDEomics` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.24/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://hte123.github.io/TiDEomics) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.24/biocthis)*.
