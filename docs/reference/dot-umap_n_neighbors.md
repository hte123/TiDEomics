# UMAP neighbors number

Select UMAP n_neighbors parameter based on the number of samples For
sample number \> 5, use sample number / 5, with a min of 5 and max of
15. For sample number \<= 5, use sample number - 1.

## Usage

``` r
.umap_n_neighbors(sample_n)
```

## Arguments

- sample_n:

  Number of samples

## Value

UMAP n_neighbors parameter
