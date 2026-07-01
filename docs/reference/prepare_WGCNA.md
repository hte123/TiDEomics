# Prepare data and choose power for WGCNA

Prepare input data for WGCNA (format, QC), then choose the appropriate
soft thresholding power for network construction, by analysing scale
free topology with different soft thresholding powers.

## Usage

``` r
prepare_WGCNA(
  se_obj,
  assay = 2,
  networkType = "signed",
  RsquaredCut = 0.8,
  MeanConnectivity = 100,
  powers = NULL,
  fontsize = 8,
  ...
)
```

## Arguments

- se_obj:

  A SummarizedExperiment object. Data normalised to time point 0 can be
  in the second assay slot, created by
  [`normalise_to_start()`](https://hte123.github.io/TiDEomics/reference/normalise_to_start.md).
  Features can be pre-filtered, e.g. by residual variance calculated by
  [`decomp_variance()`](https://hte123.github.io/TiDEomics/reference/decomp_variance.md),
  to remove noisy features before running WGCNA.

- assay:

  Which assay slot of the SummarizedExperiment object to use for WGCNA
  input (default is 2, which is where the time 0 normalised data is
  stored by
  [`normalise_to_start()`](https://hte123.github.io/TiDEomics/reference/normalise_to_start.md))

- networkType:

  (Optional) Parameter of
  [`WGCNA::pickSoftThreshold()`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html)
  (default is "signed")

- RsquaredCut:

  (Optional) Parameter of
  [`WGCNA::pickSoftThreshold()`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html)
  (default is 0.8)

- MeanConnectivity:

  (Optional) Line of mean connectivity (default is 100)

- powers:

  (Optional) Parameter of
  [`WGCNA::pickSoftThreshold()`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html)
  (default: NULL, auto-assigned as `seq(1, 20)`)

- fontsize:

  Base font size for diagnostic plots (default: 8).

- ...:

  Additional parameters to be passed to
  [`WGCNA::pickSoftThreshold()`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html)

## Value

A list containing results of the scale-free topology fit indices with
different powers, suggested power, network type and prepared input data
used for reuse in
[`run_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_WGCNA.md)

## References

https://github.com/edo98811/WGCNA_official_documentation/

## Examples

``` r
data(example_obj)
example_obj <- normalise_to_start(example_obj)
#> Normalising to group baseline at each feature's first non-NA time point.

wgcna_input <- prepare_WGCNA(example_obj, assay = 2, powers = seq(1, 30),
    networkType = "signed", RsquaredCut = 0.8)
#> Allowing multi-threading with up to 24 threads.
#> pickSoftThreshold: will use block size 100.
#>  pickSoftThreshold: calculating connectivity for given powers...
#>    ..working on genes 1 through 100 of 100
#>    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
#> 1      1   0.8430  3.800         0.9140 52.5000   53.3000 61.000
#> 2      2   0.0646 -0.315         0.1410 30.5000   30.8000 42.500
#> 3      3   0.6650 -1.280         0.5850 19.0000   18.6000 31.600
#> 4      4   0.7580 -1.150         0.8160 12.5000   11.8000 24.500
#> 5      5   0.6820 -1.050         0.7270  8.6000    7.7300 19.500
#> 6      6   0.7190 -1.130         0.7330  6.1400    5.2300 15.900
#> 7      7   0.7580 -1.100         0.8040  4.5100    3.6800 13.100
#> 8      8   0.8680 -1.140         0.9360  3.4000    2.6200 11.000
#> 9      9   0.8240 -1.260         0.8680  2.6200    2.0100  9.350
#> 10    10   0.8460 -1.220         0.8770  2.0600    1.5000  8.000
#> 11    11   0.8670 -1.160         0.8860  1.6400    1.2300  6.900
#> 12    12   0.7820 -1.250         0.7270  1.3200    1.0000  5.990
#> 13    13   0.8750 -1.290         0.8720  1.0800    0.7600  5.240
#> 14    14   0.8790 -1.240         0.8780  0.8940    0.5780  4.610
#> 15    15   0.8480 -1.210         0.8220  0.7450    0.4460  4.080
#> 16    16   0.8010 -1.240         0.7450  0.6270    0.3560  3.620
#> 17    17   0.8350 -1.230         0.7900  0.5310    0.2740  3.240
#> 18    18   0.1200 -1.940        -0.1080  0.4530    0.2120  2.900
#> 19    19   0.1220 -1.870        -0.1070  0.3890    0.1670  2.610
#> 20    20   0.8880 -1.200         0.8590  0.3360    0.1340  2.360
#> 21    21   0.8730 -1.210         0.8370  0.2920    0.1080  2.140
#> 22    22   0.8080 -1.170         0.7550  0.2540    0.0886  1.950
#> 23    23   0.8450 -1.190         0.8030  0.2230    0.0739  1.780
#> 24    24   0.1450 -2.560        -0.0914  0.1960    0.0618  1.620
#> 25    25   0.1460 -2.540        -0.0918  0.1730    0.0512  1.490
#> 26    26   0.1510 -2.540        -0.0875  0.1540    0.0423  1.370
#> 27    27   0.0977 -1.490        -0.0975  0.1370    0.0344  1.260
#> 28    28   0.9180 -1.150         0.8950  0.1220    0.0280  1.160
#> 29    29   0.9220 -1.120         0.9020  0.1090    0.0225  1.080
#> 30    30   0.8760 -1.120         0.8400  0.0979    0.0181  0.999
wgcna_input$fitIndices
#>    Power   SFT.R.sq      slope truncated.R.sq     mean.k.   median.k.
#> 1      1 0.84286571  3.8040087     0.91423798 52.51960238 53.34184052
#> 2      2 0.06456817 -0.3146023     0.14136303 30.53003771 30.76853194
#> 3      3 0.66481255 -1.2767899     0.58472439 19.01674430 18.57511438
#> 4      4 0.75838466 -1.1499238     0.81633829 12.51048871 11.84379473
#> 5      5 0.68163955 -1.0462208     0.72716272  8.60319272  7.72652181
#> 6      6 0.71878451 -1.1308681     0.73335186  6.13592098  5.22626809
#> 7      7 0.75754120 -1.1024293     0.80426538  4.51077219  3.68053165
#> 8      8 0.86835973 -1.1446400     0.93573615  3.40121169  2.61507055
#> 9      9 0.82448245 -1.2633108     0.86803680  2.62007559  2.01016163
#> 10    10 0.84575084 -1.2179895     0.87657268  2.05548331  1.50432291
#> 11    11 0.86662570 -1.1568222     0.88604377  1.63803619  1.22596921
#> 12    12 0.78183532 -1.2522205     0.72720615  1.32325157  1.00069379
#> 13    13 0.87501515 -1.2939303     0.87217562  1.08177336  0.75980830
#> 14    14 0.87892749 -1.2364530     0.87771012  0.89371838  0.57809789
#> 15    15 0.84752610 -1.2051610     0.82190396  0.74530426  0.44613190
#> 16    16 0.80123884 -1.2403673     0.74533663  0.62677815  0.35607325
#> 17    17 0.83484572 -1.2265512     0.79043207  0.53111009  0.27379134
#> 18    18 0.11955170 -1.9402328    -0.10826258  0.45314831  0.21187911
#> 19    19 0.12163397 -1.8714360    -0.10737965  0.38906091  0.16718735
#> 20    20 0.88781004 -1.1996277     0.85914411  0.33595920  0.13438680
#> 21    21 0.87332488 -1.2077457     0.83713926  0.29163891  0.10825081
#> 22    22 0.80763163 -1.1709049     0.75520878  0.25439934  0.08860988
#> 23    23 0.84529468 -1.1899539     0.80308953  0.22291500  0.07387197
#> 24    24 0.14533913 -2.5600205    -0.09137579  0.19614317  0.06181385
#> 25    25 0.14636875 -2.5395493    -0.09184163  0.17325658  0.05118420
#> 26    26 0.15064567 -2.5419407    -0.08750052  0.15359366  0.04232804
#> 27    27 0.09768514 -1.4862196    -0.09748448  0.13662144  0.03441214
#> 28    28 0.91763379 -1.1487339     0.89520184  0.12190765  0.02803675
#> 29    29 0.92239635 -1.1225358     0.90195533  0.10909937  0.02253314
#> 30    30 0.87570167 -1.1205509     0.84048182  0.09790681  0.01807656
#>        max.k.
#> 1  61.0086938
#> 2  42.4575218
#> 3  31.5655844
#> 4  24.4616573
#> 5  19.5066573
#> 6  15.8859320
#> 7  13.1490520
#> 8  11.0267391
#> 9   9.3479355
#> 10  7.9986024
#> 11  6.8997128
#> 12  5.9947653
#> 13  5.2423238
#> 14  4.6113587
#> 15  4.0782261
#> 16  3.6246430
#> 17  3.2362938
#> 18  2.9018468
#> 19  2.6122467
#> 20  2.3601963
#> 21  2.1397703
#> 22  1.9461248
#> 23  1.7752746
#> 24  1.6239227
#> 25  1.4893261
#> 26  1.3691917
#> 27  1.2615932
#> 28  1.1649050
#> 29  1.0777494
#> 30  0.9989543
picked_power <- wgcna_input$powerEstimate
```
