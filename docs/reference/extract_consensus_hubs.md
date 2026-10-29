# Extract hub features from consensus WGCNA network

Identifies hub features within each consensus module using module
membership (kME). Consensus hubs are ranked using
[`WGCNA::consensusKME()`](https://rdrr.io/pkg/WGCNA/man/consensusKME.html)
for signed and signed-hybrid networks, and absolute kME for unsigned
networks.

## Usage

``` r
extract_consensus_hubs(consensus, top_n = 1, exclude_grey = TRUE)
```

## Arguments

- consensus:

  An output object from
  [`run_consensus_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_consensus_WGCNA.md).

- top_n:

  Number of hub features to return per module, for each dataset and
  consensus. Default: 1.

- exclude_grey:

  Logical; if `TRUE`, exclude the grey/unassigned module. Default:
  `TRUE`.

## Value

A data frame with columns `Feature`, `Module`, `Set`, `Rank`, `Score`,
and `ScoreType`. `Set` identifies the consensus hub or a
dataset-specific hub set, for example `"Consensus hub"` or `"RNA hub"`.
Duplicate features across sets are allowed.

## Requirements

The `consensus` object must be the output of
[`run_consensus_WGCNA()`](https://hte123.github.io/TiDEomics/reference/run_consensus_WGCNA.md),
which augments the raw
[`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html)
output with `multiExpr` and `parameters`.

## Examples

``` r
data(example_cons)
hubs <- extract_consensus_hubs(example_cons, top_n = 1)
utils::head(hubs)
#>   Feature Module           Set Rank     Score     ScoreType
#> 1    Helb      1 Consensus hub    1 0.8369927 consensus.kME
#> 2    Lcp1      3 Consensus hub    1 0.8820006 consensus.kME
#> 3  Epsti1      2 Consensus hub    1 0.8997157 consensus.kME
#> 4    Phax      1      Exp1 hub    1 0.9476520           kME
#> 5    Lcp1      3      Exp1 hub    1 0.8820006           kME
#> 6  Epsti1      2      Exp1 hub    1 0.9316782           kME
```
