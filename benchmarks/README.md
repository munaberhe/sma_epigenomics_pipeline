# Parameter Benchmarking

These scripts were used to select DMR calling parameters. They are not part of the main pipeline and do not need to be rerun — parameters are locked in scripts/06_dmrcaller_by_chr.R.

## Scripts

- **benchmark_permutations.R** — read-count permutation null model. Shuffles readsM/readsN by random position index across 6 window sizes (100-2000bp), 20 permutations each.
- **parameter_benchmark_archie.R** — label-swap null model. Swaps ASO_VPA and ASO_CTRL conditions entirely.
- **benchmark_compare.R** — combines both null model results and produces comparison plots. Optimal binSize = 300bp by z-score criterion.

## Criterion

z = (D_obs - mu_null) / sigma_null

Both null models agreed: binSize=300bp maximises signal-to-noise on chr1.
