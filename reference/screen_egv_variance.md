# Screen environmental predictors for insufficient variance

Calculates predictor standard deviations separately among presences,
background records, and their combination, including spatial CV training
partitions and (optionally) an independent test set.

## Usage

``` r
screen_egv_variance(
  train,
  test = NULL,
  occ_folds = NULL,
  bg_folds = NULL,
  env = NULL,
  env_sample_n = 10000L,
  seed = NULL,
  sd_tol = 0,
  check_presences = TRUE,
  check_background = TRUE,
  check_combined = TRUE,
  presence_affects_selection = TRUE,
  include_test_bg = FALSE,
  check_test = FALSE,
  test_affects_selection = FALSE,
  verbose = TRUE
)
```

## Arguments

- train:

  An `SDMtune` `SWD` object with presence and background data.

- test:

  An optional independent-test `SWD` object.

- occ_folds:

  An optional fold identifiers for TRAINING presences, in the order of
  `train@pa == 1`. A held-out fold is excluded in each CV check.

- bg_folds:

  An optional fold identifiers for TRAINING background records, in the
  order of `train@pa == 0`. If omitted, all training backgrounds are
  reused in every CV partition (and the diagnostics state this
  explicitly).

- env:

  An optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  representing the environmental domain.

- env_sample_n:

  Whole integer \>= 2 for random raster-cell sampling.

- seed:

  NULL or one non-negative integer for reproducible `env` sampling;
  existing R random-number-generator state is restored.

- sd_tol:

  Finite non-negative standard-deviation cutoff. Values \<= `sd_tol`, or
  fewer than two finite observations, fail a check.

- check_presences:

  Check presence-only subsets, default TRUE.

- check_background:

  Check background-only subsets, default TRUE.

- check_combined:

  Check joined presence/background subsets, default TRUE.

- presence_affects_selection:

  If FALSE, presence-only failures are diagnostic-only. Defaults to TRUE
  for compatibility with previous screening.

- include_test_bg:

  Whether to additionally check training and independent test
  backgrounds pooled together. Defaults to FALSE. The pooled check is
  diagnostic-only unless `test_affects_selection = TRUE`.

- check_test:

  Whether to diagnose the independent test set separately
  (presences/background/combined according to the `check_*` switches).
  Defaults to FALSE. Its checks are diagnostic-only by default.

- test_affects_selection:

  Explicitly allow checks involving independent test records to
  influence `keep`/`remove`. Defaults to FALSE to avoid test-informed
  preprocessing. Setting TRUE forfeits full test independence.

- verbose:

  Print a summary (TRUE/FALSE).

## Value

A named list containing:

- keep:

  Predictors retained after variance screening.

- remove:

  Predictors excluded by selection-relevant variance checks.

- removed_by_check:

  Predictors removed by each selection-relevant check.

- diagnostics:

  A data frame containing results of all performed variance checks,
  including `failed` and `used_for_selection` logical columns. The
  latter identifies whether each check contributes to predictor
  selection.

## Details

Checks performed on training observations use only `train@data`.
Supplying `env` adds an independent random environmental-domain check;
it does not subsample or modify the supplied background observations.

CV checks exclude the held-out fold's presences and, when `bg_folds` is
supplied, its background records. Without `bg_folds`, the full training
background is shared across CV partitions. CV checks never inspect
`test`.

When `check_presences = FALSE`, presence-only checks are omitted
entirely. When `presence_affects_selection = FALSE`, presence-only
checks are performed but do not influence predictor selection.

Each row in the returned `diagnostics` data frame contains a logical
column named `used_for_selection`. This column indicates whether the
corresponding variance check contributes to predictor selection. A
predictor is included in `remove` if it fails at least one check for
which `used_for_selection` is `TRUE`. Checks marked `FALSE` are retained
for diagnostic purposes but do not influence `keep` or `remove`.

`test_affects_selection = TRUE` is an explicit opt-in and is generally
inappropriate if `test` is meant as a final independent evaluation set.
Select predictors *within training data for each CV fold* if estimating
the predictive performance of an entire selection/model-fitting
pipeline.

## See also

[`SDMtune::prepareSWD()`](https://consbiol-unibern.github.io/SDMtune/reference/prepareSWD.html),
[`terra::spatSample()`](https://rspatial.github.io/terra/reference/sample.html)

## Examples

``` r
# Self-contained example with package dependencies, no external files.
r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10,
                 ymin = 0, ymax = 10)
terra::values(r) <- seq_len(terra::ncell(r))
names(r) <- "gradient"
r2 <- r
terra::values(r2) <- rep(seq_len(10), each = 10)
names(r2) <- "stripes"
env <- c(r, r2)
p <- data.frame(x = c(1.5, 3.5, 5.5, 7.5),
                y = c(1.5, 3.5, 5.5, 7.5))
a <- data.frame(x = c(0.5, 2.5, 4.5, 6.5, 8.5),
                y = c(8.5, 6.5, 4.5, 2.5, 0.5))
train <- SDMtune::prepareSWD("Example species", env, p = p, a = a,
                             verbose = FALSE)
result <- screen_egv_variance(train, occ_folds = c(1, 1, 2, 2),
                              env = env, env_sample_n = 50,
                              seed = 1, verbose = FALSE)
result$keep
#> [1] "gradient" "stripes" 
head(result$diagnostics)
#>      check subset data_source      group fold used_for_selection variable
#> 1 training    all    training  presences <NA>               TRUE gradient
#> 2 training    all    training  presences <NA>               TRUE  stripes
#> 3 training    all    training background <NA>               TRUE gradient
#> 4 training    all    training background <NA>               TRUE  stripes
#> 5 training    all    training   combined <NA>               TRUE gradient
#> 6 training    all    training   combined <NA>               TRUE  stripes
#>   n_finite        sd failed reason
#> 1        4 23.237900  FALSE passed
#> 2        4  2.581989  FALSE passed
#> 3        5 34.785054  FALSE passed
#> 4        5  3.162278  FALSE passed
#> 5        9 28.416544  FALSE passed
#> 6        9  2.738613  FALSE passed
```
