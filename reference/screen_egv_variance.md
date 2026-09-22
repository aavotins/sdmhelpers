# Screen Environmental Predictors for Zero Variance

Identifies environmental predictors that have zero or insufficient
variance in data subsets relevant to species distribution model
training, cross-validation, background characterization, and optionally
independent testing.

## Usage

``` r
screen_egv_variance(
  train,
  test = NULL,
  occ_folds = NULL,
  env = NULL,
  env_sample_n = 10000L,
  seed = NULL,
  sd_tol = 0,
  include_test_bg = TRUE,
  check_test = FALSE,
  verbose = TRUE
)
```

## Arguments

- train:

  An `SWD` object, typically created with
  [`SDMtune::prepareSWD()`](https://consbiol-unibern.github.io/SDMtune/reference/prepareSWD.html),
  containing the training presence and absence/background records.

- test:

  An optional `SWD` object containing independent testing data. When
  supplied, its background records can be included in the pooled
  background variance check. The complete independent testing dataset
  can additionally be checked by setting `check_test = TRUE`.

- occ_folds:

  An optional vector giving the cross-validation fold assigned to each
  presence record in `train`. Its length must equal the number of
  presence records in `train`. For each fold, predictor variance is
  evaluated using the presence records belonging to all other folds,
  corresponding to the training partition when that fold is held out.

- env:

  An optional
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  containing the environmental predictors. When supplied, `env_sample_n`
  random raster cells are sampled and predictor variance is evaluated
  across the sampled environmental domain.

- env_sample_n:

  Positive integer. Number of cells randomly sampled from `env` when
  `env` is supplied. Defaults to `10000`.

- seed:

  Optional integer used for reproducible sampling of `env`. The previous
  state of R's random-number generator is restored before the function
  exits.

- sd_tol:

  Non-negative numeric value defining the minimum acceptable standard
  deviation. A predictor is flagged when its standard deviation is less
  than or equal to `sd_tol`. The default, `0`, removes only predictors
  with exactly zero variance, together with predictors for which
  variance cannot be estimated.

- include_test_bg:

  Logical. If `TRUE` and `test` is supplied, background records from the
  independent testing dataset are pooled with training background
  records for the background variance check. Defaults to `TRUE`.

- check_test:

  Logical. If `TRUE`, predictor variance is additionally evaluated
  across all records in the independent testing `SWD` object. Defaults
  to `FALSE`.

- verbose:

  Logical. If `TRUE`, print a short summary of the screening results.
  Defaults to `TRUE`.

## Value

A named list with the following components:

- keep:

  Character vector containing predictors that passed all requested
  variance checks.

- remove:

  Character vector containing predictors that failed at least one
  variance check.

- removed_by_check:

  Named list showing which predictors were flagged by each individual
  check.

- diagnostics:

  A `data.frame` containing the check, subset, predictor name, number of
  finite observations, standard deviation, failure status, and reason
  for every predictor evaluated.

## Details

The function is intended as a preprocessing check for environmental
predictors used in species distribution modelling. A predictor is
retained only when it passes every variance check that is requested.

By default, the following subsets are evaluated:

- all presence and absence/background records in the training data;

- all training presence records;

- pooled training background records and, when `include_test_bg = TRUE`,
  independent testing background records;

- each cross-validation presence training partition when `occ_folds` is
  supplied; and

- a random sample of the environmental domain when `env` is supplied.

For the cross-validation check, the function does not evaluate the
records within a held-out fold alone. For each fold `k`, it evaluates
the records for which `occ_folds != k`. These are the presence records
available for model fitting when fold `k` is used for validation.

If `check_test = TRUE`, all records in `test` are also evaluated as a
separate subset.

Standard deviations are calculated using finite values only. Predictors
having fewer than two finite observations in a checked subset are
flagged because their variance cannot be estimated.

Predictor names are taken from `train@data`. The same predictors must be
present in `test@data`, when `test` is supplied, and in `env`, when
`env` is supplied.

The environmental raster itself is not modified. Instead, the function
returns vectors containing the predictors that passed and failed
screening. These can subsequently be used to subset an
environmental-variable table, an `SWD` object, or a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

## See also

[`SDMtune::prepareSWD()`](https://consbiol-unibern.github.io/SDMtune/reference/prepareSWD.html),
[`terra::spatSample()`](https://rspatial.github.io/terra/reference/sample.html)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- screen_egv_variance(
  train = trenin_dati,
  test = testa_dati,
  occ_folds = block_folds$occs.grp,
  env = vide,
  env_sample_n = 10000,
  seed = 1
)

result$keep
result$remove
result$removed_by_check

# Retain only predictors passing all checks
videi <- videi[
  videi$layername %in% result$keep,
]

vide <- terra::rast(
  paste0(
    "./RasterGrids_100m/2024/Scaled/",
    videi$new_filename
  )
)
} # }
```
