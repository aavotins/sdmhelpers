# Changelog

## sdmhelpers 0.1.0

### First stable release

- Introduced and documented tools supporting species distribution
  modelling.
- Improved predictor variance screening across training and
  cross-validation subsets.
- Refined spatial kernel-density and temporal overlap functions.
- Added utilities for independent spatial partitioning and model
  evaluation.
- Expanded automated tests and the introductory vignette.

## sdmhelpers 0.0.0.9000

### Initial development version

- Added phenological weighting with
  [`doy_kde_weights()`](https://aavotins.github.io/sdmhelpers/reference/doy_kde_weights.md).
- Added seasonal and diel overlap estimation with
  [`torus_overlap_many()`](https://aavotins.github.io/sdmhelpers/reference/torus_overlap_many.md).
- Added fold-selection functions.
- Added TSS evaluation functions.
- Added fast raster-based spatial KDE.
- Added
  [`phenology_overlap_weights()`](https://aavotins.github.io/sdmhelpers/reference/phenology_overlap_weights.md)
  to calculate kernel-density overlap between a target species’
  day-of-year distribution and multiple comparison species, with
  normalized overlap-based weights.
- Renamed `sfKDEterraff()` to
  [`kde_surface()`](https://aavotins.github.io/sdmhelpers/reference/kde_surface.md).
- Added
  [`screen_egv_variance()`](https://aavotins.github.io/sdmhelpers/reference/screen_egv_variance.md)
  to select only variable with sufficient variability.
- Created vignette.
