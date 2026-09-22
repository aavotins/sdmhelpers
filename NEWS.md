# sdmhelpers 0.0.0.9000

## Initial development version

* Added phenological weighting with `doy_kde_weights()`.
* Added seasonal and diel overlap estimation with `torus_overlap_many()`.
* Added fold-selection functions.
* Added TSS evaluation functions.
* Added fast raster-based spatial KDE.
* Added `phenology_overlap_weights()` to calculate kernel-density overlap
  between a target species' day-of-year distribution and multiple comparison
  species, with normalized overlap-based weights.
* Renamed `sfKDEterraff()` to `kde_surface()`.
* Added `screen_egv_variance()` to select only variable with sufficient variability.
* Created vignette.
