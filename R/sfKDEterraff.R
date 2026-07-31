#' Fast Gaussian KDE from sf POINT geometries
#'
#' Computes a weighted or unweighted Gaussian kernel density surface by first
#' rasterizing point observations to a reference raster and then applying
#' Gaussian smoothing with `fastfocal::fastfocal()`.
#'
#' This is a binned KDE because point coordinates are first aggregated to raster
#' cells. It can be considerably faster and less memory demanding than
#' evaluating a continuous kernel density directly at every output cell,
#' particularly for large rasters or large bandwidths.
#'
#' @param x An `sf` or `sfc` object containing only single-part `POINT`
#'   geometries. The coordinates of `x` must use the same coordinate reference
#'   system as `ref`. The function does not reproject the points.
#'
#' @param weight_field Optional character string giving the name of a numeric
#'   column in `x` containing point weights. If `NULL`, every point is assigned
#'   weight 1. Weights belonging to points in the same raster cell are summed
#'   before smoothing. Weights must be finite numeric values without missing
#'   values.
#'
#' @param sigma A single finite numeric value greater than zero. Gaussian kernel
#'   standard deviation, or bandwidth, expressed in the map units of `ref`.
#'   For example, if `ref` uses a projected coordinate reference system measured
#'   in metres, `sigma = 3000` represents a bandwidth of 3000 metres.
#'
#' @param ref A `terra::SpatRaster` defining the extent, resolution, origin,
#'   dimensions, and coordinate reference system of the output raster.
#'
#' @param normalize Character string specifying the output scaling. One of:
#'
#'   `"none"` returns the Gaussian-smoothed rasterized counts or weight sums
#'   without additional scaling.
#'
#'   `"pdf"` divides the smoothed surface by its raster integral so that it
#'   integrates to approximately 1 over its non-missing support. Values have
#'   units of inverse map-area units, such as per square metre when the raster
#'   coordinate units are metres.
#'
#'   `"intensity"` divides the smoothed values by raster-cell area. For
#'   unweighted data, the result may be interpreted as smoothed points per unit
#'   area.
#'
#'   `"meanG"` divides the surface by its global non-missing mean. The returned
#'   relative surface therefore has mean 1 and may be useful as an effort,
#'   availability, or sampling-bias layer.
#'
#' @param mask Logical. If `TRUE`, the KDE is masked using the non-`NA`
#'   footprint of `ref`. Masking is applied after smoothing but before
#'   normalization. Therefore, `"pdf"` and `"meanG"` normalization are
#'   calculated over the retained raster footprint.
#'
#' @param engine Character string passed to the `engine` argument of
#'   `fastfocal::fastfocal()`. Common options include `"auto"`, `"fft"`, and
#'   `"cpp"`, depending on the installed version of `fastfocal`. The default,
#'   `"auto"`, allows `fastfocal` to select an appropriate backend.
#'
#' @return A single-layer `terra::SpatRaster` aligned with `ref`.
#'
#'   With `normalize = "none"`, values are smoothed rasterized counts or weight
#'   sums.
#'
#'   With `normalize = "pdf"`, values form a density surface that integrates to
#'   approximately 1.
#'
#'   With `normalize = "intensity"`, values are smoothed counts or weights per
#'   unit area.
#'
#'   With `normalize = "meanG"`, values form a relative surface with global
#'   non-missing mean 1.
#'
#' @details
#' The function performs the following operations:
#'
#' 1. Converts `x` to a `terra::SpatVector`.
#'
#' 2. Assigns each point either weight 1 or the value stored in
#'    `weight_field`.
#'
#' 3. Rasterizes the points to `ref`, summing point weights within each raster
#'    cell and assigning zero to cells without points.
#'
#' 4. Applies Gaussian smoothing using `fastfocal::fastfocal()`.
#'
#' 5. Optionally masks the result using `ref`.
#'
#' 6. Applies the requested normalization.
#'
#' Because points are rasterized before smoothing, the function calculates a
#' binned approximation to a continuous KDE. The approximation generally
#' improves as raster resolution becomes finer relative to `sigma`. Results may
#' differ from continuous-coordinate KDE estimators when the bandwidth is small
#' relative to raster-cell size.
#'
#' The function assumes that `x` and `ref` use the same coordinate reference
#' system. No CRS comparison or coordinate transformation is performed. Points
#' outside the extent of `ref` do not contribute to the output.
#'
#' Cell area is calculated as `prod(terra::res(ref))`. Therefore, `"pdf"` and
#' `"intensity"` assume a regular projected raster with linear coordinate
#' units. These normalizations should generally not be used with
#' longitude-latitude rasters because cell area varies spatially.
#'
#' Masking is performed after Gaussian smoothing. Consequently, the mask
#' determines which cells are retained and included in normalization, but it
#' does not prevent smoothing across internal `NA` boundaries in `ref`.
#'
#' Gaussian kernels near the raster extent may lose part of their mass because
#' the kernel extends beyond the output grid. The `"pdf"` option renormalizes
#' the retained surface so that its raster integral is approximately 1.
#'
#' For `"pdf"` and `"meanG"` normalization, the input should contain at least
#' one point with a nonzero effective weight inside the extent of `ref`.
#' Otherwise, the normalization denominator may be zero or non-finite.
#'
#' Negative weights are permitted by the current input checks, but they can
#' produce negative raster values and may invalidate probability-density or
#' relative-bias interpretations.
#'
#' @section Dependencies:
#' This function requires the `sf`, `terra`, and `fastfocal` packages.
#'
#' @seealso
#' `fastfocal::fastfocal()`, `terra::rasterize()`, `terra::mask()`,
#' `terra::global()`
#'
#' @examples
#' \dontrun{
#' # Create a projected raster template with 100-m cells
#' ref <- terra::rast(
#'   xmin = 0,
#'   xmax = 10000,
#'   ymin = 0,
#'   ymax = 10000,
#'   resolution = 100,
#'   crs = "EPSG:3857"
#' )
#'
#' # Give the template non-missing values so it can also be used as a mask
#' terra::values(ref) <- 1
#'
#' # Create an irregular raster footprint
#' xy <- terra::xyFromCell(ref, seq_len(terra::ncell(ref)))
#' distance_from_centre <- sqrt(
#'   (xy[, 1] - 5000)^2 +
#'   (xy[, 2] - 5000)^2
#' )
#' ref[distance_from_centre > 4500] <- NA
#'
#' # Generate example points in the same projected CRS
#' set.seed(123)
#'
#' dat <- data.frame(
#'   x = stats::rnorm(200, mean = 5000, sd = 1500),
#'   y = stats::rnorm(200, mean = 5000, sd = 1500),
#'   weight = stats::runif(200, min = 0.5, max = 2)
#' )
#'
#' pts <- sf::st_as_sf(
#'   dat,
#'   coords = c("x", "y"),
#'   crs = 3857
#' )
#'
#' # Unweighted KDE normalized as a probability density
#' kde_pdf <- sfKDEterraff(
#'   x = pts,
#'   sigma = 750,
#'   ref = ref,
#'   normalize = "pdf",
#'   mask = TRUE
#' )
#'
#' terra::plot(kde_pdf)
#'
#' # Check that the raster density integrates to approximately 1
#' cell_area <- prod(terra::res(kde_pdf))
#'
#' terra::global(
#'   kde_pdf,
#'   fun = "sum",
#'   na.rm = TRUE
#' )[1, 1] * cell_area
#'
#' # Weighted KDE standardized to global mean 1
#' kde_bias <- sfKDEterraff(
#'   x = pts,
#'   weight_field = "weight",
#'   sigma = 750,
#'   ref = ref,
#'   normalize = "meanG",
#'   mask = TRUE
#' )
#'
#' terra::global(
#'   kde_bias,
#'   fun = "mean",
#'   na.rm = TRUE
#' )
#'
#' terra::plot(kde_bias)
#'
#' # Smoothed point intensity per square map unit
#' kde_intensity <- sfKDEterraff(
#'   x = pts,
#'   sigma = 750,
#'   ref = ref,
#'   normalize = "intensity",
#'   mask = TRUE,
#'   engine = "auto"
#' )
#'
#' terra::plot(kde_intensity)
#' }
#'
#' @export
sfKDEterraff <- function(
    x,
    weight_field = NULL,
    sigma,
    ref,
    normalize = c("none", "pdf", "intensity", "meanG"),
    mask = FALSE,
    engine = "auto"
) {
  normalize <- match.arg(normalize)
  
  # --- checks ---
  
  stopifnot(inherits(ref, "SpatRaster"))
  stopifnot(inherits(x, c("sf", "sfc")))
  if (unique(as.character(sf::st_geometry_type(x))) != "POINT")
    stop("`x` must be single-part POINT geometry")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
    stop("`sigma` must be a single numeric value > 0")
  
  # --- weights ---
  
  if (!is.null(weight_field)) {
    if (!inherits(x, "sf")) x <- sf::st_as_sf(x)
    if (!weight_field %in% names(x)) stop("`weight_field` not found in `x`")
    w <- x[[weight_field]]
    if (!is.numeric(w) || anyNA(w) || any(!is.finite(w))) stop("weights must be finite numeric, no NA")
  } else {
    w <- rep(1, nrow(sf::st_as_sf(x)))
  }
  
  # --- rasterize (sum weights per cell) ---
  
  v <- terra::vect(sf::st_as_sf(x))
  v$w <- w
  r_counts <- terra::rasterize(v, ref, field = "w", fun = "sum", background = 0)
  
  # --- Gaussian smoothing (d = sigma in map units) ---
  
  # fastfocal supports gaussian windows and auto-selects FFT for big kernels.
  
  r <- fastfocal::fastfocal(
    x      = r_counts,
    d      = sigma,
    w      = "gaussian",
    fun    = "sum",
    engine = engine
  )
  
  # mask BEFORE normalization so PDF integrates to 1 over the masked footprint
  
  if (mask) r <- terra::mask(r, ref)
  
  # --- normalization ---
  
  cell_area <- prod(terra::res(ref))
  
  if (normalize == "pdf") {
    total_mass <- terra::global(r, "sum", na.rm = TRUE)[1, 1] * cell_area
    r <- r / total_mass
  } else if (normalize == "intensity") {
    # points per unit area
    r <- r / cell_area
  } else if (normalize == "meanG") {
    m <- terra::global(r, "mean", na.rm = TRUE)[1, 1]
    r <- r / m
  }
  
  r
}
