#' Example ecogeographical-variable raster
#'
#' A small multilayer raster included with `sdmhelpers` for examples and
#' demonstrations of spatial workflows.
#'
#' @details
#' The raster is stored as a GeoTIFF file under the package's `extdata`
#' directory. It is not loaded automatically as an R data object.
#'
#' Locate the installed file with:
#'
#' ```r
#' egv_path <- system.file(
#'   "extdata",
#'   "example_egv.tif",
#'   package = "sdmhelpers"
#' )
#' ```
#'
#' Read it as a `terra::SpatRaster` with:
#'
#' ```r
#' example_egv <- terra::rast(egv_path)
#' ```
#'
#' @format A GeoTIFF file containing a multilayer raster in the coordinate
#'   reference system used by the package's spatial examples. Layers are
#'   described in `data(egv_names)`
#'
#' @source A reduced subset of ecogeographical variables prepared for
#'   `sdmhelpers` examples.
#'
#' @examples
#' egv_path <- system.file("extdata",
#' "example_egv.tif",
#' package = "sdmhelpers")
#'
#' if (nzchar(egv_path) && requireNamespace("terra", quietly = TRUE)) {
#'   example_egv <- terra::rast(egv_path)
#'   example_egv
#' }
#'
#' @name example_egv
NULL
