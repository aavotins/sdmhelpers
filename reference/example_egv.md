# Example ecogeographical-variable raster

A small multilayer raster included with `sdmhelpers` for examples and
demonstrations of spatial workflows.

## Format

A GeoTIFF file containing a multilayer raster in the coordinate
reference system used by the package's spatial examples. Layers are
described in `data(egv_names)`

## Source

A reduced subset of ecogeographical variables prepared for `sdmhelpers`
examples.

## Details

The raster is stored as a GeoTIFF file under the package's `extdata`
directory. It is not loaded automatically as an R data object.

Locate the installed file with:

    egv_path <- system.file(
      "extdata",
      "example_egv.tif",
      package = "sdmhelpers"
    )

Read it as a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with:

    example_egv <- terra::rast(egv_path)

## Examples

``` r
egv_path <- system.file("extdata",
"example_egv.tif",
package = "sdmhelpers")

if (nzchar(egv_path) && requireNamespace("terra", quietly = TRUE)) {
  example_egv <- terra::rast(egv_path)
  example_egv
}
#> class       : SpatRaster
#> size        : 650, 650, 5  (nrow, ncol, nlyr)
#> resolution  : 100, 100  (x, y)
#> extent      : 455000, 520000, 280000, 345000  (xmin, xmax, ymin, ymax)
#> coord. ref. : LKS-92 / Latvia TM (EPSG:3059)
#> source      : example_egv.tif
#> names       :   egv_280,   egv_293,   egv_302,   egv_385,   egv_400
#> min values  : -0.466794, -0.719917, -0.397964, -1.461203,  -0.37327
#> max values  :   3.70894,  7.054723,  8.359298,  5.206886, 11.300995
```
