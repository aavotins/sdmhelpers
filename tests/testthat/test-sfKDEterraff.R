make_kde_fixture <- function() {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("fastfocal")

  ref <- terra::rast(
    xmin = 0, xmax = 100,
    ymin = 0, ymax = 100,
    resolution = 10,
    crs = "EPSG:3059"
  )
  terra::values(ref) <- 1

  points <- sf::st_as_sf(
    data.frame(
      x = c(35, 45, 55, 65),
      y = c(45, 55, 45, 55),
      weight = c(1, 2, 3, 4)
    ),
    coords = c("x", "y"),
    crs = 3059
  )

  list(ref = ref, points = points)
}

test_that("sfKDEterraff returns an aligned SpatRaster", {
  x <- make_kde_fixture()
  result <- sfKDEterraff(x$points, sigma = 15, ref = x$ref)

  expect_s4_class(result, "SpatRaster")
  expect_equal(as.vector(terra::ext(result)),as.vector(terra::ext(x$ref)))
  expect_equal(terra::res(result), terra::res(x$ref))
  expect_equal(terra::ncell(result), terra::ncell(x$ref))
  expect_true(all(terra::values(result) >= 0, na.rm = TRUE))
})

test_that("normalization options have their documented scale", {
  x <- make_kde_fixture()
  none <- sfKDEterraff(x$points, sigma = 15, ref = x$ref, normalize = "none")
  intensity <- sfKDEterraff(x$points, sigma = 15, ref = x$ref, normalize = "intensity")
  pdf <- sfKDEterraff(x$points, sigma = 15, ref = x$ref, normalize = "pdf")
  mean_g <- sfKDEterraff(x$points, sigma = 15, ref = x$ref, normalize = "meanG")

  cell_area <- prod(terra::res(x$ref))
  expect_equal(
    terra::values(intensity),
    terra::values(none) / cell_area,
    tolerance = 1e-10
  )
  expect_equal(
    terra::global(pdf, "sum", na.rm = TRUE)[1, 1] * cell_area,
    1,
    tolerance = 1e-8
  )
  expect_equal(
    terra::global(mean_g, "mean", na.rm = TRUE)[1, 1],
    1,
    tolerance = 1e-8
  )
})

test_that("a numeric weight field changes the result", {
  x <- make_kde_fixture()
  unweighted <- sfKDEterraff(x$points, sigma = 15, ref = x$ref)
  weighted <- sfKDEterraff(
    x$points, weight_field = "weight", sigma = 15, ref = x$ref
  )

  expect_gt(
    terra::global(weighted, "sum", na.rm = TRUE)[1, 1],
    terra::global(unweighted, "sum", na.rm = TRUE)[1, 1]
  )
})

test_that("masking respects NA cells in the reference raster", {
  x <- make_kde_fixture()
  values <- terra::values(x$ref)
  values[1:10] <- NA
  terra::values(x$ref) <- values

  result <- sfKDEterraff(
    x$points, sigma = 15, ref = x$ref, mask = TRUE
  )

  expect_true(all(is.na(terra::values(result)[1:10])))
})

test_that("sfKDEterraff validates geometry, sigma, weights, and normalization", {
  x <- make_kde_fixture()
  polygon <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(x$points)))

  expect_error(sfKDEterraff(polygon, sigma = 15, ref = x$ref), "POINT")
  expect_error(sfKDEterraff(x$points, sigma = 0, ref = x$ref), "sigma")
  expect_error(
    sfKDEterraff(x$points, weight_field = "missing", sigma = 15, ref = x$ref),
    "weight_field"
  )

  bad_weights <- x$points
  bad_weights$weight[1] <- NA_real_
  expect_error(
    sfKDEterraff(bad_weights, weight_field = "weight", sigma = 15, ref = x$ref),
    "weights"
  )

  expect_error(
    sfKDEterraff(x$points, sigma = 15, ref = x$ref, normalize = "bad"),
    "arg"
  )
})
