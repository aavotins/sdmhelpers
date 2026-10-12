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

test_that("kde_surface returns an aligned SpatRaster", {
  x <- make_kde_fixture()
  result <- kde_surface(x$points, sigma = 15, ref = x$ref)

  expect_s4_class(result, "SpatRaster")
  expect_equal(as.vector(terra::ext(result)),as.vector(terra::ext(x$ref)))
  expect_equal(terra::res(result), terra::res(x$ref))
  expect_equal(terra::ncell(result), terra::ncell(x$ref))
  expect_true(all(terra::values(result) >= 0, na.rm = TRUE))
})

test_that("normalization options have their documented scale", {
  x <- make_kde_fixture()
  none <- kde_surface(x$points, sigma = 15, ref = x$ref, normalize = "none")
  intensity <- kde_surface(x$points, sigma = 15, ref = x$ref, normalize = "intensity")
  pdf <- kde_surface(x$points, sigma = 15, ref = x$ref, normalize = "pdf")
  mean_g <- kde_surface(x$points, sigma = 15, ref = x$ref, normalize = "meanG")

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
  unweighted <- kde_surface(x$points, sigma = 15, ref = x$ref)
  weighted <- kde_surface(
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

  result <- kde_surface(
    x$points, sigma = 15, ref = x$ref, mask = TRUE
  )

  expect_true(all(is.na(terra::values(result)[1:10])))
})

test_that("kde_surface validates geometry, sigma, weights, and normalization", {
  x <- make_kde_fixture()
  polygon <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(x$points)))

  expect_error(
    kde_surface(polygon, sigma = 15, ref = x$ref),
    "POINT"
  )

  expect_error(
    kde_surface(x$points, sigma = 0, ref = x$ref),
    "sigma"
  )

  expect_error(
    kde_surface(
      x$points,
      weight_field = "missing",
      sigma = 15,
      ref = x$ref
    ),
    "weight_field"
  )

  bad_weights <- x$points
  bad_weights$weight[1] <- NA_real_

  expect_error(
    kde_surface(
      bad_weights,
      weight_field = "weight",
      sigma = 15,
      ref = x$ref
    ),
    "Weights must be finite numeric"
  )

  expect_error(
    kde_surface(
      x$points,
      sigma = 15,
      ref = x$ref,
      normalize = "bad"
    ),
    "arg"
  )
})

# -----------------------------------------------------------------------------
# Regression tests for input validations (Issues 21, 22, 23, mask argument)
# -----------------------------------------------------------------------------

test_that("non-square raster cells are rejected before KDE calculation", {
  x <- make_kde_fixture()

  rectangular_ref <- terra::rast(
    xmin = 0, xmax = 100,
    ymin = 0, ymax = 200,
    resolution = c(10, 20),
    crs = "EPSG:3059"
  )
  terra::values(rectangular_ref) <- 1

  expect_equal(terra::res(rectangular_ref), c(10, 20))
  expect_error(
    kde_surface(x$points, sigma = 15, ref = rectangular_ref),
    "square raster cells"
  )

  # The same code continues to accept square projected cells.
  expect_s4_class(
    kde_surface(x$points, sigma = 15, ref = x$ref),
    "SpatRaster"
  )
})

test_that("zero effective weights cannot be normalized as PDF or meanG", {
  x <- make_kde_fixture()
  x$points$weight <- 0

  expect_error(
    kde_surface(
      x$points, weight_field = "weight", sigma = 15,
      ref = x$ref, normalize = "pdf"
    ),
    "non-positive or non-finite total mass"
  )

  expect_error(
    kde_surface(
      x$points, weight_field = "weight", sigma = 15,
      ref = x$ref, normalize = "meanG"
    ),
    "non-positive or non-finite mean"
  )

  # A zero-valued KDE is still a valid unnormalized or intensity surface.
  none <- kde_surface(
    x$points, weight_field = "weight", sigma = 15,
    ref = x$ref, normalize = "none"
  )
  intensity <- kde_surface(
    x$points, weight_field = "weight", sigma = 15,
    ref = x$ref, normalize = "intensity"
  )
  expect_s4_class(none, "SpatRaster")
  expect_true(any(is.finite(terra::values(none))))
  expect_true(all(terra::values(none) == 0, na.rm = TRUE))
  expect_equal(
    terra::values(intensity),
    terra::values(none) / prod(terra::res(x$ref)),
    tolerance = 1e-10
  )
})

test_that("masking can produce an empty retained KDE support", {
  x <- make_kde_fixture()
  # Mask all cells after smoothing, so normalization is undefined.
  terra::values(x$ref) <- NA_real_

  expect_error(
    kde_surface(
      x$points, sigma = 15, ref = x$ref,
      normalize = "pdf", mask = TRUE
    ),
    "non-positive or non-finite total mass"
  )

  expect_error(
    kde_surface(
      x$points, sigma = 15, ref = x$ref,
      normalize = "meanG", mask = TRUE
    ),
    "non-positive or non-finite mean"
  )
})

test_that("invalid backend names and types are rejected", {
  x <- make_kde_fixture()

  for (bad_engine in list(
    "fftt", "unknown", 1, NA_character_,
    character(0), c("fft", "cpp")
  )) {
    expect_error(
      kde_surface(
        x$points, sigma = 15, ref = x$ref,
        engine = bad_engine
      ),
      "`engine` must be one of"
    )
  }

  # A valid engine remains accepted.
  expect_s4_class(
    kde_surface(x$points, sigma = 15, ref = x$ref, engine = "auto"),
    "SpatRaster"
  )
})

test_that("mask must be one non-missing logical value", {
  x <- make_kde_fixture()

  for (bad_mask in list(
    1, "TRUE", NA, logical(0), c(TRUE, FALSE)
  )) {
    expect_error(
      kde_surface(
        x$points, sigma = 15, ref = x$ref,
        mask = bad_mask
      ),
      "`mask` must be TRUE or FALSE"
    )
  }

  expect_s4_class(
    kde_surface(x$points, sigma = 15, ref = x$ref, mask = TRUE),
    "SpatRaster"
  )
})
