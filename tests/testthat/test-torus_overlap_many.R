make_torus_data <- function() {
  skip_if_not_installed("hms")

  reference <- data.frame(
    date = as.Date(c("2020-05-01", "2020-05-03", "2020-05-05", "2020-05-07", "2020-05-09")),
    time = hms::as_hms(c("06:00:00", "06:30:00", "07:00:00", "07:30:00", "08:00:00"))
  )

  comparison <- rbind(
    data.frame(
      group = "similar",
      date = reference$date,
      time = reference$time
    ),
    data.frame(
      group = "different",
      date = as.Date(c("2020-11-01", "2020-11-03", "2020-11-05", "2020-11-07", "2020-11-09")),
      time = hms::as_hms(c("18:00:00", "18:30:00", "19:00:00", "19:30:00", "20:00:00"))
    )
  )

  list(reference = reference, comparison = comparison)
}

test_that("torus_overlap_many returns one result per comparison group", {
  x <- make_torus_data()
  result <- torus_overlap_many(
    reference_data = x$reference,
    comparison_data = x$comparison,
    comparison_group_col = "group",
    reference_date_col = "date",
    reference_time_col = "time",
    min_n = 5,
    date_resolution_days = 15,
    time_resolution_minutes = 60
  )

  expect_s3_class(result, "torus_overlap_many")
  expect_named(result, c("overlap", "settings"))
  expect_equal(nrow(result$overlap), 2L)
  expect_named(
    result$overlap,
    c("reference", "comparison", "n_reference", "n_comparison", "overlap", "status")
  )
  expect_true(all(result$overlap$status == "OK"))
  expect_true(all(result$overlap$overlap >= 0 & result$overlap$overlap <= 1))

  similar <- result$overlap$overlap[result$overlap$comparison == "similar"]
  different <- result$overlap$overlap[result$overlap$comparison == "different"]
  expect_equal(similar, 1, tolerance = 1e-8)
  expect_gt(similar, different)
})

test_that("groups below min_n are retained with NA overlap", {
  x <- make_torus_data()
  x$comparison <- rbind(
    x$comparison,
    data.frame(
      group = "small",
      date = as.Date(c("2020-05-01", "2020-05-02")),
      time = hms::as_hms(c("06:00:00", "07:00:00"))
    )
  )

  result <- torus_overlap_many(
    x$reference, x$comparison, "group", "date", "time",
    min_n = 5, date_resolution_days = 15,
    time_resolution_minutes = 60
  )

  small <- result$overlap[result$overlap$comparison == "small", ]
  expect_equal(small$n_comparison, 2L)
  expect_true(is.na(small$overlap))
  expect_match(small$status, "Too few observations")
})

test_that("return_density adds a well-formed long table", {
  x <- make_torus_data()
  result <- torus_overlap_many(
    x$reference, x$comparison, "group", "date", "time",
    min_n = 5, date_resolution_days = 30,
    time_resolution_minutes = 120,
    return_density = TRUE
  )

  expect_named(result, c("overlap", "settings", "density"))
  expect_s3_class(result$density, "data.frame")
  expect_gt(nrow(result$density), 0L)
  expect_true(all(result$density$reference_density >= 0))
  expect_true(all(result$density$comparison_density >= 0))
  expect_true(all(result$density$shared_density >= 0))
})

test_that("invalid rows are removed and recorded in settings", {
  x <- make_torus_data()
  x$comparison$date[1] <- NA

  result <- torus_overlap_many(
    x$reference, x$comparison, "group", "date", "time",
    min_n = 4, date_resolution_days = 30,
    time_resolution_minutes = 120
  )

  expect_equal(
    result$overlap$n_comparison[result$overlap$comparison == "similar"],
    4L
  )
})

test_that("torus_overlap_many validates data and parameter inputs", {
  x <- make_torus_data()
  expect_error(
    torus_overlap_many(1, x$comparison, "group", "date", "time"),
    "reference_data.*data frame"
  )
  expect_error(
    torus_overlap_many(x$reference, x$comparison, "missing", "date", "time"),
    "missing from.*comparison_data"
  )

  bad_reference <- x$reference
  bad_reference$date <- as.character(bad_reference$date)
  expect_error(
    torus_overlap_many(bad_reference, x$comparison, "group", "date", "time"),
    "Date"
  )

  expect_error(
    torus_overlap_many(
      x$reference, x$comparison, "group", "date", "time",
      season_bw_days = 0
    ),
    "season_bw_days"
  )
  expect_error(
    torus_overlap_many(
      x$reference, x$comparison, "group", "date", "time",
      min_n = 1
    ),
    "min_n"
  )
})
