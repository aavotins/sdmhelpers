make_activity_dates <- function() {
  as.Date(c(
    "2020-04-20", "2020-05-01", "2020-05-10", "2020-05-20",
    "2020-06-01", "2021-04-25", "2021-05-05", "2021-05-15",
    "2021-05-25", "2021-06-05"
  ))
}

test_that("doy_kde_weights returns a complete, bounded result", {
  B_dates <- as.Date(c("2022-04-01", "2022-05-10", "2022-07-01"))
  result <- doy_kde_weights(make_activity_dates(), B_dates)

  expect_named(
    result,
    c("weights_raw", "weights", "B_doy", "density", "scale", "n_days")
  )
  expect_length(result$weights_raw, length(B_dates))
  expect_length(result$weights, length(B_dates))
  expect_equal(result$B_doy, as.POSIXlt(B_dates)$yday + 1L)
  expect_true(all(is.finite(result$weights_raw)))
  expect_true(all(result$weights >= 0 & result$weights <= 1))
  expect_equal(result$scale, "minmax_A")
  expect_equal(result$n_days, 365L)
})

test_that("different input lengths do not produce recycling warnings", {
  expect_no_warning(
    result <- doy_kde_weights(
      A_dates = make_activity_dates(),
      B_dates = as.Date(c("2021-05-05", "2021-06-05"))
    )
  )
  expect_length(result$weights, 2L)
})

test_that("day 366 triggers a 366-day seasonal cycle", {
  A_dates <- as.Date(c(
    "2020-03-10",
    "2020-06-10",
    "2020-12-31"
  ))

  B_dates <- as.Date(c(
    "2021-03-01",
    "2021-04-01"
  ))

  expect_no_warning(
    result <- doy_kde_weights(
      A_dates = A_dates,
      B_dates = B_dates
    )
  )

  expect_equal(result$n_days, 366L)
  expect_length(result$density$x, 366L)
})

test_that("n_days and density-grid size can be supplied explicitly", {
  result <- doy_kde_weights(
    make_activity_dates(),
    as.Date("2022-05-10"),
    n_days = 365,
    n = 100
  )

  expect_equal(result$n_days, 365L)
  expect_length(result$density$x, 100L)
  expect_length(result$density$y, 100L)
})

test_that("percentile scaling respects eps clamping", {
  result <- doy_kde_weights(
    make_activity_dates(),
    as.Date(c("2022-01-01", "2022-05-10", "2022-12-31")),
    scale = "percentile_A",
    eps = 0.1
  )

  expect_equal(result$scale, "percentile_A")
  expect_true(all(result$weights >= 0.1 & result$weights <= 0.9))
})

test_that("doy_kde_weights validates arguments", {
  expect_error(doy_kde_weights(NULL, as.Date("2020-05-10")), "A_dates.*non-empty")
  expect_error(doy_kde_weights(make_activity_dates(), character()), "B_dates.*non-empty")
  expect_error(
    doy_kde_weights(c("2020-05-01", NA), "2020-05-10"),
    "A_dates.*NA"
  )
  expect_error(
    doy_kde_weights(make_activity_dates(), "2020-05-10", n_days = 364),
    "365 or 366"
  )
  expect_error(
    doy_kde_weights(make_activity_dates(), "2020-05-10", n = 9),
    "n.*>= 10"
  )
  expect_error(
    doy_kde_weights(
      make_activity_dates(), "2020-05-10",
      scale = "percentile_A", eps = 0.5
    ),
    "eps"
  )
})

test_that("n_days can be explicitly set to 366", {
  result <- doy_kde_weights(
    A_dates = as.Date(c(
      "2021-03-01",
      "2021-04-01",
      "2021-05-01"
    )),
    B_dates = as.Date(c(
      "2022-03-01",
      "2022-04-01"
    )),
    n_days = 366L
  )

  expect_equal(result$n_days, 366L)
  expect_length(result$density$x, 366L)
})
