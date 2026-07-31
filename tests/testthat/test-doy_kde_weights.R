test_that("different input lengths do not produce recycling warnings", {
  A_dates <- as.Date(c(
    "2020-05-01",
    "2020-05-10",
    "2020-05-20",
    "2020-05-30",
    "2020-06-10"
  ))

  B_dates <- as.Date(c(
    "2021-05-05",
    "2021-06-05"
  ))

  expect_no_warning(
    result <- doy_kde_weights(
      A_dates = A_dates,
      B_dates = B_dates
    )
  )

  expect_length(result$weights, length(B_dates))
})

test_that("leap-day observations use a 366-day seasonal cycle", {
  A_dates <- as.Date(c(
    "2020-02-29",
    "2020-03-10",
    "2020-04-10"
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

  expect_length(result$weights, 2L)
})

test_that("doy_kde_weights returns one weight per B date", {
  A_dates <- as.Date(c(
    "2020-05-01",
    "2020-05-10",
    "2020-05-20",
    "2021-05-05",
    "2021-05-15"
  ))

  B_dates <- as.Date(c(
    "2022-04-01",
    "2022-05-10",
    "2022-07-01"
  ))

  result <- doy_kde_weights(
    A_dates = A_dates,
    B_dates = B_dates
  )

  expect_length(result$weights, length(B_dates))
  expect_true(all(result$weights >= 0))
  expect_true(all(result$weights <= 1))
})

test_that("doy_kde_weights rejects invalid dates", {
  expect_error(
    doy_kde_weights(
      A_dates = c("2020-05-01", NA),
      B_dates = "2020-05-10"
    )
  )
})
