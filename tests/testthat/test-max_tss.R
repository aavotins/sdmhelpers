test_that("max_tss identifies complete separation", {
  result <- max_tss(
    p_pos = c(0.7, 0.8, 0.9),
    p_neg = c(0.1, 0.2, 0.3)
  )

  expect_named(result, c("tss", "threshold", "sensitivity", "specificity"))
  expect_equal(unname(result["tss"]), 1)
  expect_equal(unname(result["sensitivity"]), 1)
  expect_equal(unname(result["specificity"]), 1)
  expect_equal(unname(result["threshold"]), 0.7)
})

test_that("max_tss removes non-finite predictions", {
  expect_no_warning(
    result <- max_tss(
      p_pos = c(0.8, 0.9, NA, Inf),
      p_neg = c(0.1, 0.2, NA, -Inf)
    )
  )

  expect_true(all(is.finite(result)))
  expect_equal(unname(result["tss"]), 1)
})

test_that("max_tss returns missing results when one filtered class is empty", {
  result <- max_tss(
    p_pos = c(NA_real_, Inf),
    p_neg = c(0.1, 0.2)
  )

  expect_true(all(is.na(result)))
})

test_that("max_tss handles tied predictions consistently", {
  result <- max_tss(
    p_pos = c(0.8, 0.5),
    p_neg = c(0.5, 0.2)
  )

  expect_equal(unname(result["threshold"]), 0.8)
  expect_equal(unname(result["tss"]), 0.5)
})

test_that("max_tss validates its inputs", {
  expect_error(max_tss("high", c(0.1, 0.2)), "p_pos.*numeric")
  expect_error(max_tss(c(0.8, 0.9), factor(c("low", "high"))), "p_neg.*numeric")
})
