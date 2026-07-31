test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("max_tss identifies complete separation", {
  result <- max_tss(
    p_pos = c(0.7, 0.8, 0.9),
    p_neg = c(0.1, 0.2, 0.3)
  )

  expect_equal(unname(result["tss"]), 1)
  expect_equal(unname(result["sensitivity"]), 1)
  expect_equal(unname(result["specificity"]), 1)
})

test_that("max_tss removes non-finite predictions", {
  result <- max_tss(
    p_pos = c(0.8, 0.9, NA, Inf),
    p_neg = c(0.1, 0.2, NA, -Inf)
  )

  expect_true(is.finite(result["tss"]))
})

test_that("max_tss returns missing results when one class is empty", {
  result <- max_tss(
    p_pos = numeric(),
    p_neg = c(0.1, 0.2)
  )

  expect_true(all(is.na(result)))
})


