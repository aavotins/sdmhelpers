test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("fold selection is reproducible with a seed", {
  folds <- rep(letters[1:8], each = 10)

  first <- select_folds_by_quota(
    folds = folds,
    seed = 123,
    print_report = FALSE
  )

  second <- select_folds_by_quota(
    folds = folds,
    seed = 123,
    print_report = FALSE
  )

  expect_equal(first, second)
})
