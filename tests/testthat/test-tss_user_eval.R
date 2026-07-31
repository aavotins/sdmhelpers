make_enmeval_vars <- function(validation_bg = "partition") {
  list(
    occs.train.pred = c(0.8, 0.9),
    occs.val.pred = c(0.7, 0.85),
    bg.train.pred = c(0.1, 0.2),
    bg.val.pred = c(0.15, 0.3),
    other.settings = list(validation.bg = validation_bg)
  )
}

test_that("tss_user_eval returns the ENMeval-compatible columns", {
  result <- tss_user_eval(make_enmeval_vars())

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_named(
    result,
    c("tss.val", "thr.val", "sens.val", "spec.val", "tss.train")
  )
  expect_equal(result$tss.val, 1)
  expect_equal(result$tss.train, 1)
  expect_equal(result$sens.val, 1)
  expect_equal(result$spec.val, 1)
})

test_that("full validation background combines training and validation background", {
  vars <- make_enmeval_vars("full")
  vars$bg.train.pred <- c(0.1, 0.75)

  result <- tss_user_eval(vars)
  expected <- max_tss(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))

  expect_equal(result$tss.val, unname(expected["tss"]))
  expect_equal(result$thr.val, unname(expected["threshold"]))
  expect_equal(result$sens.val, unname(expected["sensitivity"]))
  expect_equal(result$spec.val, unname(expected["specificity"]))
})

test_that("other validation settings use validation background only", {
  vars <- make_enmeval_vars(validation_bg = "partition")
  result <- tss_user_eval(vars)
  expected <- max_tss(vars$occs.val.pred, vars$bg.val.pred)

  expect_equal(result$tss.val, unname(expected["tss"]))
})

test_that("tss_user_eval propagates missing-class results", {
  vars <- make_enmeval_vars()
  vars$occs.val.pred <- NA_real_
  result <- tss_user_eval(vars)

  expect_true(all(is.na(result[c("tss.val", "thr.val", "sens.val", "spec.val")])))
  expect_true(is.finite(result$tss.train))
})

test_that("tss_user_eval validates its input structure", {
  expect_error(tss_user_eval(1), "vars.*list")
  expect_error(tss_user_eval(list()), "missing.*occs.train.pred")

  vars <- make_enmeval_vars()
  vars$bg.val.pred <- factor(c("low", "high"))
  expect_error(tss_user_eval(vars), "must be numeric: bg\\.val\\.pred")
})
