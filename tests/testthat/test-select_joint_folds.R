make_joint_folds <- function() {
  list(
    pres = rep(letters[1:8], each = 5),
    bg = rep(letters[1:8], each = 20)
  )
}

test_that("joint selection is reproducible and uses the same fold IDs", {
  x <- make_joint_folds()
  first <- select_joint_folds(
    x$pres, x$bg, target_prop = 0.25,
    pres_bounds = c(0.20, 0.30), bg_bounds = c(0.20, 0.30),
    seed = 123, print_report = FALSE
  )
  second <- select_joint_folds(
    x$pres, x$bg, target_prop = 0.25,
    pres_bounds = c(0.20, 0.30), bg_bounds = c(0.20, 0.30),
    seed = 123, print_report = FALSE
  )

  expect_equal(first, second)
  expect_named(
    first,
    c("selected_fold_ids", "pres_selected_idx", "bg_selected_idx", "summary")
  )
  expect_length(first$pres_selected_idx, length(x$pres))
  expect_length(first$bg_selected_idx, length(x$bg))
  expect_setequal(unique(x$pres[first$pres_selected_idx]), first$selected_fold_ids)
  expect_setequal(unique(x$bg[first$bg_selected_idx]), first$selected_fold_ids)
  expect_equal(unname(first$summary$prop_pres_selected), 0.25)
  expect_equal(unname(first$summary$prop_bg_selected), 0.25)
})

test_that("summary counts agree with selection vectors", {
  x <- make_joint_folds()
  result <- select_joint_folds(
    x$pres, x$bg, seed = 42, print_report = FALSE
  )

  expect_equal(sum(result$pres_selected_idx), unname(result$summary$n_pres_selected))
  expect_equal(sum(result$bg_selected_idx), unname(result$summary$n_bg_selected))
  expect_equal(unname(result$summary$n_pres_in_scope), length(x$pres))
  expect_equal(unname(result$summary$n_bg_in_scope), length(x$bg))
  expect_true(unname(result$summary$prop_pres_selected) >= result$summary$bounds_pres[1])
  expect_true(unname(result$summary$prop_pres_selected) <= result$summary$bounds_pres[2])
  expect_true(unname(result$summary$prop_bg_selected) >= result$summary$bounds_bg[1])
  expect_true(unname(result$summary$prop_bg_selected) <= result$summary$bounds_bg[2])
})

test_that("missing IDs are excluded when include_na is false", {
  x <- make_joint_folds()
  x$pres[c(1, 2)] <- NA
  x$bg[c(1, 2, 3)] <- NA

  result <- select_joint_folds(
    x$pres, x$bg,
    include_na = FALSE, seed = 1, print_report = FALSE
  )

  expect_false(any(result$pres_selected_idx[is.na(x$pres)]))
  expect_false(any(result$bg_selected_idx[is.na(x$bg)]))
  expect_equal(unname(result$summary$n_pres_in_scope), sum(!is.na(x$pres)))
  expect_equal(unname(result$summary$n_bg_in_scope), sum(!is.na(x$bg)))
})

test_that("optional report can be suppressed or printed", {
  x <- make_joint_folds()
  expect_silent(select_joint_folds(x$pres, x$bg, seed = 1, print_report = FALSE))
  expect_output(
    select_joint_folds(x$pres, x$bg, seed = 1, print_report = TRUE),
    "Joint fold selection"
  )
})

test_that("select_joint_folds validates inputs", {
  x <- make_joint_folds()
  expect_error(select_joint_folds(list(1), x$bg), "pres_folds")
  expect_error(select_joint_folds(x$pres, list(1)), "bg_folds")
  expect_error(select_joint_folds(x$pres, x$bg, target_prop = 0), "target_prop")
  expect_error(select_joint_folds(x$pres, x$bg, pres_bounds = c(0.4, 0.2)), "pres_bounds")
  expect_error(select_joint_folds(x$pres, x$bg, bg_bounds = c(-0.1, 0.3)), "bg_bounds")
  expect_error(select_joint_folds(x$pres, x$bg, max_tries = 0), "max_tries")
  expect_error(select_joint_folds(x$pres, x$bg, bg_match_weight = -1), "bg_match_weight")
})
