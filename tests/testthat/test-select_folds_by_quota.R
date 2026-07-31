test_that("fold selection is reproducible and satisfies an exact quota", {
  folds <- rep(letters[1:8], each = 10)

  first <- select_folds_by_quota(
    folds, target_prop = 0.25, bounds = c(0.20, 0.30),
    seed = 123, print_report = FALSE
  )
  second <- select_folds_by_quota(
    folds, target_prop = 0.25, bounds = c(0.20, 0.30),
    seed = 123, print_report = FALSE
  )

  expect_equal(first, second)
  expect_s3_class(first, "fold_quota_selection")
  expect_named(first, c("selected_fold_ids", "selected_idx", "summary"))
  expect_length(first$selected_idx, length(folds))
  expect_equal(sum(first$selected_idx), first$summary$n_rows_selected)
  expect_equal(first$summary$prop_selected, 0.25)
  expect_true(first$summary$within_bounds)
  expect_true(first$summary$exact_target)
  expect_equal(first$summary$n_groups_selected, 2L)
})

test_that("whole folds rather than individual rows are selected", {
  folds <- rep(c("large", "medium", "small"), times = c(6, 3, 1))
  result <- select_folds_by_quota(
    folds, target_prop = 0.3, bounds = c(0.2, 0.4),
    seed = 1, print_report = FALSE
  )

  for (id in unique(folds)) {
    selected <- result$selected_idx[folds == id]
    expect_true(all(selected) || !any(selected))
  }
  expect_setequal(unique(folds[result$selected_idx]), result$selected_fold_ids)
})

test_that("missing folds are excluded or treated as a group", {
  folds <- c(rep("a", 4), rep("b", 4), NA, NA)

  excluded <- select_folds_by_quota(
    folds, target_prop = 0.5, bounds = c(0.4, 0.6),
    include_na = FALSE, seed = 2, print_report = FALSE
  )
  expect_false(any(excluded$selected_idx[is.na(folds)]))
  expect_equal(excluded$summary$n_rows_total, 8L)

  included <- select_folds_by_quota(
    folds, target_prop = 0.2, bounds = c(0.2, 0.2),
    include_na = TRUE, seed = 2, print_report = FALSE
  )
  expect_equal(included$summary$n_rows_total, 10L)
  expect_true("<NA>" %in% included$selected_fold_ids)
  expect_true(all(included$selected_idx[is.na(folds)]))
})

test_that("optional report can be suppressed or printed", {
  folds <- rep(letters[1:4], each = 5)
  expect_silent(select_folds_by_quota(folds, seed = 1, print_report = FALSE))
  expect_output(
    select_folds_by_quota(folds, seed = 1, print_report = TRUE),
    "Fold selection summary"
  )
})

test_that("select_folds_by_quota validates inputs", {
  expect_error(select_folds_by_quota(NULL), "folds")
  expect_error(select_folds_by_quota(list(a = 1)), "atomic")
  expect_error(select_folds_by_quota(character()),"at least one value")
  expect_error(select_folds_by_quota(letters, target_prop = 1), "target_prop")
  expect_error(select_folds_by_quota(letters, bounds = c(0.3, 0.2)), "bounds")
  expect_error(select_folds_by_quota(letters, max_tries = 0), "max_tries")
  expect_error(select_folds_by_quota(letters, include_na = NA), "include_na")
  expect_error(select_folds_by_quota(letters, seed = c(1, 2)), "seed")
})
