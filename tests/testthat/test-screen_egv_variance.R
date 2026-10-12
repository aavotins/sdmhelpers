# Tests for screen_egv_variance(), including legacy checks and revised API.
# The primary fixtures are synthetic to guarantee known variance patterns.
# A separate integration section uses the package example data.

make_variance_swd <- function(data, pa, species = "Example") {
  testthat::skip_if_not_installed("SDMtune")
  data <- as.data.frame(data)
  methods::new(
    "SWD",
    species = species,
    coords = data.frame(
      X = seq_len(nrow(data)),
      Y = seq_len(nrow(data))
    ),
    data = data,
    pa = as.numeric(pa)
  )
}

make_variance_fixture <- function() {
  train <- make_variance_swd(
    data.frame(
      A = c(rep(0, 4), 1:4),
      B = c(1:4, rep(0, 4)),
      C = rep(1, 8),
      D = 1:8
    ),
    pa = c(rep(1, 4), rep(0, 4))
  )
  test <- make_variance_swd(
    data.frame(
      A = c(rep(0, 4), 1:4),
      B = c(1:4, rep(0, 4)),
      C = rep(1, 8),
      D = rep(0, 8)
    ),
    pa = c(rep(1, 4), rep(0, 4))
  )
  list(train = train, test = test)
}

make_variance_raster <- function() {
  testthat::skip_if_not_installed("terra")
  env <- terra::rast(
    nrows = 10, ncols = 10, nlyrs = 4,
    xmin = 0, xmax = 10, ymin = 0, ymax = 10
  )
  names(env) <- c("A", "B", "C", "D")
  terra::values(env) <- cbind(
    seq_len(100),
    rep(seq_len(10), each = 10),
    rep(1, 100),
    rev(seq_len(100))
  )
  env
}

variance_rows <- function(result, variable, check = NULL,
                          group = NULL, fold = NULL) {
  d <- result$diagnostics
  selected <- d$variable == variable
  if (!is.null(check)) selected <- selected & d$check == check
  if (!is.null(group)) selected <- selected & d$group == group
  if (!is.null(fold)) {
    selected <- selected & !is.na(d$fold) & d$fold == as.character(fold)
  }
  d[selected, , drop = FALSE]
}

test_that("varying predictors are retained and output schema is stable", {
  train <- make_variance_swd(
    data.frame(var1 = 1:6, var2 = 6:1),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  result <- screen_egv_variance(train, verbose = FALSE)

  expect_setequal(result$keep, c("var1", "var2"))
  expect_length(result$remove, 0L)
  expect_length(result$removed_by_check, 0L)
  expect_true(all(!result$diagnostics$failed))
  expect_true(all(result$diagnostics$used_for_selection))
  expect_named(result, c("keep", "remove", "removed_by_check", "diagnostics"))
  expect_true(all(c(
    "check", "subset", "data_source", "group", "fold",
    "used_for_selection", "variable", "n_finite", "sd",
    "failed", "reason"
  ) %in% names(result$diagnostics)))
  expect_setequal(unique(result$diagnostics$group),
                  c("presences", "background", "combined"))
})

test_that("constant predictors are removed and reported", {
  train <- make_variance_swd(
    data.frame(varying = 1:6, constant = rep(1, 6)),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  result <- screen_egv_variance(train, verbose = FALSE)
  expect_equal(result$keep, "varying")
  expect_equal(result$remove, "constant")
  expect_true(all(variance_rows(result, "constant", "training")$failed))
  expect_true(length(result$removed_by_check) > 0L)
})

test_that("zero variance among presences is detected", {
  train <- make_variance_swd(
    data.frame(good = 1:6,
               presence_constant = c(5, 5, 5, 1, 2, 3)),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  result <- screen_egv_variance(train, verbose = FALSE)
  d <- variance_rows(result, "presence_constant", "training", "presences")
  expect_equal(nrow(d), 1L)
  expect_true(d$failed[[1]])
  expect_equal(d$sd[[1]], 0)
  expect_true(d$used_for_selection[[1]])
  expect_true("presence_constant" %in% result$remove)
  expect_false(variance_rows(result, "presence_constant",
                             "training", "combined")$failed[[1]])
})

test_that("presence-only and background-only failures differ from combined", {
  f <- make_variance_fixture()
  out <- screen_egv_variance(f$train, verbose = FALSE)
  expect_true(variance_rows(out, "A", "training", "presences")$failed[[1]])
  expect_false(variance_rows(out, "A", "training", "combined")$failed[[1]])
  expect_true(variance_rows(out, "B", "training", "background")$failed[[1]])
  expect_false(variance_rows(out, "B", "training", "combined")$failed[[1]])
  expect_true(all(variance_rows(out, "C", "training")$failed))
  expect_setequal(out$remove, c("A", "B", "C"))
  expect_equal(out$keep, "D")
})

test_that("disabled groups are omitted, and diagnostics-only presences do not select", {
  f <- make_variance_fixture()
  combined <- screen_egv_variance(
    f$train,
    check_presences = FALSE, check_background = FALSE,
    verbose = FALSE
  )
  expect_equal(combined$keep, c("A", "B", "D"))
  expect_setequal(unique(combined$diagnostics$group), "combined")

  diagnostics_only <- screen_egv_variance(
    f$train,
    presence_affects_selection = FALSE, check_background = FALSE,
    verbose = FALSE
  )
  expect_equal(diagnostics_only$keep, c("A", "B", "D"))
  pres <- subset(diagnostics_only$diagnostics, group == "presences")
  expect_true(any(pres$failed))
  expect_false(any(pres$used_for_selection))
  expect_false("A" %in% diagnostics_only$remove)
  expect_false(any(vapply(diagnostics_only$removed_by_check,
                          function(x) "A" %in% x, logical(1))))
})

test_that("CV checks identify variance failure in a training presence partition", {
  train <- make_variance_swd(
    data.frame(
      good = c(1:6, 10:12),
      fold_constant = c(1, 1, 2, 2, 2, 2, 3, 4, 5)
    ),
    pa = c(rep(1, 6), rep(0, 3))
  )
  out <- screen_egv_variance(
    train, occ_folds = c(1, 1, 2, 2, 3, 3), verbose = FALSE
  )
  d <- variance_rows(out, "fold_constant", "cv_training", "presences", 1)
  expect_equal(nrow(d), 1L)
  expect_true(d$failed[[1]])
  expect_equal(d$sd[[1]], 0)
  expect_match(d$subset[[1]], "held_out_fold_1", fixed = TRUE)
  expect_true("fold_constant" %in% out$remove)
})

test_that("CV background can be shared or withheld by its own folds", {
  f <- make_variance_fixture()
  shared <- screen_egv_variance(
    f$train, occ_folds = c(1, 1, 2, 2),
    check_presences = FALSE, verbose = FALSE
  )
  partitioned <- screen_egv_variance(
    f$train, occ_folds = c(1, 1, 2, 2),
    bg_folds = c(1, 1, 2, 2),
    check_presences = FALSE, verbose = FALSE
  )
  shared_cv <- subset(shared$diagnostics, check == "cv_training")
  part_cv <- subset(partitioned$diagnostics, check == "cv_training")
  expect_true(all(grepl("shared_background", shared_cv$subset, fixed = TRUE)))
  expect_false(any(grepl("shared_background", part_cv$subset, fixed = TRUE)))
  expect_setequal(unique(part_cv$fold), c("1", "2"))
  combined <- subset(part_cv, group == "combined" & variable == "D")
  background <- subset(part_cv, group == "background" & variable == "D")
  expect_equal(combined$n_finite, c(4L, 4L))
  expect_equal(background$n_finite, c(2L, 2L))
})

test_that("pooled test backgrounds are an additional diagnostic check", {
  train <- make_variance_swd(
    data.frame(var1 = c(1, 2, 3, 10, 10, 10)),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  test <- make_variance_swd(
    data.frame(var1 = c(4, 5, 20, 20, 20)),
    pa = c(1, 1, 0, 0, 0)
  )
  without <- screen_egv_variance(
    train, test = test, include_test_bg = FALSE, verbose = FALSE
  )
  pooled <- screen_egv_variance(
    train, test = test, include_test_bg = TRUE, verbose = FALSE
  )
  expect_equal(pooled$keep, without$keep)
  expect_true("var1" %in% pooled$remove)
  base_bg <- variance_rows(pooled, "var1", "training", "background")
  expect_true(base_bg$failed[[1]])
  pooled_bg <- variance_rows(pooled, "var1", "pooled_background", "background")
  expect_equal(nrow(pooled_bg), 1L)
  expect_false(pooled_bg$failed[[1]])
  expect_false(pooled_bg$used_for_selection[[1]])
  expect_true(pooled_bg$sd[[1]] > 0)

  opt_in <- screen_egv_variance(
    train, test = test, include_test_bg = TRUE,
    test_affects_selection = TRUE, verbose = FALSE
  )
  expect_true(variance_rows(opt_in, "var1", "pooled_background",
                            "background")$used_for_selection[[1]])
  # A passing pooled check cannot cancel a failed training-only check.
  expect_true("var1" %in% opt_in$remove)
})

test_that("independent-test variance checks are diagnostic-only by default", {
  f <- make_variance_fixture()
  args <- list(check_presences = FALSE, check_background = FALSE,
               verbose = FALSE)
  base <- do.call(screen_egv_variance, c(list(train = f$train), args))
  out <- do.call(screen_egv_variance, c(
    list(train = f$train, test = f$test,
         check_test = TRUE, include_test_bg = TRUE), args
  ))
  expect_equal(out$keep, base$keep)
  test_rows <- subset(out$diagnostics, check == "independent_test")
  expect_true(nrow(test_rows) > 0L)
  expect_false(any(test_rows$used_for_selection))
  expect_false(any(subset(out$diagnostics,
                          check == "pooled_background")$used_for_selection))
  d <- variance_rows(out, "D", "independent_test", "combined")
  expect_equal(nrow(d), 1L)
  expect_true(d$failed[[1]])
  expect_false(d$used_for_selection[[1]])
  expect_true("D" %in% out$keep)
})

test_that("explicit independent-test opt-in can alter selection", {
  f <- make_variance_fixture()
  out <- screen_egv_variance(
    f$train, test = f$test,
    check_test = TRUE, test_affects_selection = TRUE,
    check_presences = FALSE, check_background = FALSE,
    verbose = FALSE
  )
  d <- variance_rows(out, "D", "independent_test", "combined")
  expect_true(d$failed[[1]])
  expect_true(d$used_for_selection[[1]])
  expect_true("D" %in% out$remove)
  expect_false("D" %in% out$keep)
})

test_that("independent-test records are not checked by default", {
  f <- make_variance_fixture()
  without_test <- screen_egv_variance(f$train, verbose = FALSE)
  with_test <- screen_egv_variance(
    f$train, test = f$test, verbose = FALSE
  )
  expect_identical(with_test$keep, without_test$keep)
  expect_identical(with_test$remove, without_test$remove)
  expect_false(any(with_test$diagnostics$check == "independent_test"))
  expect_false(any(with_test$diagnostics$check == "pooled_background"))
})

test_that("independent test can be diagnosed in all three groups", {
  f <- make_variance_fixture()
  out <- screen_egv_variance(
    f$train, test = f$test, check_test = TRUE, verbose = FALSE
  )
  test_rows <- subset(out$diagnostics, check == "independent_test")
  expect_setequal(unique(test_rows$group),
                  c("presences", "background", "combined"))
  expect_false(any(test_rows$used_for_selection))
  expect_true(all(variance_rows(out, "D", "independent_test")$failed))
  expect_true("D" %in% out$keep)
})

test_that("sd_tol screens small but nonzero standard deviations", {
  train <- make_variance_swd(
    data.frame(
      good = 1:6,
      tiny = 1 + (0:5) * 1e-10
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  exact <- screen_egv_variance(train, sd_tol = 0, verbose = FALSE)
  tolerant <- screen_egv_variance(train, sd_tol = 1e-8, verbose = FALSE)
  expect_true("tiny" %in% exact$keep)
  expect_true("tiny" %in% tolerant$remove)
})

test_that("fewer than two finite observations are reported", {
  train <- make_variance_swd(
    data.frame(
      good = 1:6,
      sparse = c(NA_real_, Inf, 2, 3, 4, 5)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  out <- screen_egv_variance(train, verbose = FALSE)
  d <- variance_rows(out, "sparse", "training", "presences")
  expect_equal(d$n_finite[[1]], 1L)
  expect_true(is.na(d$sd[[1]]))
  expect_true(d$failed[[1]])
  expect_equal(d$reason[[1]], "fewer_than_2_finite_values")
  expect_true("sparse" %in% out$remove)
})

test_that("raster sampling is separate from supplied training backgrounds", {
  f <- make_variance_fixture()
  env <- make_variance_raster()
  out <- screen_egv_variance(
    f$train, env = env, env_sample_n = 30, seed = 2, verbose = FALSE
  )
  expect_true(any(out$diagnostics$check == "environment"))
  expect_true(any(out$diagnostics$check == "training" &
                    out$diagnostics$group == "background"))
  env_c <- variance_rows(out, "C", "environment", "environment")
  expect_equal(nrow(env_c), 1L)
  expect_true(env_c$failed[[1]])
  expect_equal(env_c$sd[[1]], 0)
  expect_true(env_c$used_for_selection[[1]])
})

test_that("environment-only zero variance leads to removal", {
  testthat::skip_if_not_installed("terra")
  train <- make_variance_swd(
    data.frame(varying = 1:6, env_constant = 2:7),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  env <- terra::rast(nrows = 10, ncols = 10, nlyrs = 2,
                     xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  names(env) <- c("varying", "env_constant")
  terra::values(env) <- cbind(seq_len(100), rep(5, 100))
  without <- screen_egv_variance(train, verbose = FALSE)
  with <- screen_egv_variance(
    train, env = env, env_sample_n = 50, seed = 1,
    verbose = FALSE
  )
  expect_true("env_constant" %in% without$keep)
  expect_true("env_constant" %in% with$remove)
  expect_true(variance_rows(with, "env_constant", "environment")$failed[[1]])
})

test_that("environmental sample with seed is reproducible and restores RNG", {
  f <- make_variance_fixture()
  env <- make_variance_raster()
  local({
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(23)
    before <- get(".Random.seed", envir = .GlobalEnv)
    out1 <- screen_egv_variance(
      f$train, env = env, env_sample_n = 30, seed = 2,
      verbose = FALSE
    )
    expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
    out2 <- screen_egv_variance(
      f$train, env = env, env_sample_n = 30, seed = 2,
      verbose = FALSE
    )
    expect_identical(out1$diagnostics, out2$diagnostics)
    expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  })
})

test_that("numeric arguments reject fractional and nonfinite values", {
  f <- make_variance_fixture()
  run <- function(...) screen_egv_variance(f$train, verbose = FALSE, ...)
  expect_error(run(env_sample_n = 100.9), "env_sample_n")
  expect_error(run(env_sample_n = 1L), "env_sample_n")
  expect_error(run(env_sample_n = Inf), "env_sample_n")
  expect_error(run(env_sample_n = .Machine$integer.max + 1), "env_sample_n")
  expect_error(run(seed = 1.1), "seed")
  expect_error(run(seed = "1"), "seed")
  expect_error(run(seed = -1), "seed")
  expect_error(run(seed = NA_real_), "seed")
  expect_error(run(sd_tol = -1), "sd_tol")
  expect_error(run(sd_tol = Inf), "sd_tol")
  expect_error(run(sd_tol = NA_real_), "sd_tol")
  expect_error(run(sd_tol = "0"), "sd_tol")
  expect_silent(run(env_sample_n = 100, seed = 1))
})

test_that("logical arguments require scalar TRUE or FALSE", {
  f <- make_variance_fixture()
  run <- function(...) screen_egv_variance(f$train, verbose = FALSE, ...)
  expect_error(run(include_test_bg = 1), "include_test_bg")
  expect_error(run(include_test_bg = NA), "include_test_bg")
  expect_error(run(check_test = "TRUE"), "check_test")
  expect_error(run(check_combined = NA), "check_combined")
  expect_error(run(check_background = c(TRUE, FALSE)), "check_background")
  expect_error(run(check_presences = NULL), "check_presences")
  expect_error(run(presence_affects_selection = 1), "presence_affects_selection")
  expect_error(run(test_affects_selection = "FALSE"), "test_affects_selection")
  # Important: test 'verbose' directly, without passing a duplicate argument.
  expect_error(screen_egv_variance(f$train, verbose = NA), "verbose")
  expect_error(run(check_presences = FALSE, check_background = FALSE,
                   check_combined = FALSE), "At least one")
})

test_that("missing test and invalid SWD inputs report useful errors", {
  f <- make_variance_fixture()
  expect_error(screen_egv_variance(data.frame(x = 1:5)), "train")
  expect_error(screen_egv_variance(f$train, test = data.frame(x = 1:5)), "test")
  expect_error(screen_egv_variance(f$train, check_test = TRUE), "requires `test`")
  expect_error(screen_egv_variance(f$train, include_test_bg = TRUE),
               "requires `test`")
  mismatched <- make_variance_swd(
    data.frame(other = 1:6), pa = c(1, 1, 1, 0, 0, 0)
  )
  expect_error(screen_egv_variance(f$train, test = mismatched),
               "missing from `test@data`")
  non_numeric <- make_variance_swd(
    data.frame(num = 1:6, category = letters[1:6]),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  expect_error(screen_egv_variance(non_numeric), "Non-numeric")
})

test_that("CV fold arguments require matching non-missing assignments", {
  f <- make_variance_fixture()
  run <- function(...) screen_egv_variance(f$train, verbose = FALSE, ...)
  expect_error(run(occ_folds = c(1, 2)), "occ_folds")
  expect_error(run(occ_folds = c(1, NA, 2, 2)), "occ_folds")
  expect_error(run(occ_folds = rep(1, 4)), "at least two folds")
  expect_error(run(bg_folds = c(1, 1, 2, 2)), "requires `occ_folds`")
  expect_error(run(occ_folds = c(1, 1, 2, 2),
                   bg_folds = c(1, NA, 2, 2)), "bg_folds")
  expect_error(run(occ_folds = c(1, 1, 2, 2),
                   bg_folds = c(1, 2)), "bg_folds")
  expect_error(run(occ_folds = c(1, 1, 2, 2),
                   bg_folds = c(1, 1, 3, 3)), "bg_folds")
})


# ---- Additional edge cases -----------------------------------------------

test_that("variance screening preserves training predictor order", {
  train <- make_variance_swd(
    data.frame(z = 1:6, a = rep(2, 6), m = 6:1),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  out <- screen_egv_variance(train, verbose = FALSE)
  expect_identical(out$keep, c("z", "m"))
  expect_identical(out$remove, "a")
})

test_that("background-only checks can be disabled without losing diagnostics for other groups", {
  f <- make_variance_fixture()
  out <- screen_egv_variance(
    f$train, check_background = FALSE,
    presence_affects_selection = FALSE, verbose = FALSE
  )
  expect_false(any(out$diagnostics$group == "background"))
  expect_setequal(out$remove, "C")
  expect_setequal(out$keep, c("A", "B", "D"))
})

test_that("a failing test-only presence check is diagnostic unless opted in", {
  train <- make_variance_swd(
    data.frame(x = c(1, 2, 3, 4, 5, 6)),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  test <- make_variance_swd(
    data.frame(x = c(10, 10, 10, 1, 2, 3)),
    pa = c(1, 1, 1, 0, 0, 0)
  )
  diagnostic <- screen_egv_variance(
    train, test = test, check_test = TRUE,
    verbose = FALSE
  )
  d <- variance_rows(diagnostic, "x", "independent_test", "presences")
  expect_equal(nrow(d), 1L)
  expect_true(d$failed[[1]])
  expect_false(d$used_for_selection[[1]])
  expect_true("x" %in% diagnostic$keep)

  selected <- screen_egv_variance(
    train, test = test, check_test = TRUE,
    test_affects_selection = TRUE, verbose = FALSE
  )
  expect_true("x" %in% selected$remove)
  expect_true(variance_rows(selected, "x", "independent_test",
                            "presences")$used_for_selection[[1]])
})

test_that("CV empty training presence partitions are diagnosed rather than crashing", {
  train <- make_variance_swd(
    data.frame(x = c(1, 2, 3, 4, 5)),
    pa = c(1, 1, 0, 0, 0)
  )
  out <- screen_egv_variance(
    train, occ_folds = c("one", "two"), verbose = FALSE
  )
  cv <- subset(out$diagnostics,
               check == "cv_training" & group == "presences")
  expect_equal(nrow(cv), 2L)
  expect_true(all(cv$failed))
  expect_true(all(cv$n_finite == 1L))
  expect_true(all(cv$reason == "fewer_than_2_finite_values"))
})

# ---- Package-data integration tests ---------------------------------------
# The test below confirms that the environmental raster distributed with
# sdmhelpers can be used without relying on private absolute file paths.

test_that("bundled example EGV raster works in screening", {
  testthat::skip_if_not_installed("terra")
  tif <- system.file("extdata", "example_egv.tif", package = "sdmhelpers")
  if (!nzchar(tif)) testthat::skip("example_egv.tif not installed in package")

  env <- terra::rast(tif)
  expect_equal(terra::nlyr(env), 5L)
  expect_setequal(names(env),
                  c("egv_280", "egv_293", "egv_302", "egv_385", "egv_400"))

  # The random raster sample is independent of the SWD background records.
  # Use sample values for a realistic SWD fixture while retaining deterministic
  # unit tests above for exact failure patterns.
  sample <- terra::spatSample(
    env, size = 80, method = "random", na.rm = TRUE,
    values = TRUE, xy = FALSE, as.df = TRUE
  )
  sample <- as.data.frame(sample)
  sample <- sample[stats::complete.cases(sample), names(env), drop = FALSE]
  expect_gte(nrow(sample), 6L)
  n <- nrow(sample)
  n_pres <- floor(n / 2L)
  train <- make_variance_swd(
    sample, pa = c(rep(1, n_pres), rep(0, n - n_pres))
  )

  out <- screen_egv_variance(
    train = train, env = env, env_sample_n = 100L,
    seed = 1L, verbose = FALSE
  )
  expect_setequal(c(out$keep, out$remove), names(env))
  expect_equal(nrow(subset(out$diagnostics, check == "environment")),
               terra::nlyr(env))
  expect_true(any(out$diagnostics$check == "training" &
                    out$diagnostics$group == "background"))
})

test_that("bundled presence and background coordinates support SWD screening", {
  testthat::skip_if_not_installed("terra")
  tif <- system.file("extdata", "example_egv.tif", package = "sdmhelpers")
  if (!nzchar(tif)) testthat::skip("example EGV raster not installed")

  e <- new.env(parent = emptyenv())
  utils::data(list = c("example_presences", "example_background"),
              package = "sdmhelpers", envir = e)
  if (!exists("example_presences", envir = e, inherits = FALSE) ||
      !exists("example_background", envir = e, inherits = FALSE)) {
    testthat::skip("example location datasets not installed")
  }
  p <- e$example_presences
  b <- e$example_background
  expect_true(all(c("x", "y") %in% names(p)))
  expect_true(all(c("x", "y") %in% names(b)))

  env <- terra::rast(tif)
  sample_at <- function(locations) {
    extracted <- terra::extract(env, as.matrix(locations[, c("x", "y")]))
    extracted <- as.data.frame(extracted)[, names(env), drop = FALSE]
    extracted[stats::complete.cases(extracted), , drop = FALSE]
  }
  pres_vals <- sample_at(p)
  bg_vals <- sample_at(b)
  if (nrow(pres_vals) < 3L || nrow(bg_vals) < 3L) {
    testthat::skip("too few example points overlap complete EGV cells")
  }
  pres_vals <- utils::head(pres_vals, 20L)
  bg_vals <- utils::head(bg_vals, 20L)
  train <- make_variance_swd(
    rbind(pres_vals, bg_vals),
    pa = c(rep(1, nrow(pres_vals)), rep(0, nrow(bg_vals)))
  )
  out <- screen_egv_variance(train, verbose = FALSE)
  expect_setequal(c(out$keep, out$remove), names(env))
  expect_setequal(unique(out$diagnostics$group),
                  c("presences", "background", "combined"))
  expect_equal(nrow(out$diagnostics), 3L * terra::nlyr(env))
})
