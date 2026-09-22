make_test_swd <- function(data, pa, species = "test") {

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


test_that("predictors with variance in all basic subsets are retained", {

  train <- make_test_swd(
    data = data.frame(
      var1 = c(1, 2, 3, 4, 5, 6),
      var2 = c(6, 5, 4, 3, 2, 1)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  result <- screen_egv_variance(
    train = train,
    verbose = FALSE
  )

  expect_setequal(
    result$keep,
    c("var1", "var2")
  )

  expect_length(result$remove, 0)

  expect_true(
    all(!result$diagnostics$failed)
  )
})


test_that("zero-variance predictors are removed", {

  train <- make_test_swd(
    data = data.frame(
      varying = c(1, 2, 3, 4, 5, 6),
      constant = rep(1, 6)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  result <- screen_egv_variance(
    train = train,
    verbose = FALSE
  )

  expect_equal(
    result$keep,
    "varying"
  )

  expect_equal(
    result$remove,
    "constant"
  )

  expect_true(
    any(
      result$diagnostics$variable == "constant" &
        result$diagnostics$failed
    )
  )
})


test_that("zero variance among presences is detected", {

  train <- make_test_swd(
    data = data.frame(
      good = c(1, 2, 3, 4, 5, 6),
      presence_constant = c(5, 5, 5, 1, 2, 3)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  result <- screen_egv_variance(
    train = train,
    verbose = FALSE
  )

  expect_equal(
    result$remove,
    "presence_constant"
  )

  diagnostic <- result$diagnostics[
    result$diagnostics$variable == "presence_constant" &
      result$diagnostics$check == "training" &
      result$diagnostics$subset == "presences",
    ,
    drop = FALSE
  ]

  expect_equal(nrow(diagnostic), 1)
  expect_true(diagnostic$failed)
  expect_equal(diagnostic$sd, 0)
})


test_that("zero variance in a CV training partition is detected", {

  train <- make_test_swd(
    data = data.frame(
      good = c(
        1, 2, 3, 4, 5, 6,
        10, 11, 12
      ),
      fold_constant = c(
        1, 1, 2, 2, 2, 2,
        3, 4, 5
      )
    ),
    pa = c(
      rep(1, 6),
      rep(0, 3)
    )
  )

  folds <- c(
    1, 1,
    2, 2,
    3, 3
  )

  result <- screen_egv_variance(
    train = train,
    occ_folds = folds,
    verbose = FALSE
  )

  expect_true(
    "fold_constant" %in% result$remove
  )

  diagnostic <- result$diagnostics[
    result$diagnostics$variable == "fold_constant" &
      result$diagnostics$check == "cv_training_presences" &
      result$diagnostics$subset == "held_out_fold_1",
    ,
    drop = FALSE
  ]

  expect_equal(nrow(diagnostic), 1)
  expect_true(diagnostic$failed)
  expect_equal(diagnostic$sd, 0)
})


test_that("independent-test background can be included in pooled background check", {

  train <- make_test_swd(
    data = data.frame(
      var1 = c(
        1, 2, 3,  # presences
        10, 10, 10 # training background
      )
    ),
    pa = c(
      1, 1, 1,
      0, 0, 0
    )
  )

  test <- make_test_swd(
    data = data.frame(
      var1 = c(
        4, 5,       # test presences
        20, 20, 20  # test background
      )
    ),
    pa = c(
      1, 1,
      0, 0, 0
    )
  )

  result_without_test_bg <- screen_egv_variance(
    train = train,
    test = test,
    include_test_bg = FALSE,
    verbose = FALSE
  )

  result_with_test_bg <- screen_egv_variance(
    train = train,
    test = test,
    include_test_bg = TRUE,
    verbose = FALSE
  )

  expect_true(
    "var1" %in% result_without_test_bg$remove
  )

  expect_true(
    "var1" %in% result_with_test_bg$keep
  )
})

test_that("zero variance in independent testing data is detected when requested", {

  train <- make_test_swd(
    data = data.frame(
      var1 = c(1, 2, 3, 4, 5, 6)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  test <- make_test_swd(
    data = data.frame(
      var1 = rep(10, 6)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  result_no_check <- screen_egv_variance(
    train = train,
    test = test,
    include_test_bg = FALSE,
    check_test = FALSE,
    verbose = FALSE
  )

  result_check <- screen_egv_variance(
    train = train,
    test = test,
    include_test_bg = FALSE,
    check_test = TRUE,
    verbose = FALSE
  )

  expect_true(
    "var1" %in% result_no_check$keep
  )

  expect_true(
    "var1" %in% result_check$remove
  )

  expect_true(
    any(
      result_check$diagnostics$check == "independent_test" &
        result_check$diagnostics$variable == "var1" &
        result_check$diagnostics$failed
    )
  )
})


test_that("sd_tol flags predictors with very small variance", {

  train <- make_test_swd(
    data = data.frame(
      good = c(1, 2, 3, 4, 5, 6),
      tiny = c(
        1,
        1 + 1e-10,
        1 + 2e-10,
        1 + 3e-10,
        1 + 4e-10,
        1 + 5e-10
      )
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  exact <- screen_egv_variance(
    train = train,
    sd_tol = 0,
    verbose = FALSE
  )

  tolerant <- screen_egv_variance(
    train = train,
    sd_tol = 1e-8,
    verbose = FALSE
  )

  expect_true(
    "tiny" %in% exact$keep
  )

  expect_true(
    "tiny" %in% tolerant$remove
  )
})

test_that("zero variance in environmental domain is detected", {

  train <- make_test_swd(
    data = data.frame(
      good = c(1, 2, 3, 4, 5, 6),
      env_constant = c(2, 3, 4, 5, 6, 7)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  env <- terra::rast(
    nrows = 10,
    ncols = 10,
    xmin = 0,
    xmax = 10,
    ymin = 0,
    ymax = 10,
    nlyrs = 2
  )

  names(env) <- c(
    "good",
    "env_constant"
  )

  terra::values(env) <- cbind(
    seq_len(terra::ncell(env)),
    rep(5, terra::ncell(env))
  )

  result <- screen_egv_variance(
    train = train,
    env = env,
    env_sample_n = 50,
    seed = 1,
    verbose = FALSE
  )

  expect_true(
    "good" %in% result$keep
  )

  expect_true(
    "env_constant" %in% result$remove
  )

  expect_true(
    any(
      result$diagnostics$variable == "env_constant" &
        result$diagnostics$check == "environment" &
        result$diagnostics$failed
    )
  )
})

test_that("invalid inputs produce informative errors", {

  train <- make_test_swd(
    data = data.frame(
      var1 = c(1, 2, 3, 4, 5, 6)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  expect_error(
    screen_egv_variance(
      train = data.frame(x = 1:5)
    ),
    "`train` must be an SDMtune `SWD` object"
  )

  expect_error(
    screen_egv_variance(
      train = train,
      sd_tol = -1
    ),
    "`sd_tol` must be one finite, non-negative number"
  )

  expect_error(
    screen_egv_variance(
      train = train,
      occ_folds = c(1, 2)
    ),
    "`occ_folds` must have one value for each presence record"
  )

  expect_error(
    screen_egv_variance(
      train = train,
      check_test = TRUE
    ),
    "`check_test = TRUE` requires a `test` SWD object"
  )
})

test_that("occ_folds cannot contain missing values", {

  train <- make_test_swd(
    data = data.frame(
      var1 = c(1, 2, 3, 4, 5, 6)
    ),
    pa = c(1, 1, 1, 0, 0, 0)
  )

  expect_error(
    screen_egv_variance(
      train = train,
      occ_folds = c(1, NA, 2)
    ),
    "`occ_folds` must not contain missing values"
  )
})

