test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("phenology_overlap_weights returns expected structure", {

  target <- data.frame(
    CODE = rep(c("A", "B"), each = 20),
    DoY = c(
      seq(100, 138, length.out = 20),
      seq(200, 238, length.out = 20)
    )
  )

  group <- data.frame(
    CODE = rep(c("A", "C", "D"), each = 20),
    DoY = c(
      seq(100, 138, length.out = 20),
      seq(120, 158, length.out = 20),
      seq(250, 288, length.out = 20)
    )
  )

  res <- phenology_overlap_weights(
    target = target,
    group = group,
    target_code = "A"
  )

  expect_s3_class(res, "data.frame")
  expect_named(res, c("CODE", "overlap", "weight"))

  expect_equal(nrow(res), 3)
  expect_equal(res$CODE, c("A", "C", "D"))

  expect_type(res$overlap, "double")
  expect_type(res$weight, "double")
})

test_that("target species is retained when present in comparison group", {

  target <- data.frame(
    CODE = rep("A", 20),
    DoY = seq(100, 138, length.out = 20)
  )

  group <- data.frame(
    CODE = rep(c("A", "B"), each = 20),
    DoY = c(
      seq(100, 138, length.out = 20),
      seq(180, 218, length.out = 20)
    )
  )

  res <- phenology_overlap_weights(
    target = target,
    group = group,
    target_code = "A"
  )

  expect_true("A" %in% res$CODE)
  expect_equal(nrow(res), 2)
})

test_that("weights are normalized to sum to one", {

  target <- data.frame(
    CODE = rep("A", 30),
    DoY = seq(100, 160, length.out = 30)
  )

  group <- data.frame(
    CODE = rep(c("A", "B", "C"), each = 30),
    DoY = c(
      seq(100, 160, length.out = 30),
      seq(120, 180, length.out = 30),
      seq(180, 240, length.out = 30)
    )
  )

  res <- phenology_overlap_weights(
    target = target,
    group = group,
    target_code = "A"
  )

  expect_equal(
    sum(res$weight, na.rm = TRUE),
    1,
    tolerance = 1e-8
  )
})

test_that("identical distributions have maximal overlap", {

  doy <- seq(100, 160, length.out = 30)

  target <- data.frame(
    CODE = rep("A", length(doy)),
    DoY = doy
  )

  group <- data.frame(
    CODE = rep(c("A", "B"), each = length(doy)),
    DoY = c(
      doy,
      seq(220, 280, length.out = length(doy))
    )
  )

  res <- phenology_overlap_weights(
    target = target,
    group = group,
    target_code = "A"
  )

  overlap_A <- res$overlap[res$CODE == "A"]
  overlap_B <- res$overlap[res$CODE == "B"]

  expect_true(overlap_A > overlap_B)
})

test_that("non-finite DoY values are removed when na_rm is TRUE", {

  target <- data.frame(
    CODE = rep("A", 12),
    DoY = c(100:109, NA, Inf)
  )

  group <- data.frame(
    CODE = rep("B", 12),
    DoY = c(105:114, NA, Inf)
  )

  expect_no_error(
    phenology_overlap_weights(
      target = target,
      group = group,
      target_code = "A",
      na_rm = TRUE
    )
  )
})

test_that("non-finite target values cause an error when na_rm is FALSE", {

  target <- data.frame(
    CODE = rep("A", 10),
    DoY = c(100:108, NA)
  )

  group <- data.frame(
    CODE = rep("B", 10),
    DoY = 110:119
  )

  expect_error(
    phenology_overlap_weights(
      target = target,
      group = group,
      target_code = "A",
      na_rm = FALSE
    ),
    "Non-finite values found"
  )
})

test_that("comparison groups with too few observations receive NA", {

  target <- data.frame(
    CODE = rep("A", 20),
    DoY = 100:119
  )

  group <- data.frame(
    CODE = c(rep("B", 20), "C"),
    DoY = c(110:129, 150)
  )

  res <- phenology_overlap_weights(
    target = target,
    group = group,
    target_code = "A",
    min_n = 2
  )

  expect_true(is.na(res$overlap[res$CODE == "C"]))
  expect_true(is.na(res$weight[res$CODE == "C"]))

  expect_true(is.finite(res$overlap[res$CODE == "B"]))
})


test_that("invalid inputs generate informative errors", {

  target <- data.frame(
    CODE = rep("A", 10),
    DoY = 100:109
  )

  group <- data.frame(
    CODE = rep("B", 10),
    DoY = 110:119
  )

  expect_error(
    phenology_overlap_weights(
      target = target,
      group = group,
      target_code = "X"
    ),
    "not found"
  )

  expect_error(
    phenology_overlap_weights(
      target = target,
      group = group,
      target_code = "A",
      target_doy = "wrong"
    ),
    "not found"
  )

  expect_error(
    phenology_overlap_weights(
      target = target,
      group = group,
      target_code = "A",
      group_id = "wrong"
    ),
    "not found"
  )

  expect_error(
    phenology_overlap_weights(
      target = target,
      group = group,
      target_code = "A",
      min_n = 1
    ),
    "must be an integer >= 2"
  )
})


test_that("custom column names are supported", {

  target <- data.frame(
    species = rep("AAA", 20),
    doy = 100:119
  )

  group <- data.frame(
    species = rep(c("AAA", "BBB"), each = 20),
    doy = c(100:119, 130:149)
  )

  res <- phenology_overlap_weights(
    target = target,
    group = group,
    target_code = "AAA",
    target_id = "species",
    group_id = "species",
    target_doy = "doy",
    group_doy = "doy"
  )

  expect_named(res, c("species", "overlap", "weight"))
  expect_equal(res$species, c("AAA", "BBB"))
})

