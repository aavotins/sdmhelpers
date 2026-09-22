#' Screen Environmental Predictors for Zero Variance
#'
#' Identifies environmental predictors that have zero or insufficient variance
#' in data subsets relevant to species distribution model training,
#' cross-validation, background characterization, and optionally independent
#' testing.
#'
#' @param train An `SWD` object, typically created with
#'   [SDMtune::prepareSWD()], containing the training presence and
#'   absence/background records.
#'
#' @param test An optional `SWD` object containing independent testing data.
#'   When supplied, its background records can be included in the pooled
#'   background variance check. The complete independent testing dataset can
#'   additionally be checked by setting `check_test = TRUE`.
#'
#' @param occ_folds An optional vector giving the cross-validation fold assigned
#'   to each presence record in `train`. Its length must equal the number of
#'   presence records in `train`. For each fold, predictor variance is evaluated
#'   using the presence records belonging to all other folds, corresponding to
#'   the training partition when that fold is held out.
#'
#' @param env An optional [terra::SpatRaster] containing the environmental
#'   predictors. When supplied, `env_sample_n` random raster cells are sampled
#'   and predictor variance is evaluated across the sampled environmental
#'   domain.
#'
#' @param env_sample_n Positive integer. Number of cells randomly sampled from
#'   `env` when `env` is supplied. Defaults to `10000`.
#'
#' @param seed Optional integer used for reproducible sampling of `env`.
#'   The previous state of R's random-number generator is restored before the
#'   function exits.
#'
#' @param sd_tol Non-negative numeric value defining the minimum acceptable
#'   standard deviation. A predictor is flagged when its standard deviation is
#'   less than or equal to `sd_tol`. The default, `0`, removes only predictors
#'   with exactly zero variance, together with predictors for which variance
#'   cannot be estimated.
#'
#' @param include_test_bg Logical. If `TRUE` and `test` is supplied,
#'   background records from the independent testing dataset are pooled with
#'   training background records for the background variance check. Defaults
#'   to `TRUE`.
#'
#' @param check_test Logical. If `TRUE`, predictor variance is additionally
#'   evaluated across all records in the independent testing `SWD` object.
#'   Defaults to `FALSE`.
#'
#' @param verbose Logical. If `TRUE`, print a short summary of the screening
#'   results. Defaults to `TRUE`.
#'
#' @details
#' The function is intended as a preprocessing check for environmental
#' predictors used in species distribution modelling. A predictor is retained
#' only when it passes every variance check that is requested.
#'
#' By default, the following subsets are evaluated:
#'
#' \itemize{
#'   \item all presence and absence/background records in the training data;
#'   \item all training presence records;
#'   \item pooled training background records and, when `include_test_bg = TRUE`,
#'     independent testing background records;
#'   \item each cross-validation presence training partition when
#'     `occ_folds` is supplied; and
#'   \item a random sample of the environmental domain when `env` is supplied.
#' }
#'
#' For the cross-validation check, the function does not evaluate the records
#' within a held-out fold alone. For each fold `k`, it evaluates the records for
#' which `occ_folds != k`. These are the presence records available for model
#' fitting when fold `k` is used for validation.
#'
#' If `check_test = TRUE`, all records in `test` are also evaluated as a
#' separate subset.
#'
#' Standard deviations are calculated using finite values only. Predictors
#' having fewer than two finite observations in a checked subset are flagged
#' because their variance cannot be estimated.
#'
#' Predictor names are taken from `train@data`. The same predictors must be
#' present in `test@data`, when `test` is supplied, and in `env`, when `env`
#' is supplied.
#'
#' The environmental raster itself is not modified. Instead, the function
#' returns vectors containing the predictors that passed and failed screening.
#' These can subsequently be used to subset an environmental-variable table,
#' an `SWD` object, or a [terra::SpatRaster].
#'
#' @return A named list with the following components:
#'
#' \describe{
#'   \item{keep}{
#'     Character vector containing predictors that passed all requested
#'     variance checks.
#'   }
#'   \item{remove}{
#'     Character vector containing predictors that failed at least one
#'     variance check.
#'   }
#'   \item{removed_by_check}{
#'     Named list showing which predictors were flagged by each individual
#'     check.
#'   }
#'   \item{diagnostics}{
#'     A `data.frame` containing the check, subset, predictor name, number of
#'     finite observations, standard deviation, failure status, and reason for
#'     every predictor evaluated.
#'   }
#' }
#'
#' @importClassesFrom SDMtune SWD
#'
#' @seealso
#' [SDMtune::prepareSWD()],
#' [terra::spatSample()]
#'
#' @examples
#' \dontrun{
#' result <- screen_egv_variance(
#'   train = trenin_dati,
#'   test = testa_dati,
#'   occ_folds = block_folds$occs.grp,
#'   env = vide,
#'   env_sample_n = 10000,
#'   seed = 1
#' )
#'
#' result$keep
#' result$remove
#' result$removed_by_check
#'
#' # Retain only predictors passing all checks
#' videi <- videi[
#'   videi$layername %in% result$keep,
#' ]
#'
#' vide <- terra::rast(
#'   paste0(
#'     "./RasterGrids_100m/2024/Scaled/",
#'     videi$new_filename
#'   )
#' )
#' }
#'
#' @export
screen_egv_variance <- function(
    train,
    test = NULL,
    occ_folds = NULL,
    env = NULL,
    env_sample_n = 10000L,
    seed = NULL,
    sd_tol = 0,
    include_test_bg = TRUE,
    check_test = FALSE,
    verbose = TRUE
) {

  # ---- Validate inputs ----

  if (!methods::is(train, "SWD")) {
    stop("`train` must be an SDMtune `SWD` object.")
  }

  if (!is.null(test) && !methods::is(test, "SWD")) {
    stop("`test` must be an SDMtune `SWD` object or NULL.")
  }

  if (length(sd_tol) != 1L ||
      !is.numeric(sd_tol) ||
      !is.finite(sd_tol) ||
      sd_tol < 0) {
    stop("`sd_tol` must be one finite, non-negative number.")
  }

  if (length(env_sample_n) != 1L ||
      !is.numeric(env_sample_n) ||
      !is.finite(env_sample_n) ||
      env_sample_n < 2) {
    stop("`env_sample_n` must be an integer >= 2.")
  }

  env_sample_n <- as.integer(env_sample_n)

  # ---- Extract training data ----

  train_data <- as.data.frame(methods::slot(train, "data"))
  train_pa <- methods::slot(train, "pa")

  if (nrow(train_data) != length(train_pa)) {
    stop("The dimensions of `train@data` and `train@pa` do not agree.")
  }

  vars <- names(train_data)

  if (length(vars) == 0L) {
    stop("No environmental predictors were found in `train@data`.")
  }

  if (anyDuplicated(vars)) {
    stop("Predictor names in `train@data` must be unique.")
  }

  is_numeric <- vapply(train_data, is.numeric, logical(1))

  if (!all(is_numeric)) {
    stop(
      "All predictors must be numeric. Non-numeric predictors: ",
      paste(vars[!is_numeric], collapse = ", ")
    )
  }

  # ---- Extract testing data ----

  if (!is.null(test)) {

    test_data <- as.data.frame(methods::slot(test, "data"))
    test_pa <- methods::slot(test, "pa")

    if (nrow(test_data) != length(test_pa)) {
      stop("The dimensions of `test@data` and `test@pa` do not agree.")
    }

    missing_test <- setdiff(vars, names(test_data))

    if (length(missing_test) > 0L) {
      stop(
        "The following training predictors are missing from `test@data`: ",
        paste(missing_test, collapse = ", ")
      )
    }

    test_data <- test_data[, vars, drop = FALSE]
  }

  # ---- Helper for one variance check ----

  check_variance <- function(x, check, subset) {

    x <- as.data.frame(x)[, vars, drop = FALSE]

    n_finite <- vapply(
      x,
      function(z) sum(is.finite(z)),
      integer(1)
    )

    sds <- vapply(
      x,
      function(z) {

        z <- z[is.finite(z)]

        if (length(z) < 2L) {
          return(NA_real_)
        }

        stats::sd(z)
      },
      numeric(1)
    )

    failed <- n_finite < 2L |
      !is.finite(sds) |
      sds <= sd_tol

    reason <- ifelse(
      n_finite < 2L,
      "fewer_than_2_finite_values",
      ifelse(
        !is.finite(sds),
        "undefined_sd",
        ifelse(
          sds <= sd_tol,
          "sd_below_or_equal_to_tolerance",
          "passed"
        )
      )
    )

    data.frame(
      check = check,
      subset = subset,
      variable = vars,
      n_finite = n_finite,
      sd = sds,
      failed = failed,
      reason = reason,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- list()

  # ---- Complete training data ----

  diagnostics[[length(diagnostics) + 1L]] <- check_variance(
    train_data,
    check = "training",
    subset = "all"
  )

  # ---- Training presences ----

  pres_idx <- train_pa == 1

  if (!any(pres_idx)) {
    stop("No presence records were found in `train`.")
  }

  pres_data <- train_data[pres_idx, , drop = FALSE]

  diagnostics[[length(diagnostics) + 1L]] <- check_variance(
    pres_data,
    check = "training",
    subset = "presences"
  )

  # ---- Background data ----

  bg_idx <- train_pa == 0

  bg_data <- train_data[bg_idx, , drop = FALSE]
  bg_subset <- "training"

  if (!is.null(test) && isTRUE(include_test_bg)) {

    test_bg <- test_data[test_pa == 0, , drop = FALSE]

    if (nrow(test_bg) > 0L) {
      bg_data <- rbind(bg_data, test_bg)
      bg_subset <- "training + independent_test"
    }
  }

  if (nrow(bg_data) == 0L) {
    stop("No absence/background records were available for checking.")
  }

  diagnostics[[length(diagnostics) + 1L]] <- check_variance(
    bg_data,
    check = "background",
    subset = bg_subset
  )

  # ---- Independent testing data ----

  if (isTRUE(check_test)) {

    if (is.null(test)) {
      stop("`check_test = TRUE` requires a `test` SWD object.")
    }

    diagnostics[[length(diagnostics) + 1L]] <- check_variance(
      test_data,
      check = "independent_test",
      subset = "all"
    )
  }

  # ---- Cross-validation training partitions ----

  if (!is.null(occ_folds)) {

    if (length(occ_folds) != nrow(pres_data)) {
      stop(
        "`occ_folds` must have one value for each presence record in `train`."
      )
    }

    if (anyNA(occ_folds)) {
      stop("`occ_folds` must not contain missing values.")
    }

    fold_ids <- sort(unique(occ_folds))

    if (length(fold_ids) < 2L) {
      stop("`occ_folds` must contain at least two different folds.")
    }

    for (k in fold_ids) {

      train_idx <- occ_folds != k

      diagnostics[[length(diagnostics) + 1L]] <- check_variance(
        pres_data[train_idx, , drop = FALSE],
        check = "cv_training_presences",
        subset = paste0("held_out_fold_", k)
      )
    }
  }

  # ---- Environmental-domain sample ----

  if (!is.null(env)) {

    if (!inherits(env, "SpatRaster")) {
      stop("`env` must be a terra `SpatRaster` or NULL.")
    }

    missing_env <- setdiff(vars, names(env))

    if (length(missing_env) > 0L) {
      stop(
        "The following training predictors are missing from `env`: ",
        paste(missing_env, collapse = ", ")
      )
    }

    env_use <- env[[vars]]

    sample_env <- function() {
      terra::spatSample(
        env_use,
        size = env_sample_n,
        method = "random",
        na.rm = TRUE,
        values = TRUE,
        xy = FALSE,
        as.df = TRUE
      )
    }

    if (is.null(seed)) {

      env_values <- sample_env()

    } else {

      had_seed <- exists(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )

      if (had_seed) {
        old_seed <- get(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      }

      on.exit(
        {
          if (had_seed) {
            assign(
              ".Random.seed",
              old_seed,
              envir = .GlobalEnv
            )
          } else if (exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
          )) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        },
        add = TRUE
      )

      set.seed(seed)
      env_values <- sample_env()
    }

    env_values <- env_values[, vars, drop = FALSE]

    diagnostics[[length(diagnostics) + 1L]] <- check_variance(
      env_values,
      check = "environment",
      subset = "random_sample"
    )
  }

  # ---- Combine results ----

  diagnostics <- do.call(rbind, diagnostics)

  bad_rows <- diagnostics[
    diagnostics$failed,
    ,
    drop = FALSE
  ]

  remove <- unique(bad_rows$variable)
  keep <- setdiff(vars, remove)

  if (nrow(bad_rows) > 0L) {

    check_names <- paste(
      bad_rows$check,
      bad_rows$subset,
      sep = ": "
    )

    removed_by_check <- split(
      bad_rows$variable,
      check_names
    )

    removed_by_check <- lapply(
      removed_by_check,
      unique
    )

  } else {

    removed_by_check <- list()
  }

  if (isTRUE(verbose)) {

    message(
      "Predictors checked: ", length(vars),
      "; retained: ", length(keep),
      "; removed: ", length(remove), "."
    )

    if (length(remove) > 0L) {
      message(
        "Removed predictors: ",
        paste(remove, collapse = ", ")
      )
    }
  }

  list(
    keep = keep,
    remove = remove,
    removed_by_check = removed_by_check,
    diagnostics = diagnostics
  )
}
