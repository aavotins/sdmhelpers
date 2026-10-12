#' Screen environmental predictors for insufficient variance
#'
#' Calculates predictor standard deviations separately among presences,
#' background records, and their combination, including spatial CV training
#' partitions and (optionally) an independent test set.
#'
#' @param train An `SDMtune` `SWD` object with presence and background data.
#' @param test An optional independent-test `SWD` object.
#' @param occ_folds An optional fold identifiers for TRAINING presences, in the
#'   order of `train@pa == 1`. A held-out fold is excluded in each CV check.
#' @param bg_folds An optional fold identifiers for TRAINING background records,
#'   in the order of `train@pa == 0`. If omitted, all training backgrounds are
#'   reused in every CV partition (and the diagnostics state this explicitly).
#' @param env An optional `terra::SpatRaster` representing the environmental domain.
#' @param env_sample_n Whole integer >= 2 for random raster-cell sampling.
#' @param seed NULL or one non-negative integer for reproducible `env` sampling;
#'   existing R random-number-generator state is restored.
#' @param sd_tol Finite non-negative standard-deviation cutoff. Values <=
#'   `sd_tol`, or fewer than two finite observations, fail a check.
#' @param check_presences Check presence-only subsets, default TRUE.
#' @param check_background Check background-only subsets, default TRUE.
#' @param check_combined Check joined presence/background subsets, default TRUE.
#' @param presence_affects_selection If FALSE, presence-only failures are
#'   diagnostic-only. Defaults to TRUE for compatibility with previous screening.
#' @param include_test_bg Whether to additionally check training and independent
#'   test backgrounds pooled together. Defaults to FALSE. The pooled check is
#'   diagnostic-only unless `test_affects_selection = TRUE`.
#' @param check_test Whether to diagnose the independent test set separately
#'   (presences/background/combined according to the `check_*` switches).
#'   Defaults to FALSE. Its checks are diagnostic-only by default.
#' @param test_affects_selection Explicitly allow checks involving independent
#'   test records to influence `keep`/`remove`. Defaults to FALSE to avoid
#'   test-informed preprocessing. Setting TRUE forfeits full test independence.
#' @param verbose Print a summary (TRUE/FALSE).
#'
#' @details
#' Checks performed on training observations use only `train@data`.
#' Supplying `env` adds an independent random environmental-domain check; it
#' does not subsample or modify the supplied background observations.
#'
#' CV checks exclude the held-out fold's presences and, when `bg_folds` is
#' supplied, its background records. Without `bg_folds`, the full training
#' background is shared across CV partitions. CV checks never inspect `test`.
#'
#' When `check_presences = FALSE`, presence-only checks are omitted entirely.
#' When `presence_affects_selection = FALSE`, presence-only checks are
#' performed but do not influence predictor selection.
#'
#' Each row in the returned `diagnostics` data frame contains a logical
#' column named `used_for_selection`. This column indicates whether the
#' corresponding variance check contributes to predictor selection.
#' A predictor is included in `remove` if it fails at least one check
#' for which `used_for_selection` is `TRUE`. Checks marked `FALSE` are
#' retained for diagnostic purposes but do not influence `keep` or `remove`.
#'
#' `test_affects_selection = TRUE` is an explicit opt-in and is generally
#' inappropriate if `test` is meant as a final independent evaluation set.
#' Select predictors *within training data for each CV fold* if estimating
#' the predictive performance of an entire selection/model-fitting pipeline.
#'
#' @return A named list containing:
#' \describe{
#'   \item{keep}{Predictors retained after variance screening.}
#'   \item{remove}{Predictors excluded by selection-relevant variance checks.}
#'   \item{removed_by_check}{Predictors removed by each selection-relevant check.}
#'   \item{diagnostics}{A data frame containing results of all performed
#'     variance checks, including `failed` and `used_for_selection` logical
#'     columns. The latter identifies whether each check contributes to
#'     predictor selection.}
#' }
#'
#' @importClassesFrom SDMtune SWD
#' @seealso [SDMtune::prepareSWD()], [terra::spatSample()]
#' @examples
#' # Self-contained example with package dependencies, no external files.
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10,
#'                  ymin = 0, ymax = 10)
#' terra::values(r) <- seq_len(terra::ncell(r))
#' names(r) <- "gradient"
#' r2 <- r
#' terra::values(r2) <- rep(seq_len(10), each = 10)
#' names(r2) <- "stripes"
#' env <- c(r, r2)
#' p <- data.frame(x = c(1.5, 3.5, 5.5, 7.5),
#'                 y = c(1.5, 3.5, 5.5, 7.5))
#' a <- data.frame(x = c(0.5, 2.5, 4.5, 6.5, 8.5),
#'                 y = c(8.5, 6.5, 4.5, 2.5, 0.5))
#' train <- SDMtune::prepareSWD("Example species", env, p = p, a = a,
#'                              verbose = FALSE)
#' result <- screen_egv_variance(train, occ_folds = c(1, 1, 2, 2),
#'                               env = env, env_sample_n = 50,
#'                               seed = 1, verbose = FALSE)
#' result$keep
#' head(result$diagnostics)
#'
#' @export
screen_egv_variance <- function(
    train, test = NULL, occ_folds = NULL, bg_folds = NULL, env = NULL,
    env_sample_n = 10000L, seed = NULL, sd_tol = 0,
    check_presences = TRUE, check_background = TRUE, check_combined = TRUE,
    presence_affects_selection = TRUE,
    include_test_bg = FALSE, check_test = FALSE,
    test_affects_selection = FALSE, verbose = TRUE
) {
  fail <- function(...) stop(..., call. = FALSE)
  flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      fail("`", name, "` must be TRUE or FALSE.")
    }
  }
  whole_number <- function(x, name, min_value) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < min_value || x != floor(x) || x > .Machine$integer.max) {
      fail("`", name, "` must be a whole integer between ",
           min_value, " and .Machine$integer.max.")
    }
    as.integer(x)
  }
  for (nm in c("check_presences", "check_background", "check_combined",
               "presence_affects_selection", "include_test_bg", "check_test",
               "test_affects_selection", "verbose")) {
    flag(get(nm, inherits = FALSE), nm)
  }
  if (!any(c(check_presences, check_background, check_combined))) {
    fail("At least one of `check_presences`, `check_background`, or ",
         "`check_combined` must be TRUE.")
  }
  env_sample_n <- whole_number(env_sample_n, "env_sample_n", 2L)
  if (!is.null(seed)) seed <- whole_number(seed, "seed", 0L)
  if (!is.numeric(sd_tol) || length(sd_tol) != 1L ||
      !is.finite(sd_tol) || sd_tol < 0) {
    fail("`sd_tol` must be one finite non-negative number.")
  }
  if (!methods::is(train, "SWD")) fail("`train` must be an SDMtune `SWD` object.")
  if (!is.null(test) && !methods::is(test, "SWD")) {
    fail("`test` must be an SDMtune `SWD` object or NULL.")
  }
  if (check_test && is.null(test)) fail("`check_test = TRUE` requires `test`.")
  if (include_test_bg && is.null(test)) {
    fail("`include_test_bg = TRUE` requires `test`.")
  }
  if (!is.null(bg_folds) && is.null(occ_folds)) {
    fail("`bg_folds` requires `occ_folds`.")
  }

<<<<<<< HEAD
  extract <- function(obj, name, vars = NULL) {
    x <- as.data.frame(methods::slot(obj, "data"))
    pa <- methods::slot(obj, "pa")
    if (nrow(x) != length(pa)) {
      fail("Dimensions of `", name, "@data` and `", name, "@pa` disagree.")
=======
  if (length(sd_tol) != 1L ||
      !is.numeric(sd_tol) ||
      !is.finite(sd_tol) ||
      sd_tol < 0) {
    stop("`sd_tol` must be one finite, non-negative number.")
  }

  if (
    !is.numeric(env_sample_n) ||
    length(env_sample_n) != 1L ||
    is.na(env_sample_n) ||
    !is.finite(env_sample_n) ||
    env_sample_n < 2 ||
    env_sample_n != floor(env_sample_n)
  ) {
    stop("`env_sample_n` must be an integer >= 2.")
  }

  env_sample_n <- as.integer(env_sample_n)


  if (!is.null(seed)) {
    if (
      !is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed) ||
      !is.finite(seed) ||
      seed != floor(seed)
    ) {
      stop(
        "`seed` must be NULL or one finite integer.",
        call. = FALSE
      )
    }
  }

  if (
    !is.logical(include_test_bg) ||
    length(include_test_bg) != 1L ||
    is.na(include_test_bg)
  ) {
    ...
  }

  if (
    !is.logical(check_test) ||
    length(check_test) != 1L ||
    is.na(check_test)
  ) {
    ...
  }

  if (
    !is.logical(verbose) ||
    length(verbose) != 1L ||
    is.na(verbose)
  ) {
    ...
  }



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
>>>>>>> 60cdeaf24ae5cbdf34138878f015d1200dd64263
    }
    if (anyNA(pa) || !all(pa %in% c(0, 1))) {
      fail("`", name, "@pa` must contain only 0 and 1 without NA values.")
    }
    if (is.null(vars)) {
      vars <- names(x)
      if (!length(vars) || anyDuplicated(vars)) {
        fail("`train@data` must have non-empty, unique predictor names.")
      }
    } else {
      missing <- setdiff(vars, names(x))
      if (length(missing)) {
        fail("Predictors missing from `", name, "@data`: ",
             paste(missing, collapse = ", "))
      }
      x <- x[, vars, drop = FALSE]
    }
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      fail("Non-numeric predictors in `", name, "@data`: ",
           paste(names(x)[!numeric_cols], collapse = ", "))
    }
    list(data = x, pa = pa, vars = vars)
  }
  tr <- extract(train, "train")
  vars <- tr$vars
  te <- if (!is.null(test)) extract(test, "test", vars) else NULL
  train_p <- tr$data[tr$pa == 1, , drop = FALSE]
  train_b <- tr$data[tr$pa == 0, , drop = FALSE]
  if (!nrow(train_p) || !nrow(train_b)) {
    fail("`train` must contain both presence and background records.")
  }
  if (!is.null(occ_folds)) {
    if (length(occ_folds) != nrow(train_p) || anyNA(occ_folds) ||
        !(is.numeric(occ_folds) || is.character(occ_folds) || is.factor(occ_folds))) {
      fail("`occ_folds` must provide one non-missing numeric or character fold ",
           "identifier per training presence.")
    }
    folds <- unique(as.character(occ_folds))
    if (length(folds) < 2L) fail("`occ_folds` must contain at least two folds.")
    occ_folds <- as.character(occ_folds)
    if (!is.null(bg_folds)) {
      if (length(bg_folds) != nrow(train_b) || anyNA(bg_folds) ||
          !(is.numeric(bg_folds) || is.character(bg_folds) || is.factor(bg_folds))) {
        fail("`bg_folds` must provide one non-missing numeric or character fold ",
             "identifier per training background record.")
      }
      bg_folds <- as.character(bg_folds)
      if (!all(bg_folds %in% folds)) {
        fail("All `bg_folds` identifiers must occur in `occ_folds`.")
      }
    }
  }

  diagnostics <- list()
  check_one <- function(x, check, subset, data_source, group, fold = NA_character_,
                        select = TRUE) {
    n_finite <- vapply(x, function(z) sum(is.finite(z)), integer(1))
    sds <- vapply(x, function(z) {
      z <- z[is.finite(z)]
      if (length(z) < 2L) return(NA_real_)
      stats::sd(z)
    }, numeric(1))
    failed <- n_finite < 2L | !is.finite(sds) | sds <= sd_tol
    reason <- rep("passed", length(vars))
    reason[n_finite < 2L] <- "fewer_than_2_finite_values"
    reason[n_finite >= 2L & !is.finite(sds)] <- "undefined_sd"
    reason[n_finite >= 2L & is.finite(sds) & sds <= sd_tol] <-
      "sd_below_or_equal_to_tolerance"
    diagnostics[[length(diagnostics) + 1L]] <<- data.frame(
      check = check, subset = subset, data_source = data_source,
      group = group, fold = fold, used_for_selection = select,
      variable = vars, n_finite = n_finite, sd = sds,
      failed = failed, reason = reason, row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  check_groups <- function(p, b, check, subset, data_source,
                           fold = NA_character_, select = TRUE) {
    if (check_presences) check_one(
      p, check, subset, data_source, "presences", fold,
      select && presence_affects_selection
    )
    if (check_background) check_one(
      b, check, subset, data_source, "background", fold, select
    )
    if (check_combined) check_one(
      rbind(p, b), check, subset, data_source, "combined", fold, select
    )
  }

  check_groups(train_p, train_b, "training", "all", "training")
  if (!is.null(occ_folds)) {
    for (k in folds) {
      p <- train_p[occ_folds != k, , drop = FALSE]
      b <- if (is.null(bg_folds)) train_b else
        train_b[bg_folds != k, , drop = FALSE]
      check_groups(p, b, "cv_training", paste0(
        "held_out_fold_", k,
        if (is.null(bg_folds)) "_shared_background" else ""
      ), "training", fold = k)
    }
  }
  if (!is.null(te)) {
    test_p <- te$data[te$pa == 1, , drop = FALSE]
    test_b <- te$data[te$pa == 0, , drop = FALSE]
    if (check_test) {
      check_groups(test_p, test_b, "independent_test", "all",
                   "independent_test", select = test_affects_selection)
    }
    if (include_test_bg) {
      check_one(rbind(train_b, test_b), "pooled_background",
                "training_plus_independent_test", "training_plus_independent_test",
                "background", select = test_affects_selection)
    }
  }
  if (!is.null(env)) {
    if (!inherits(env, "SpatRaster")) {
      fail("`env` must be a terra `SpatRaster` or NULL.")
    }
    missing <- setdiff(vars, names(env))
    if (length(missing)) {
      fail("Predictors missing from `env`: ", paste(missing, collapse = ", "))
    }
    sample_env <- function() terra::spatSample(
      env[[vars]], size = env_sample_n, method = "random", na.rm = TRUE,
      values = TRUE, xy = FALSE, as.df = TRUE
    )
    if (is.null(seed)) {
      env_values <- sample_env()
    } else {
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
      set.seed(seed)
      env_values <- sample_env()
    }
    env_values <- as.data.frame(env_values)[, vars, drop = FALSE]
    check_one(env_values, "environment", "random_sample",
              "environment", "environment")
  }
  diagnostics <- do.call(rbind, diagnostics)
  bad <- diagnostics[diagnostics$failed & diagnostics$used_for_selection,
                     , drop = FALSE]
  remove <- vars[vars %in% bad$variable]
  keep <- vars[!vars %in% remove]
  if (nrow(bad)) {
    keys <- paste(bad$check, bad$subset, bad$group, sep = ": ")
    removed_by_check <- lapply(split(bad$variable, keys), unique)
  } else {
    removed_by_check <- list()
  }
  if (verbose) {
    message("Predictors checked: ", length(vars), "; retained: ",
            length(keep), "; removed: ", length(remove), ".")
    if (length(remove)) message("Removed predictors: ", paste(remove, collapse = ", "))
  }
  list(keep = keep, remove = remove, removed_by_check = removed_by_check,
       diagnostics = diagnostics)
}
