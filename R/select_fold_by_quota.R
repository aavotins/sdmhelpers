#' Select Complete Fold Groups to Approximate a Target Sample Proportion
#'
#' Selects complete fold groups so that the number of rows belonging to the
#' selected folds is as close as possible to a requested proportion of all
#' rows. Fold groups are never split between the selected and unselected sets.
#'
#' @description
#' This function is intended for creating an independent testing subset when
#' observations have already been assigned to spatial, temporal, environmental,
#' or other grouped folds. It searches for a subset of complete fold groups
#' whose combined number of observations approximates `target_prop`.
#'
#' Candidate subsets are generated using repeated random orderings of the fold
#' groups. For each ordering, groups are greedily added while their cumulative
#' size does not exceed the upper bound. If the resulting subset is smaller
#' than the lower bound, one additional group is selected to improve the
#' solution.
#'
#' Because folds are indivisible, an exact match to the target proportion may
#' not be possible. The function therefore returns the best solution found
#' during at most `max_tries` randomized attempts. In some cases, particularly
#' when one or more folds are very large, no combination may satisfy
#' `bounds`. The returned summary indicates whether the selected solution lies
#' within the requested bounds.
#'
#' Missing fold identifiers are excluded by default. When `include_na = TRUE`,
#' all missing identifiers are treated as one additional fold group.
#'
#' @param folds A vector of fold identifiers, with one value per observation.
#'   Fold identifiers may be numeric, character, factor, logical, or another
#'   atomic vector type that can be processed by [base::table()].
#'
#' @param target_prop A single numeric value strictly between `0` and `1`
#'   giving the desired proportion of rows assigned to the selected folds.
#'   It must fall within the interval specified by `bounds`. The default is
#'   `0.25`.
#'
#' @param bounds A numeric vector of length two giving the minimum and maximum
#'   acceptable proportions of rows assigned to the selected folds. Values
#'   must satisfy
#'   `0 <= bounds[1] <= target_prop <= bounds[2] <= 1`.
#'   The default is `c(0.20, 0.30)`.
#'
#' @param max_tries A positive integer giving the maximum number of randomized
#'   fold orderings to evaluate. The search may stop earlier if it finds a
#'   solution containing exactly the target number of rows and lying within
#'   `bounds`. The default is `2000`.
#'
#' @param include_na Logical. If `FALSE`, observations with missing fold
#'   identifiers are excluded from the fold-frequency table and are never
#'   selected. If `TRUE`, all missing fold identifiers are treated as a single
#'   fold group represented by `"<NA>"` in `selected_fold_ids`.
#'
#' @param seed An optional single integer used to initialize the random-number
#'   generator before the search. Use a fixed value to make the selected folds
#'   reproducible. If `NULL`, the current random-number generator state is
#'   used.
#'
#' @param print_report Logical. If `TRUE`, print a summary of the selected
#'   solution to the console. The returned object is invisible regardless of
#'   this setting.
#'
#' @return
#' Invisibly returns an object of class `"fold_quota_selection"`, which is a
#' list with the following components:
#'
#' - `selected_fold_ids`: A character vector containing the identifiers of the selected fold
#'     groups. If missing values are included and selected, they are represented
#'     by `"<NA>"`.
#'
#' - `selected_idx`: A logical vector with the same length as `folds`. Values are `TRUE` for
#'     rows belonging to selected fold groups and `FALSE` otherwise.
#'
#' - `summary`: A list containing:
#'
#'  - `n_rows_total`: Number of rows included in the selection problem. Missing
#'  fold identifiers are omitted when `include_na = FALSE`.
#'
#'  - `n_rows_selected`: Number of rows belonging to the selected folds.
#'
#'  - `prop_selected`: Proportion of included rows belonging to the selected folds.
#'
#'  - `bounds`: Requested lower and upper proportional bounds.
#'
#'  - `row_bounds`: Integer lower and upper bounds used during selection.
#'
#'  - `target_prop`: Requested target proportion.
#'
#'  - `target_n`: Target number of selected rows after rounding.
#'
#'  - `within_bounds`: Logical value indicating whether the returned solution
#'  lies within the requested integer bounds.
#'
#'  - `exact_target`: Logical value indicating whether the selected row count equals `target_n`.
#'
#'  - `n_groups`: Total number of fold groups considered.
#'
#'  - `n_groups_selected`: Number of fold groups selected.
#'
#'  - `iterations_attempted`: Number of randomized attempts performed.
#'
#'  - `best_iteration`: Attempt on which the best returned solution was first found.
#'
#' @details
#' The row-count constraints are calculated as:
#'
#' \preformatted{
#' lower   = ceiling(bounds[1] * N)
#' upper   = floor(bounds[2] * N)
#' target  = round(target_prop * N)
#' }
#'
#' where `N` is the number of observations included in the selection problem.
#'
#' The procedure is a randomized heuristic rather than an exhaustive
#' combinatorial search. Consequently, it does not guarantee the globally
#' optimal combination of folds. Increasing `max_tries` can improve the
#' probability of finding a good combination when many fold groups are
#' available.
#'
#' The same fold identifiers and seed produce the same result under the same
#' R version and random-number generator settings. Calling this function with
#' a non-`NULL` `seed` changes the global random-number generator state.
#'
#' Fold identifiers in `selected_fold_ids` are returned as character values
#' because [base::table()] uses group labels as names. Use `selected_idx` when
#' subsetting the original data, because it preserves the direct relationship
#' with the input rows and avoids type-conversion issues.
#'
#' @section Missing fold identifiers:
#' When `include_na = FALSE`, rows with missing fold identifiers:
#'
#' - are not included in `n_rows_total`;
#' - do not affect the target or bounds;
#' - receive `FALSE` in `selected_idx`.
#'
#' When `include_na = TRUE`, all missing fold identifiers form one complete
#' group. Therefore, either all or none of the missing-fold rows are selected.
#'
#' @section Interpretation:
#' This function selects complete fold groups based only on their sizes. It does
#' not assess whether the selected observations are environmentally,
#' geographically, temporally, or statistically independent. Independence must
#' arise from the procedure used to define `folds`.
#'
#' @seealso
#' [base::table()], [base::sample.int()], [base::set.seed()]
#'
#' @examples
#' ## Example with equally sized folds
#' folds <- rep(letters[1:10], each = 10)
#'
#' selection <- select_folds_by_quota(
#'   folds,
#'   target_prop = 0.25,
#'   bounds = c(0.20, 0.30),
#'   seed = 123,
#'   print_report = FALSE
#' )
#'
#' selection$selected_fold_ids
#' selection$summary
#' table(selection$selected_idx)
#'
#' ## Use selected_idx to divide a data frame
#' dat <- data.frame(
#'   occurrence_id = seq_along(folds),
#'   fold = folds
#' )
#'
#' independent_test <- dat[selection$selected_idx, , drop = FALSE]
#' model_training <- dat[!selection$selected_idx, , drop = FALSE]
#'
#' ## Verify that complete folds were selected
#' unique(independent_test$fold)
#' intersect(
#'   unique(independent_test$fold),
#'   unique(model_training$fold)
#' )
#'
#' ## Unequal fold sizes
#' unequal_folds <- rep(
#'   paste0("fold_", 1:8),
#'   times = c(5, 8, 12, 15, 18, 20, 7, 15)
#' )
#'
#' unequal_selection <- select_folds_by_quota(
#'   unequal_folds,
#'   target_prop = 0.25,
#'   bounds = c(0.20, 0.30),
#'   max_tries = 5000,
#'   seed = 42,
#'   print_report = FALSE
#' )
#'
#' unequal_selection$summary$prop_selected
#' unequal_selection$summary$within_bounds
#'
#' ## Missing fold identifiers excluded
#' folds_with_na <- c("A", "A", "B", "B", "C", NA, NA)
#'
#' without_na <- select_folds_by_quota(
#'   folds_with_na,
#'   target_prop = 0.4,
#'   bounds = c(0.3, 0.5),
#'   include_na = FALSE,
#'   seed = 1,
#'   print_report = FALSE
#' )
#'
#' without_na$selected_idx
#'
#' ## Treat all missing identifiers as one fold
#' with_na <- select_folds_by_quota(
#'   folds_with_na,
#'   target_prop = 0.4,
#'   bounds = c(0.3, 0.5),
#'   include_na = TRUE,
#'   seed = 1,
#'   print_report = FALSE
#' )
#'
#' with_na$selected_fold_ids
#' with_na$selected_idx
#'
#' @export
select_folds_by_quota <- function(folds,
                                  target_prop = 0.25,
                                  bounds = c(0.20, 0.30),
                                  max_tries = 2000L,
                                  include_na = FALSE,
                                  seed = NULL,
                                  print_report = TRUE) {

  # Validate folds ----------------------------------------------------------
  if (is.null(folds)) {
    stop("`folds` must not be NULL.", call. = FALSE)
  }

  if (!is.atomic(folds) || is.list(folds)) {
    stop("`folds` must be an atomic vector.", call. = FALSE)
  }

  if (length(folds) == 0L) {
    stop("`folds` must contain at least one value.", call. = FALSE)
  }

  # Validate target proportion ---------------------------------------------
  if (!is.numeric(target_prop) ||
      length(target_prop) != 1L ||
      is.na(target_prop) ||
      !is.finite(target_prop) ||
      target_prop <= 0 ||
      target_prop >= 1) {
    stop(
      "`target_prop` must be one finite numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  # Validate bounds --------------------------------------------------------
  if (!is.numeric(bounds) ||
      length(bounds) != 2L ||
      anyNA(bounds) ||
      any(!is.finite(bounds))) {
    stop(
      "`bounds` must be a finite numeric vector of length two.",
      call. = FALSE
    )
  }

  if (bounds[1L] < 0 ||
      bounds[2L] > 1 ||
      bounds[1L] > bounds[2L]) {
    stop(
      "`bounds` must satisfy 0 <= bounds[1] <= bounds[2] <= 1.",
      call. = FALSE
    )
  }

  if (target_prop < bounds[1L] || target_prop > bounds[2L]) {
    stop(
      "`target_prop` must lie within the interval specified by `bounds`.",
      call. = FALSE
    )
  }

  # Validate number of attempts --------------------------------------------
  if (!is.numeric(max_tries) ||
      length(max_tries) != 1L ||
      is.na(max_tries) ||
      !is.finite(max_tries) ||
      max_tries < 1 ||
      max_tries != floor(max_tries)) {
    stop(
      "`max_tries` must be one positive integer.",
      call. = FALSE
    )
  }

  max_tries <- as.integer(max_tries)

  # Validate logical arguments ---------------------------------------------
  if (!is.logical(include_na) ||
      length(include_na) != 1L ||
      is.na(include_na)) {
    stop(
      "`include_na` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(print_report) ||
      length(print_report) != 1L ||
      is.na(print_report)) {
    stop(
      "`print_report` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  # Validate and apply seed ------------------------------------------------
  if (!is.null(seed)) {
    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed != floor(seed)) {
      stop(
        "`seed` must be NULL or one finite integer.",
        call. = FALSE
      )
    }

    set.seed(as.integer(seed))
  }

  # Construct fold-frequency table -----------------------------------------
  if (include_na) {
    tab <- table(folds, useNA = "ifany")

    fold_names <- names(tab)
    fold_names[is.na(fold_names)] <- "<NA>"
    names(tab) <- fold_names
  } else {
    tab <- table(folds, useNA = "no")
  }

  if (length(tab) == 0L) {
    stop(
      paste0(
        "No non-missing fold identifiers were found. ",
        "Set `include_na = TRUE` to treat missing identifiers as one group."
      ),
      call. = FALSE
    )
  }

  # Calculate target row counts --------------------------------------------
  n_total <- sum(tab)

  lower <- ceiling(bounds[1L] * n_total)
  upper <- floor(bounds[2L] * n_total)
  target_n <- round(target_prop * n_total)

  # At least one row must be selected
  lower <- max(1L, lower)

  if (upper < 1L) {
    stop(
      "The upper bound is less than one row; increase the upper `bounds` value.",
      call. = FALSE
    )
  }

  if (lower > upper) {
    stop(
      paste0(
        "The proportional bounds produce an empty integer row-count interval: ",
        "lower = ", lower, " and upper = ", upper, ". ",
        "Widen `bounds` or use more observations."
      ),
      call. = FALSE
    )
  }

  sizes <- as.integer(tab)
  ids <- names(tab)
  n_groups <- length(sizes)

  best <- list(
    score = Inf,
    sel = rep(FALSE, n_groups),
    sum_n = 0L,
    iteration = NA_integer_
  )

  # Evaluate one randomly ordered greedy solution --------------------------
  try_once <- function(order_idx) {
    selected <- rep(FALSE, n_groups)
    selected_n <- 0L

    # Greedily add groups without exceeding the upper bound
    for (j in order_idx) {
      if (selected_n + sizes[j] <= upper) {
        selected[j] <- TRUE
        selected_n <- selected_n + sizes[j]
      }
    }

    # Try adding one group if the solution remains below the lower bound
    if (selected_n < lower) {
      remaining <- which(!selected)

      if (length(remaining) > 0L) {
        candidate_totals <- selected_n + sizes[remaining]

        fits <- remaining[
          candidate_totals >= lower &
            candidate_totals <= upper
        ]

        if (length(fits) > 0L) {
          totals_that_fit <- selected_n + sizes[fits]

          j <- fits[
            which.min(abs(totals_that_fit - target_n))
          ]
        } else {
          excess_penalty <- pmax(candidate_totals - upper, 0L)

          candidate_scores <-
            abs(candidate_totals - target_n) + excess_penalty

          j <- remaining[which.min(candidate_scores)]
        }

        selected[j] <- TRUE
        selected_n <- selected_n + sizes[j]
      }
    }

    score <- abs(selected_n - target_n)

    if (selected_n < lower) {
      score <- score + (lower - selected_n)
    }

    if (selected_n > upper) {
      score <- score + (selected_n - upper)
    }

    list(
      score = score,
      sel = selected,
      sum_n = selected_n
    )
  }

  # Repeated randomized search ---------------------------------------------
  iterations_attempted <- 0L

  for (attempt in seq_len(max_tries)) {
    iterations_attempted <- attempt

    candidate <- try_once(sample.int(n_groups))

    if (candidate$score < best$score) {
      best <- list(
        score = candidate$score,
        sel = candidate$sel,
        sum_n = candidate$sum_n,
        iteration = attempt
      )
    }

    exact_valid_solution <-
      best$sum_n == target_n &&
      best$sum_n >= lower &&
      best$sum_n <= upper

    if (exact_valid_solution) {
      break
    }
  }

  selected_ids <- ids[best$sel]

  # Map selected groups back to input rows ---------------------------------
  if (include_na) {
    row_ids <- as.character(folds)
    row_ids[is.na(folds)] <- "<NA>"

    selected_idx <- row_ids %in% selected_ids
  } else {
    selected_idx <- folds %in% selected_ids
  }

  within_bounds <-
    best$sum_n >= lower &&
    best$sum_n <= upper

  exact_target <- best$sum_n == target_n

  out <- list(
    selected_fold_ids = selected_ids,
    selected_idx = selected_idx,
    summary = list(
      n_rows_total = n_total,
      n_rows_selected = best$sum_n,
      prop_selected = best$sum_n / n_total,
      bounds = bounds,
      row_bounds = c(lower = lower, upper = upper),
      target_prop = target_prop,
      target_n = target_n,
      within_bounds = within_bounds,
      exact_target = exact_target,
      n_groups = n_groups,
      n_groups_selected = sum(best$sel),
      iterations_attempted = iterations_attempted,
      best_iteration = best$iteration
    )
  )

  class(out) <- c("fold_quota_selection", "list")

  # Optional console report ------------------------------------------------
  if (print_report) {
    cat("Fold selection summary\n")
    cat("----------------------\n")

    cat(sprintf(
      "Rows selected: %d / %d (%.2f%%)\n",
      out$summary$n_rows_selected,
      out$summary$n_rows_total,
      100 * out$summary$prop_selected
    ))

    cat(sprintf(
      "Allowed rows: %d-%d; target: %d\n",
      out$summary$row_bounds["lower"],
      out$summary$row_bounds["upper"],
      out$summary$target_n
    ))

    cat(sprintf(
      "Bounds: [%.0f%%, %.0f%%]; target: %.0f%%\n",
      100 * out$summary$bounds[1L],
      100 * out$summary$bounds[2L],
      100 * out$summary$target_prop
    ))

    cat(sprintf(
      "Groups selected: %d / %d\n",
      out$summary$n_groups_selected,
      out$summary$n_groups
    ))

    cat(sprintf(
      "Search attempts: %d; best solution found at attempt: %d\n",
      out$summary$iterations_attempted,
      out$summary$best_iteration
    ))

    cat(sprintf(
      "Within requested bounds: %s\n",
      if (out$summary$within_bounds) "yes" else "no"
    ))

    cat(sprintf(
      "Exact target row count: %s\n",
      if (out$summary$exact_target) "yes" else "no"
    ))

    cat("Selected fold IDs:\n")

    if (length(out$selected_fold_ids) > 0L) {
      cat("  ", paste(out$selected_fold_ids, collapse = ", "), "\n")
    } else {
      cat("  <none>\n")
    }
  }

  invisible(out)
}
