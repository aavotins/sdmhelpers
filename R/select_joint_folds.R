#' Select joint presence and background folds for independent testing
#'
#' Selects a common set of complete fold IDs for an independent testing
#' dataset. The selected folds are chosen so that the proportion of selected
#' presence records is close to a requested target, while the proportions of
#' both presence and background records remain within user-defined bounds.
#'
#' The same fold IDs are selected for both datasets. Consequently, all
#' presence and background records belonging to a selected fold are assigned
#' to the independent testing dataset.
#'
#' @param pres_folds A vector containing fold identifiers for presence
#'   records. Fold identifiers may be numeric, character, factor, or another
#'   atomic type coercible to character.
#'
#' @param bg_folds A vector containing fold identifiers for background
#'   records. Fold identifiers may be numeric, character, factor, or another
#'   atomic type coercible to character.
#'
#' @param target_prop A single numeric value strictly between `0` and `1`.
#'   The desired proportion of in-scope presence records assigned to the
#'   independent testing dataset. The default is `0.25`.
#'
#' @param pres_bounds A numeric vector of length two giving the minimum and
#'   maximum permitted proportions of presence records in the selected folds.
#'   Values must lie between `0` and `1`, and the first value must not exceed
#'   the second. The default is `c(0.20, 0.30)`.
#'
#' @param bg_bounds A numeric vector of length two giving the minimum and
#'   maximum permitted proportions of background records in the selected
#'   folds. Values must lie between `0` and `1`, and the first value must not
#'   exceed the second. The default is `c(0.20, 0.30)`.
#'
#' @param max_tries A positive whole number giving the maximum number of
#'   randomized fold orderings evaluated by the heuristic search. Larger
#'   values may improve the selected solution when fold sizes are irregular,
#'   at the cost of additional computation. The default is `3000`.
#'
#' @param bg_match_weight A non-negative numeric value controlling the
#'   importance of matching the selected background proportion to the
#'   selected presence proportion. A value of `0` disables proportional
#'   matching, although `bg_bounds` are still considered. Larger values place
#'   progressively greater emphasis on making the two selected proportions
#'   similar. The default is `1`.
#'
#' @param include_na Logical. If `FALSE`, records with missing fold IDs are
#'   excluded from fold-size calculations and cannot be selected. Their
#'   corresponding values in the returned row-selection vectors are
#'   `FALSE`. If `TRUE`, missing fold IDs are treated as one additional fold
#'   represented internally and in the returned fold IDs by `"<NA>"`.
#'   The default is `FALSE`.
#'
#' @param seed `NULL` or a single integer-like value passed to
#'   [base::set.seed()]. Supplying a seed makes the randomized search
#'   reproducible. When `NULL`, the current random-number generator state is
#'   used.
#'
#' @param print_report Logical. If `TRUE`, print a summary of the selected
#'   folds, sample sizes, proportions, iteration number, and final score.
#'   The default is `TRUE`.
#'
#' @details
#' This function is intended for spatial or otherwise grouped model
#' evaluation in which complete folds must be withheld from both presence and
#' background datasets. It does not select individual records independently.
#'
#' Fold sizes are first counted separately for the presence and background
#' data. The union of the presence and background fold IDs defines the
#' candidate folds. A fold that occurs in only one dataset therefore has size
#' zero in the other dataset.
#'
#' The function performs up to `max_tries` randomized greedy searches. For
#' each randomized ordering, candidate folds are added when doing so improves
#' the objective score or when additional presence records are needed to reach
#' the lower presence bound.
#'
#' The objective score combines:
#'
#' \enumerate{
#'   \item deviation of the selected presence count from the requested
#'     presence target;
#'   \item penalties for selected presence counts outside `pres_bounds`;
#'   \item penalties for selected background counts outside `bg_bounds`; and
#'   \item the absolute difference between selected presence and background
#'     proportions, weighted by `bg_match_weight`.
#' }
#'
#' Because complete folds are indivisible, it may be impossible to satisfy
#' all requested bounds or to attain `target_prop` exactly. The returned
#' result is the best solution encountered by the heuristic search; it is not
#' guaranteed to be the global optimum.
#'
#' Search quality may be improved by increasing `max_tries`, especially when:
#'
#' \itemize{
#'   \item the number of folds is large;
#'   \item fold sizes vary substantially;
#'   \item presence and background fold-size distributions differ; or
#'   \item the requested bounds are narrow.
#' }
#'
#' @section Interpretation of missing fold IDs:
#'
#' With `include_na = FALSE`, missing fold IDs are outside the selection
#' problem. They are included in `n_pres_total` or `n_bg_total`, but excluded
#' from `n_pres_in_scope` or `n_bg_in_scope`.
#'
#' With `include_na = TRUE`, all missing fold IDs are treated as belonging to
#' one common fold named `"<NA>"`. Therefore, selecting that fold selects every
#' record with a missing fold ID in the corresponding dataset.
#'
#' To avoid ambiguity when `include_na = TRUE`, the input fold identifiers
#' should not contain an actual, non-missing fold ID equal to `"<NA>"`.
#'
#' @section Random-number generation:
#'
#' If `seed` is supplied, the function calls [base::set.seed()] and therefore
#' changes R's global random-number generator state. The same inputs, seed,
#' and R version should produce the same result.
#'
#' @return
#' Invisibly returns a named list with the following elements:
#'
#'   - `selected_fold_ids`: A character vector containing the selected fold IDs.
#'
#'   - `pres_selected_idx`: A logical vector of length `length(pres_folds)`. Values are `TRUE` for
#'     presence records belonging to selected folds.
#'
#'   - `bg_selected_idx`: A logical vector of length `length(bg_folds)`. Values are `TRUE` for
#'     background records belonging to selected folds.
#'
#'   - `summary`: A list containing:
#'
#'       - `n_pres_total`: Total number of elements in `pres_folds`, including missing values.
#'
#'       - `n_bg_total`: Total number of elements in `bg_folds`, including missing values.
#'
#'       - `n_pres_in_scope`: Number of presence records considered during selection after
#'         applying the requested missing-value handling.
#'
#'       - `n_bg_in_scope`: Number of background records considered during selection after
#'         applying the requested missing-value handling.
#'
#'       - `n_pres_selected`: Number of selected in-scope presence records.
#'
#'       - `n_bg_selected`: Number of selected in-scope background records.
#'
#'       - `prop_pres_selected`: Proportion of in-scope presence records selected.
#'
#'       - `prop_bg_selected`: Proportion of in-scope background records selected.
#'
#'       - `target_prop`: Requested presence selection proportion.
#'
#'       - `bounds_pres`: Requested presence proportion bounds.
#'
#'       - `bounds_bg`: Requested background proportion bounds.
#'
#'       - `target_pres_abs`: Rounded target number of presence records.
#'
#'       - `target_bg_abs_range`: Integer lower and upper bounds for the selected background count.
#'
#'       - `iterations`: Iteration at which the best solution was first encountered.
#'
#'       - `score`: Score of the selected solution. Smaller values indicate
#'         better agreement with the requested target, bounds, and
#'         proportional matching criterion.
#'
#' @note
#' The function optimizes the number of records assigned to the independent
#' dataset, not the number of selected folds. A small number of large folds
#' can therefore be preferred over a larger number of small folds.
#'
#' The returned logical vectors can be used directly to separate testing and
#' training records:
#'
#' `pres_test  <- pres_data[result$pres_selected_idx, ]`
#' `pres_train <- pres_data[!result$pres_selected_idx, ]`
#' `bg_test    <- bg_data[result$bg_selected_idx, ]`
#' `bg_train   <- bg_data[!result$bg_selected_idx, ]`
#'
#' @seealso
#' [base::table()], [base::set.seed()], [base::sample.int()]
#'
#' @examples
#' ## Create fold assignments with unequal fold sizes.
#' pres_folds <- rep(
#'   paste0("fold_", 1:8),
#'   times = c(12, 15, 9, 18, 14, 11, 16, 10)
#' )
#'
#' bg_folds <- rep(
#'   paste0("fold_", 1:8),
#'   times = c(80, 95, 70, 120, 90, 75, 110, 85)
#' )
#'
#' result <- select_joint_folds(
#'   pres_folds = pres_folds,
#'   bg_folds = bg_folds,
#'   target_prop = 0.25,
#'   pres_bounds = c(0.20, 0.30),
#'   bg_bounds = c(0.20, 0.30),
#'   max_tries = 500,
#'   seed = 123,
#'   print_report = FALSE
#' )
#'
#' result$selected_fold_ids
#' result$summary
#'
#' ## Extract selected records.
#' pres_test_folds <- pres_folds[result$pres_selected_idx]
#' bg_test_folds <- bg_folds[result$bg_selected_idx]
#'
#' unique(pres_test_folds)
#' unique(bg_test_folds)
#'
#' ## Verify that both datasets use the same selected fold IDs.
#' setequal(
#'   unique(pres_test_folds),
#'   unique(bg_test_folds)
#' )
#'
#' ## Presence and background tables do not need to contain exactly the
#' ## same folds. A missing fold in one dataset receives a size of zero.
#' pres_folds2 <- rep(letters[1:6], c(7, 13, 10, 15, 9, 11))
#' bg_folds2 <- rep(letters[c(1:5, 7)], c(60, 90, 75, 100, 70, 50))
#'
#' result2 <- select_joint_folds(
#'   pres_folds = pres_folds2,
#'   bg_folds = bg_folds2,
#'   target_prop = 0.25,
#'   pres_bounds = c(0.15, 0.35),
#'   bg_bounds = c(0.15, 0.35),
#'   max_tries = 500,
#'   seed = 42,
#'   print_report = FALSE
#' )
#'
#' result2$selected_fold_ids
#'
#' ## Missing fold IDs can be excluded.
#' pres_with_na <- c(pres_folds, NA, NA)
#' bg_with_na <- c(bg_folds, NA)
#'
#' result3 <- select_joint_folds(
#'   pres_folds = pres_with_na,
#'   bg_folds = bg_with_na,
#'   include_na = FALSE,
#'   max_tries = 500,
#'   seed = 7,
#'   print_report = FALSE
#' )
#'
#' result3$pres_selected_idx[is.na(pres_with_na)]
#'
#' ## Alternatively, all missing IDs can be treated as one selectable fold.
#' result4 <- select_joint_folds(
#'   pres_folds = pres_with_na,
#'   bg_folds = bg_with_na,
#'   include_na = TRUE,
#'   max_tries = 500,
#'   seed = 7,
#'   print_report = FALSE
#' )
#'
#' result4$selected_fold_ids
#'
#' @export
select_joint_folds <- function(
    pres_folds,
    bg_folds,
    target_prop = 0.25,
    pres_bounds = c(0.20, 0.30),
    bg_bounds = c(0.20, 0.30),
    max_tries = 3000,
    bg_match_weight = 1,
    include_na = FALSE,
    seed = NULL,
    print_report = TRUE
) {
  # Input validation --------------------------------------------------------

  if (!is.atomic(pres_folds) || is.list(pres_folds)) {
    stop("`pres_folds` must be an atomic vector.", call. = FALSE)
  }

  if (!is.atomic(bg_folds) || is.list(bg_folds)) {
    stop("`bg_folds` must be an atomic vector.", call. = FALSE)
  }

  if (
    !is.numeric(target_prop) ||
    length(target_prop) != 1L ||
    is.na(target_prop) ||
    !is.finite(target_prop) ||
    target_prop <= 0 ||
    target_prop >= 1
  ) {
    stop(
      "`target_prop` must be one finite numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  validate_bounds <- function(x, name) {
    if (
      !is.numeric(x) ||
      length(x) != 2L ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x < 0 | x > 1) ||
      x[1L] > x[2L]
    ) {
      stop(
        sprintf(
          "`%s` must contain two finite proportions between 0 and 1, ",
          name
        ),
        "with the lower bound not exceeding the upper bound.",
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  validate_bounds(pres_bounds, "pres_bounds")
  validate_bounds(bg_bounds, "bg_bounds")

  if (
    !is.numeric(max_tries) ||
    length(max_tries) != 1L ||
    is.na(max_tries) ||
    !is.finite(max_tries) ||
    max_tries < 1 ||
    max_tries != floor(max_tries)
  ) {
    stop(
      "`max_tries` must be one positive whole number.",
      call. = FALSE
    )
  }

  max_tries <- as.integer(max_tries)

  if (
    !is.numeric(bg_match_weight) ||
    length(bg_match_weight) != 1L ||
    is.na(bg_match_weight) ||
    !is.finite(bg_match_weight) ||
    bg_match_weight < 0
  ) {
    stop(
      "`bg_match_weight` must be one finite, non-negative numeric value.",
      call. = FALSE
    )
  }

  if (
    !is.logical(include_na) ||
    length(include_na) != 1L ||
    is.na(include_na)
  ) {
    stop(
      "`include_na` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (
    !is.logical(print_report) ||
    length(print_report) != 1L ||
    is.na(print_report)
  ) {
    stop(
      "`print_report` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    if (
      !is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed) ||
      !is.finite(seed)
    ) {
      stop(
        "`seed` must be NULL or one finite numeric value.",
        call. = FALSE
      )
    }

    set.seed(seed)
  }

  # Prepare fold IDs -------------------------------------------------------

  pres_ids <- as.character(pres_folds)
  bg_ids <- as.character(bg_folds)

  if (!include_na) {
    pres_ids <- pres_ids[!is.na(pres_ids)]
    bg_ids <- bg_ids[!is.na(bg_ids)]
  } else {
    if (
      any(pres_ids == "<NA>", na.rm = TRUE) ||
      any(bg_ids == "<NA>", na.rm = TRUE)
    ) {
      stop(
        "An actual fold ID equal to \"<NA>\" cannot be used when ",
        "`include_na = TRUE`.",
        call. = FALSE
      )
    }

    pres_ids[is.na(pres_ids)] <- "<NA>"
    bg_ids[is.na(bg_ids)] <- "<NA>"
  }

  # Fold-size tables -------------------------------------------------------

  tab_p <- table(pres_ids, useNA = "no")
  tab_b <- table(bg_ids, useNA = "no")

  if (length(tab_p) == 0L) {
    stop(
      "No presence folds remain after applying missing-value handling.",
      call. = FALSE
    )
  }

  if (length(tab_b) == 0L) {
    stop(
      "No background folds remain after applying missing-value handling.",
      call. = FALSE
    )
  }

  all_ids <- union(names(tab_p), names(tab_b))

  if (length(all_ids) == 0L) {
    stop("No fold IDs are available for selection.", call. = FALSE)
  }

  p_sizes <- stats::setNames(integer(length(all_ids)), all_ids)
  b_sizes <- stats::setNames(integer(length(all_ids)), all_ids)

  p_sizes[names(tab_p)] <- as.integer(tab_p)
  b_sizes[names(tab_b)] <- as.integer(tab_b)

  Np_total <- sum(p_sizes)
  Nb_total <- sum(b_sizes)

  if (Np_total <= 0L) {
    stop("The in-scope presence sample is empty.", call. = FALSE)
  }

  if (Nb_total <= 0L) {
    stop("The in-scope background sample is empty.", call. = FALSE)
  }

  lower_p <- ceiling(pres_bounds[1L] * Np_total)
  upper_p <- floor(pres_bounds[2L] * Np_total)
  target_p <- round(target_prop * Np_total)

  lower_b <- ceiling(bg_bounds[1L] * Nb_total)
  upper_b <- floor(bg_bounds[2L] * Nb_total)

  if (lower_p > upper_p) {
    stop(
      "`pres_bounds` produce an empty integer count interval for the ",
      "available number of presence records.",
      call. = FALSE
    )
  }

  if (lower_b > upper_b) {
    stop(
      "`bg_bounds` produce an empty integer count interval for the ",
      "available number of background records.",
      call. = FALSE
    )
  }

  n_groups <- length(all_ids)

  # Objective function -----------------------------------------------------

  score_of <- function(sumP, sumB) {
    pres_gap <- abs(sumP - target_p)

    if (sumP < lower_p) {
      pres_gap <- pres_gap + (lower_p - sumP)
    }

    if (sumP > upper_p) {
      pres_gap <- pres_gap + (sumP - upper_p)
    }

    bg_gap <- 0

    if (sumB < lower_b) {
      bg_gap <- bg_gap + (lower_b - sumB)
    }

    if (sumB > upper_b) {
      bg_gap <- bg_gap + (sumB - upper_b)
    }

    pres_prop <- sumP / Np_total
    bg_prop <- sumB / Nb_total
    match_gap <- abs(bg_prop - pres_prop)

    pres_gap +
      bg_gap +
      bg_match_weight * match_gap * max(Np_total, Nb_total)
  }

  # One randomized greedy search ------------------------------------------

  try_once <- function(order_idx) {
    sel <- rep(FALSE, n_groups)
    sumP <- 0L
    sumB <- 0L

    for (j in order_idx) {
      p_add <- p_sizes[j]
      b_add <- b_sizes[j]

      # Presence-free folds may still improve background matching.
      if (sumP + p_add <= upper_p || p_add == 0L) {
        current_score <- score_of(sumP, sumB)
        proposed_score <- score_of(
          sumP + p_add,
          sumB + b_add
        )

        if (proposed_score <= current_score || sumP < lower_p) {
          sel[j] <- TRUE
          sumP <- sumP + p_add
          sumB <- sumB + b_add
        }
      }
    }

    # If the presence lower bound has not been reached, add the remaining
    # fold that gives the smallest resulting score.
    if (sumP < lower_p) {
      remaining <- which(!sel)

      if (length(remaining) > 0L) {
        candidate_scores <- vapply(
          remaining,
          FUN = function(j) {
            score_of(
              sumP + p_sizes[j],
              sumB + b_sizes[j]
            )
          },
          FUN.VALUE = numeric(1)
        )

        j_best <- remaining[which.min(candidate_scores)]

        sel[j_best] <- TRUE
        sumP <- sumP + p_sizes[j_best]
        sumB <- sumB + b_sizes[j_best]
      }
    }

    list(
      sel = sel,
      sumP = sumP,
      sumB = sumB,
      score = score_of(sumP, sumB)
    )
  }

  # Repeated randomized searches ------------------------------------------

  best <- list(
    score = Inf,
    sel = rep(FALSE, n_groups),
    sumP = 0L,
    sumB = 0L,
    iter = NA_integer_
  )

  tolerance <- sqrt(.Machine$double.eps)

  for (iteration in seq_len(max_tries)) {
    ordering <- sample.int(n_groups)
    candidate <- try_once(ordering)

    if (candidate$score < best$score) {
      best <- c(
        candidate,
        list(iter = iteration)
      )
    }

    pres_ok <- (
      best$sumP == target_p &&
        best$sumP >= lower_p &&
        best$sumP <= upper_p
    )

    bg_ok <- (
      best$sumB >= lower_b &&
        best$sumB <= upper_b
    )

    prop_ok <- abs(
      best$sumB / Nb_total -
        best$sumP / Np_total
    ) < tolerance

    if (pres_ok && bg_ok && prop_ok) {
      break
    }
  }

  selected_ids <- all_ids[best$sel]

  # Map selected folds to original rows -----------------------------------

  pres_ids_all <- as.character(pres_folds)
  bg_ids_all <- as.character(bg_folds)

  if (!include_na) {
    pres_idx <- (
      !is.na(pres_ids_all) &
        pres_ids_all %in% selected_ids
    )

    bg_idx <- (
      !is.na(bg_ids_all) &
        bg_ids_all %in% selected_ids
    )
  } else {
    pres_ids_all[is.na(pres_ids_all)] <- "<NA>"
    bg_ids_all[is.na(bg_ids_all)] <- "<NA>"

    pres_idx <- pres_ids_all %in% selected_ids
    bg_idx <- bg_ids_all %in% selected_ids
  }

  # Construct output -------------------------------------------------------

  selected_presence_proportion <- best$sumP / Np_total
  selected_background_proportion <- best$sumB / Nb_total

  out <- list(
    selected_fold_ids = selected_ids,
    pres_selected_idx = pres_idx,
    bg_selected_idx = bg_idx,
    summary = list(
      n_pres_total = length(pres_folds),
      n_bg_total = length(bg_folds),
      n_pres_in_scope = Np_total,
      n_bg_in_scope = Nb_total,
      n_pres_selected = best$sumP,
      n_bg_selected = best$sumB,
      prop_pres_selected = selected_presence_proportion,
      prop_bg_selected = selected_background_proportion,
      target_prop = target_prop,
      bounds_pres = pres_bounds,
      bounds_bg = bg_bounds,
      target_pres_abs = target_p,
      target_bg_abs_range = c(lower_b, upper_b),
      iterations = best$iter,
      score = best$score
    )
  )

  # Optional report --------------------------------------------------------

  if (print_report) {
    cat(
      "Joint fold selection (presence + background)\n",
      "--------------------------------------------\n",
      sep = ""
    )

    cat(sprintf(
      paste0(
        "Pres selected: %d / %d (%.2f%%) | ",
        "bounds [%.0f%%, %.0f%%], target %.0f%% (~%d)\n"
      ),
      out$summary$n_pres_selected,
      out$summary$n_pres_in_scope,
      100 * out$summary$prop_pres_selected,
      100 * out$summary$bounds_pres[1L],
      100 * out$summary$bounds_pres[2L],
      100 * out$summary$target_prop,
      out$summary$target_pres_abs
    ))

    cat(sprintf(
      paste0(
        "BG   selected: %d / %d (%.2f%%) | ",
        "bounds [%.0f%%, %.0f%%]\n"
      ),
      out$summary$n_bg_selected,
      out$summary$n_bg_in_scope,
      100 * out$summary$prop_bg_selected,
      100 * out$summary$bounds_bg[1L],
      100 * out$summary$bounds_bg[2L]
    ))

    cat(sprintf(
      "Prop match (BG - PRES): %+.2f percentage points\n",
      100 * (
        out$summary$prop_bg_selected -
          out$summary$prop_pres_selected
      )
    ))

    cat(sprintf(
      "Iterations used: %d | Final score: %.3f\n",
      out$summary$iterations,
      out$summary$score
    ))

    cat("Selected fold IDs:\n")

    if (length(out$selected_fold_ids) > 0L) {
      cat(
        "  ",
        paste(out$selected_fold_ids, collapse = ", "),
        "\n",
        sep = ""
      )
    } else {
      cat("  <none>\n")
    }
  }

  invisible(out)
}
