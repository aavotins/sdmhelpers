#' Calculate phenological overlap weights among species
#'
#' Calculates the overlap between the day-of-year distribution of a target
#' species and the day-of-year distributions of multiple comparison species.
#' Kernel-density overlap is calculated separately for each comparison species
#' using [overlapping::overlap()].
#'
#' The resulting overlap estimates are additionally normalized across all
#' comparison species so that their weights sum to 1.
#'
#' @param target A data frame containing observations of one or more potential
#'   target species.
#' @param group A data frame containing observations of the comparison species.
#' @param target_code Species code identifying the target species in `target`.
#' @param target_id Character string giving the species-code column in
#'   `target`. Defaults to `"CODE"`.
#' @param group_id Character string giving the species-code column in
#'   `group`. Defaults to `"CODE"`.
#' @param target_doy Character string giving the day-of-year column in
#'   `target`. Defaults to `"DoY"`.
#' @param group_doy Character string giving the day-of-year column in
#'   `group`. Defaults to `"DoY"`.
#' @param type Character string specifying the overlap measure passed to
#'   [overlapping::overlap()]. Either `"1"` or `"2"`. Defaults to `"2"`.
#' @param nbins Integer giving the number of equally spaced points used for
#'   comparing the estimated densities. Passed to [overlapping::overlap()].
#'   Defaults to `1024`.
#' @param boundaries Optional numeric vector of length two defining the lower
#'   and upper boundaries over which overlap is evaluated. Defaults to `NULL`.
#' @param min_n Minimum number of finite day-of-year observations required for
#'   the target species and each comparison species. Defaults to `2`.
#' @param na_rm Logical. If `TRUE`, non-finite day-of-year values are removed
#'   before estimating densities. Defaults to `TRUE`.
#' @param ... Additional arguments passed to [overlapping::overlap()] and
#'   subsequently to [stats::density()], for example `bw`, `adjust`, or
#'   `kernel`.
#'
#' @return A data frame containing the comparison-species code, the estimated
#'   overlap with the target species, and the normalized overlap weight.
#'
#' @details
#' The target species is selected from `target` using `target_code` and the
#' column specified by `target_id`.
#'
#' Each species in `group` is then compared independently with the target
#' species. The normalized weight for species \eqn{i} is
#'
#' \deqn{w_i = O_i / \sum_j O_j}
#'
#' where \eqn{O_i} is the estimated overlap between the target species and
#' comparison species \eqn{i}.
#'
#' Missing overlap estimates are excluded from the denominator. If all valid
#' overlap estimates are zero, weights are returned as `NA`.
#'
#' Note that [overlapping::overlap()] uses ordinary linear kernel-density
#' estimation. Day of year is circular, so observations near the beginning
#' and end of the year are not treated as adjacent.
#'
#' @seealso [overlapping::overlap()], [stats::density()]
#'
#' @references
#' Pastore, M. (2018). Overlapping: an R package for estimating overlapping in
#' empirical distributions. Journal of Open Source Software, 3(32), 1023.
#' \doi{10.21105/joss.01023}
#'
#' Pastore, M., & Calcagni, A. (2019). Measuring distribution similarities
#' between samples: A distribution-free overlapping index.
#' Frontiers in Psychology, 10, 1089.
#' \doi{10.3389/fpsyg.2019.01089}
#'
#' @examples
#' set.seed(123)
#'
#' target <- data.frame(
#'   CODE = rep(c("TARGET1", "TARGET2"), each = 100),
#'   DoY = c(
#'     round(rnorm(100, 140, 20)),
#'     round(rnorm(100, 220, 20))
#'   )
#' )
#'
#' group <- data.frame(
#'   CODE = rep(c("SP1", "SP2", "SP3"), each = 100),
#'   DoY = c(
#'     round(rnorm(100, 145, 20)),
#'     round(rnorm(100, 190, 20)),
#'     round(rnorm(100, 250, 20))
#'   )
#' )
#'
#' phenology_overlap_weights(
#'   target = target,
#'   group = group,
#'   target_code = "TARGET1"
#' )
#'
#' @export
phenology_overlap_weights <- function(
    target,
    group,
    target_code,
    target_id = "CODE",
    group_id = "CODE",
    target_doy = "DoY",
    group_doy = "DoY",
    type = "2",
    nbins = 1024,
    boundaries = NULL,
    min_n = 2L,
    na_rm = TRUE,
    ...
) {

  # Check data frames
  if (!is.data.frame(target)) {
    stop("`target` must be a data frame.", call. = FALSE)
  }

  if (!is.data.frame(group)) {
    stop("`group` must be a data frame.", call. = FALSE)
  }

  # Check required columns
  if (!target_id %in% names(target)) {
    stop(
      sprintf("Column '%s' not found in `target`.", target_id),
      call. = FALSE
    )
  }

  if (!target_doy %in% names(target)) {
    stop(
      sprintf("Column '%s' not found in `target`.", target_doy),
      call. = FALSE
    )
  }

  if (!group_id %in% names(group)) {
    stop(
      sprintf("Column '%s' not found in `group`.", group_id),
      call. = FALSE
    )
  }

  if (!group_doy %in% names(group)) {
    stop(
      sprintf("Column '%s' not found in `group`.", group_doy),
      call. = FALSE
    )
  }

  # Check target code
  if (length(target_code) != 1L || is.na(target_code)) {
    stop("`target_code` must contain exactly one non-missing value.",
         call. = FALSE)
  }

  if (!target_code %in% target[[target_id]]) {
    stop(
      sprintf(
        "Target species '%s' not found in column '%s' of `target`.",
        target_code,
        target_id
      ),
      call. = FALSE
    )
  }

  # Other arguments
  type <- match.arg(type, c("1", "2"))

  min_n <- as.integer(min_n)

  if (length(min_n) != 1L || is.na(min_n) || min_n < 2L) {
    stop("`min_n` must be an integer >= 2.", call. = FALSE)
  }

  if (!is.numeric(target[[target_doy]])) {
    stop("The target day-of-year column must be numeric.", call. = FALSE)
  }

  if (!is.numeric(group[[group_doy]])) {
    stop("The group day-of-year column must be numeric.", call. = FALSE)
  }

  # Select target species
  target_values <- target[
    !is.na(target[[target_id]]) &
      target[[target_id]] == target_code,
    target_doy
  ]

  # Remove non-finite target values
  if (na_rm) {

    target_values <- target_values[is.finite(target_values)]

  } else if (any(!is.finite(target_values))) {

    stop(
      "Non-finite values found in the target DoY column. ",
      "Use `na_rm = TRUE` to remove them.",
      call. = FALSE
    )
  }

  if (length(target_values) < min_n) {
    stop(
      sprintf(
        "Target species '%s' contains only %d usable observations; at least %d are required.",
        target_code,
        length(target_values),
        min_n
      ),
      call. = FALSE
    )
  }

  # Comparison species
  ids <- unique(group[[group_id]])
  ids <- ids[!is.na(ids)]

  if (length(ids) == 0L) {
    stop(
      sprintf("No valid groups found in column '%s'.", group_id),
      call. = FALSE
    )
  }

  # Calculate overlap for each comparison species
  overlap_values <- vapply(
    ids,
    FUN = function(id) {

      idx <- !is.na(group[[group_id]]) &
        group[[group_id]] == id

      values <- group[[group_doy]][idx]

      if (na_rm) {

        values <- values[is.finite(values)]

      } else if (any(!is.finite(values))) {

        return(NA_real_)
      }

      if (length(values) < min_n) {
        return(NA_real_)
      }

      result <- overlapping::overlap(
        list(
          target = target_values,
          comparison = values
        ),
        nbins = nbins,
        type = type,
        pairsOverlap = TRUE,
        plot = FALSE,
        boundaries = boundaries,
        ...
      )

      as.numeric(result$OV)
    },
    FUN.VALUE = numeric(1)
  )

  # Normalize overlap values
  total_overlap <- sum(overlap_values, na.rm = TRUE)

  if (!is.finite(total_overlap) || total_overlap <= 0) {

    weights <- rep(NA_real_, length(overlap_values))

    warning(
      "Overlap values could not be normalized because their sum is zero ",
      "or non-finite.",
      call. = FALSE
    )

  } else {

    weights <- overlap_values / total_overlap
  }

  # Output
  result <- data.frame(
    ids,
    overlap = overlap_values,
    weight = weights,
    stringsAsFactors = FALSE
  )

  names(result)[1] <- group_id

  result
}
