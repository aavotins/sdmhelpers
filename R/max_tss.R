#' Calculate the Maximum True Skill Statistic
#'
#' Calculates the maximum True Skill Statistic (TSS) across all unique
#' prediction thresholds for presence and background or absence predictions.
#'
#' @param p_pos A numeric vector containing predicted suitability values for
#'   observed presences.
#' @param p_neg A numeric vector containing predicted suitability values for
#'   background records or observed absences.
#'
#' @details
#' For a given threshold, predictions greater than or equal to the threshold are
#' classified as positive, whereas predictions below the threshold are
#' classified as negative.
#'
#' Sensitivity is the proportion of observed presences classified as positive:
#'
#' \deqn{
#' \mathrm{sensitivity} =
#' \frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FN}},
#' }
#'
#' where TP is the number of presence records with predictions greater than or
#' equal to the threshold and FN is the number of presence records with
#' predictions below the threshold.
#'
#' Specificity is calculated from the records supplied in `p_neg` as
#'
#' \deqn{
#' \mathrm{specificity} =
#' \frac{\mathrm{TN}}{\mathrm{TN} + \mathrm{FP}},
#' }
#'
#' where TN is the number of records in `p_neg` with predictions below the
#' threshold and FP is the number with predictions greater than or equal to the
#' threshold.
#'
#' When `p_neg` contains observed absences, TN and FP have their usual
#' classification interpretation. When `p_neg` contains background records,
#' however, those records are not necessarily true absences. In that case,
#' "specificity" should be interpreted operationally as the proportion of
#' background records classified as negative rather than as the proportion of
#' confirmed absences correctly classified. Consequently, the resulting TSS is
#' a presence-background discrimination measure rather than a TSS based on
#' independently observed true absences.
#'
#' The True Skill Statistic is calculated as
#'
#' \deqn{
#' \mathrm{TSS} =
#' \mathrm{sensitivity} + \mathrm{specificity} - 1.
#' }
#'
#' Thresholds are evaluated at every unique finite prediction value. Two
#' additional classifications are also considered:
#'
#' \itemize{
#' \item all records classified as negative; and
#' \item all records classified as positive.
#' }
#'
#' Both extreme classifications have a TSS of zero.
#'
#' Non-finite values, including `NA`, `NaN`, `Inf`, and `-Inf`, are removed
#' separately from `p_pos` and `p_neg.` If no finite values remain in either
#' vector, all returned values are `NA_real_`.
#'
#' When multiple thresholds produce the same maximum TSS, the highest threshold
#' is returned. This corresponds to the most restrictive classification among
#' the tied thresholds. One exception occurs when no evaluated threshold has a
#' TSS greater than zero: the function returns `Inf`, representing
#' classification of all records as negative.
#'
#' TSS is commonly described as being independent of prevalence because
#' sensitivity and specificity are calculated conditionally within the
#' positive and negative samples. Nevertheless, the value obtained here can
#' depend strongly on how `p_pos` and especially `p_neg` were sampled. This is
#' particularly important when `p_neg` represents background rather than
#' observed absences, because the composition, extent, and sampling design of
#' the background data affect the estimated specificity and therefore the
#' resulting TSS.
#'
#' @return
#' A named numeric vector of length four containing:
#'
#' - `tss`: Maximum True Skill Statistic.
#' - `threshold`: Threshold producing the maximum TSS.
#' - `sensitivity`: Sensitivity at the selected threshold.
#' - `specificity`: Specificity at the selected threshold.
#'
#' @seealso
#' [tss_user_eval()]
#'
#' @references
#' Allouche, O., Tsoar, A., and Kadmon, R. (2006). Assessing the accuracy of
#' species distribution models: prevalence, kappa and the true skill statistic
#' (TSS). *Journal of Applied Ecology*, **43**, 1223--1232.
#' [doi:10.1111/j.1365-2664.2006.01214.x](https://doi.org/10.1111/j.1365-2664.2006.01214.x)
#'
#' @examples
#' presence_predictions <- c(0.90, 0.85, 0.75, 0.60, 0.40)
#' background_predictions <- c(0.70, 0.50, 0.30, 0.20, 0.10)
#'
#' max_tss(
#'   p_pos = presence_predictions,
#'   p_neg = background_predictions
#' )
#'
#' # Non-finite values are removed
#' max_tss(
#'   p_pos = c(0.9, 0.8, NA, Inf),
#'   p_neg = c(0.4, 0.2, NA)
#' )
#'
#' # Apply the selected threshold manually
#' result <- max_tss(
#'   p_pos = presence_predictions,
#'   p_neg = background_predictions
#' )
#'
#' threshold <- unname(result["threshold"])
#' presence_classes <- presence_predictions >= threshold
#' background_classes <- background_predictions >= threshold
#'
#' @export
max_tss <- function(p_pos, p_neg) {
  if (!is.numeric(p_pos)) {
    stop("`p_pos` must be a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(p_neg)) {
    stop("`p_neg` must be a numeric vector.", call. = FALSE)
  }

  p_pos <- p_pos[is.finite(p_pos)]
  p_neg <- p_neg[is.finite(p_neg)]

  if (length(p_pos) == 0L || length(p_neg) == 0L) {
    return(
      c(
        tss = NA_real_,
        threshold = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_
      )
    )
  }

  pred <- c(p_pos, p_neg)

  obs <- c(
    rep.int(1L, length(p_pos)),
    rep.int(0L, length(p_neg))
  )

  # Sort from highest to lowest prediction.
  ordered <- order(pred, decreasing = TRUE)

  pred <- pred[ordered]
  obs <- obs[ordered]

  n_positive <- length(p_pos)
  n_negative <- length(p_neg)

  true_positive <- cumsum(obs == 1L)
  false_positive <- cumsum(obs == 0L)

  # Use the last position of each tied prediction value. This ensures that
  # all records with predictions equal to the threshold are classified in
  # the same way under the rule prediction >= threshold.
  cutpoint_index <- c(
    which(diff(pred) != 0),
    length(pred)
  )

  sensitivity <- true_positive[cutpoint_index] / n_positive

  specificity <- (
    n_negative - false_positive[cutpoint_index]
  ) / n_negative

  tss <- sensitivity + specificity - 1

  # Inf represents predicting no records as positive.
  # -Inf represents predicting every record as positive.
  all_tss <- c(
    0,
    tss,
    0
  )

  all_thresholds <- c(
    Inf,
    pred[cutpoint_index],
    -Inf
  )

  all_sensitivity <- c(
    0,
    sensitivity,
    1
  )

  all_specificity <- c(
    1,
    specificity,
    0
  )

  # which.max() returns the first maximum. Because thresholds are arranged
  # from highest to lowest, ties are resolved in favour of the highest
  # threshold.
  best <- which.max(all_tss)

  c(
    tss = all_tss[best],
    threshold = all_thresholds[best],
    sensitivity = all_sensitivity[best],
    specificity = all_specificity[best]
  )
}
