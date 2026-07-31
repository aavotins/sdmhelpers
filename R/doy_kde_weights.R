#' Calculate Day-of-Year Phenology Weights Using Circular Kernel Density
#'
#' Estimates a seasonal activity-density curve from a set of reference dates
#' and assigns phenological weights to another set of dates. Seasonality is
#' represented by day of year and estimated using a circular approximation to
#' kernel density estimation.
#'
#' @description
#' `doy_kde_weights()` converts `A_dates` and `B_dates` to day of year. A kernel
#' density curve is then estimated from `A_dates`, which are interpreted as
#' reference activity, occurrence, or phenology dates. The density of this
#' reference distribution is evaluated at each date in `B_dates`.
#'
#' Circularity is approximated by repeating the reference day-of-year values
#' one seasonal cycle before and after the original observations. This reduces
#' boundary effects between the end and beginning of the year.
#'
#' The density values associated with `B_dates` can be transformed to weights
#' using either:
#'
#'   - `"minmax_A"`: Min-max scaling relative to the complete reference density curve.
#'     The minimum density along the curve is assigned a weight of zero and
#'     the maximum density is assigned a weight of one.
#'
#'   - `"percentile_A"`: he empirical percentile rank of each density-at-`B_dates` value among
#'     all values of the reference density curve. Larger values therefore
#'     indicate dates occurring during relatively high-density portions of the
#'     estimated reference phenology.
#'
#'
#' The resulting weights may be used as relative measures of seasonal
#' availability, sampling relevance, or phenological correspondence. For
#' example, they may be supplied as observation weights in a subsequent
#' spatial kernel-density or sampling-effort analysis.
#'
#' @param A_dates A non-empty vector of reference dates representing the
#'   activity or phenological distribution from which the seasonal density
#'   curve is estimated. Values must be coercible to class \code{"Date"} by
#'   [as.Date()]. Missing or invalid dates are not allowed.
#'
#' @param B_dates A non-empty vector of dates for which phenological weights
#'   are required. Values must be coercible to class \code{"Date"} by
#'   [as.Date()]. Missing or invalid dates are not allowed.
#'
#' @param scale A character string specifying how the reference density values
#'   evaluated at `B_dates` are converted to weights. Must be one of
#'   `"minmax_A"` or `"percentile_A"`. Partial matching is supported through
#'   [match.arg()].
#'
#' @param n_days Either `NULL`, `365`, or `366`. Defines the length of the
#'   seasonal cycle used by the circular density calculation. When `NULL`,
#'   a 366-day cycle is used if any date in `A_dates` or `B_dates` falls on
#'   day 366; otherwise, a 365-day cycle is used.
#'
#'   Setting this argument explicitly can be useful when several datasets must
#'   be analysed using the same seasonal definition. Note that setting
#'   `n_days = 365` does not remove or remap observations occurring on day 366.
#'
#' @param bw Bandwidth passed to [stats::density()]. The default, `"nrd0"`,
#'   uses the corresponding automatic bandwidth-selection rule. A positive
#'   numeric bandwidth may be supplied to control the amount of seasonal
#'   smoothing directly. The bandwidth is expressed in day-of-year units.
#'
#' @param n Either `NULL` or an integer giving the number of equally spaced
#'   evaluation points used for the estimated density curve. When `NULL`,
#'   `n_days` points are used. Values smaller than 10 are not allowed.
#'
#'   Increasing `n` produces a more finely resolved density curve but does not
#'   add information beyond that contained in the input dates.
#'
#' @param eps A single numeric value in the interval
#'   \eqn{[0,\,0.5)}.
#'   Used only when `scale = "percentile_A"`.
#'   If `eps > 0`, percentile weights are clamped between
#'   `eps` and `1 - eps`.
#'   This prevents exact zero and one values before
#'   transformations such as the logit.
#'
#' @return
#' A named list containing:
#'
#'   - `weights_raw`: A numeric vector with one value per element of `B_dates`, containing the
#'     unscaled reference-density value evaluated at the corresponding day of
#'     year.
#'
#'   - `weights`: A numeric vector with one scaled phenological weight per element of
#'     `B_dates`. Values are on a zero-to-one scale, subject to clamping by
#'     `eps` when percentile scaling is used.
#'
#'   - `B_doy`: An integer vector containing the day of year derived from each element
#'     of `B_dates`.
#'
#'   - `density`: A list with numeric components `x` and `y`. Component `x` contains the
#'     day-of-year evaluation grid and component `y` contains the estimated
#'     density values derived from `A_dates`.
#'
#'   - `scale`: The scaling method used.
#'
#'   - `n_days`: The seasonal cycle length used in the calculation.
#'
#'
#' The order and length of `weights_raw`, `weights`, and `B_doy` correspond to
#' the order and length of `B_dates`.
#'
#' @details
#' Day of year is calculated as the zero-based `yday` component returned by
#' [as.POSIXlt()] plus one. Thus, January 1 is day 1 and December 31 is day 365
#' in a non-leap year or day 366 in a leap year.
#'
#' To approximate a circular density, the reference day-of-year vector
#' `A` is expanded to:
#'
#' \deqn{(A - D,\; A,\; A + D)}
#'
#' where `D` is `n_days`. A conventional one-dimensional kernel density is
#' then estimated over the interval from day 1 to day `D`. Repeating the
#' observations on both sides allows observations near the beginning of the
#' year to influence the density near the end of the year, and vice versa.
#'
#' The approach is a practical wrapped-data approximation rather than a
#' specialised circular probability-density estimator. Because three copies
#' of every reference observation are passed to [stats::density()], the
#' absolute magnitude of the returned density is affected by this
#' construction. The scaled weights remain useful for relative comparisons
#' within a result, but `weights_raw` should not be interpreted as a
#' conventional probability density integrating to one over the focal
#' day-of-year interval.
#'
#' With `"minmax_A"` scaling, weights are calculated as:
#'
#' \deqn{
#' w_i = \frac{f(B_i) - \min(f)}
#'            {\max(f) - \min(f)}
#' }
#'
#' where \eqn{f(B_i)} is the reference-density value at the day of year of the
#' \eqn{i}-th `B_dates` observation, and the minimum and maximum are calculated
#' over the full estimated reference curve.
#'
#' With `"percentile_A"` scaling, each weight is the proportion of values along
#' the estimated reference curve that are less than or equal to the
#' density-at-date value:
#'
#' \deqn{
#' w_i = \frac{1}{m}\sum_{j=1}^{m} I(f_j \leq f(B_i))
#' }
#'
#' where \eqn{m} is the number of density-grid points and
#' \eqn{\mathbf{1}\{\cdot\}} is the indicator function.
#'
#' This percentile is calculated over grid points rather than over the original
#' observations in `A_dates`. It therefore describes the relative position of
#' a date's density within the estimated annual density curve, not the
#' percentile rank of that date among the observed reference dates.
#'
#' Dates from different calendar years are pooled by day of year. Consequently,
#' the function estimates an average seasonal pattern and does not retain
#' interannual differences. In addition, calendar dates after February 28 are
#' shifted by one day in leap years relative to non-leap years because the
#' function uses literal calendar day of year. Users requiring a leap-day-free
#' or biologically standardised seasonal axis should preprocess their dates
#' before calling this function.
#'
#' @section Input validation:
#' Both date vectors must be present and non-empty. Conversion with [as.Date()]
#' must not produce missing values. The function also requires `n_days` to be
#' either 365 or 366 and `n` to be at least 10.
#'
#' For `"minmax_A"` scaling, the estimated density curve must have a finite,
#' non-degenerate range. An error is produced when its maximum is not greater
#' than its minimum.
#'
#' Validation of `eps` occurs only when `scale = "percentile_A"`, because this
#' argument is not used by min-max scaling.
#'
#' @section Bandwidth selection:
#' The bandwidth strongly affects the resulting weights. A small bandwidth
#' produces a more locally variable seasonal curve, whereas a large bandwidth
#' produces stronger smoothing across dates. Automatic bandwidth selection
#' provides a convenient default, but a biologically meaningful numeric
#' bandwidth may be preferable when the expected duration of an activity
#' period is known.
#'
#'
#' @seealso
#' [stats::density()] for kernel-density estimation,
#' [stats::approx()] for interpolation of density values, and
#' [as.Date()] for date conversion.
#'
#' @examples
#' ## Reference activity dates concentrated in spring and early summer
#' A_dates <- as.Date(c(
#'   "2015-04-28", "2015-05-07", "2015-05-18",
#'   "2016-05-03", "2016-05-21", "2016-06-04",
#'   "2017-05-11", "2017-05-26", "2017-06-09"
#' ))
#'
#' ## Inventory dates to be weighted
#' B_dates <- as.Date(c(
#'   "2018-03-15",
#'   "2018-05-15",
#'   "2018-06-15",
#'   "2018-08-15"
#' ))
#'
#' ## Min-max weights relative to the complete reference curve
#' res_minmax <- doy_kde_weights(
#'   A_dates = A_dates,
#'   B_dates = B_dates,
#'   scale = "minmax_A"
#' )
#'
#' res_minmax$weights
#' res_minmax$B_doy
#'
#' ## Percentile weights, excluding exact zero and one
#' res_percentile <- doy_kde_weights(
#'   A_dates = A_dates,
#'   B_dates = B_dates,
#'   scale = "percentile_A",
#'   eps = 0.01
#' )
#'
#' res_percentile$weights
#'
#' ## Combine dates and calculated weights
#' data.frame(
#'   inventory_date = B_dates,
#'   day_of_year = res_percentile$B_doy,
#'   raw_density = res_percentile$weights_raw,
#'   phenology_weight = res_percentile$weights
#' )
#'
#' ## Use a biologically selected bandwidth of 14 days
#' res_bw <- doy_kde_weights(
#'   A_dates = A_dates,
#'   B_dates = B_dates,
#'   scale = "minmax_A",
#'   bw = 14
#' )
#'
#' ## Inspect the estimated seasonal density curve
#' plot(
#'   res_bw$density$x,
#'   res_bw$density$y,
#'   type = "l",
#'   xlab = "Day of year",
#'   ylab = "Estimated reference density"
#' )
#' points(
#'   res_bw$B_doy,
#'   res_bw$weights_raw,
#'   pch = 19
#' )
#'
#' ## Explicitly use a 366-day cycle
#' leap_result <- doy_kde_weights(
#'   A_dates = as.Date(c("2020-02-20", "2020-02-29", "2020-03-08")),
#'   B_dates = as.Date(c("2020-02-29", "2020-03-05")),
#'   n_days = 366,
#'   bw = 5
#' )
#'
#' leap_result$n_days
#'
#' @export
doy_kde_weights <- function(
    A_dates,
    B_dates,
    scale = c("minmax_A", "percentile_A"),
    n_days = NULL,
    bw = "nrd0",
    n = NULL,
    eps = 1e-6
) {
  scale <- match.arg(scale)

  if (missing(A_dates) || is.null(A_dates) || length(A_dates) == 0L) {
    stop("`A_dates` must be provided and non-empty.")
  }
  if (missing(B_dates) || is.null(B_dates) || length(B_dates) == 0L) {
    stop("`B_dates` must be provided and non-empty.")
  }

  A_dates <- as.Date(A_dates)
  B_dates <- as.Date(B_dates)
  if (anyNA(A_dates)) stop("`A_dates` contains NA after coercion to Date.")
  if (anyNA(B_dates)) stop("`B_dates` contains NA after coercion to Date.")

  doy <- function(x) as.POSIXlt(x)$yday + 1L

  A_doy <- doy(A_dates)
  B_doy <- doy(B_dates)

  if (is.null(n_days)) {
    n_days <- if (any(A_doy == 366L) || any(B_doy == 366L)) 366L else 365L
  } else {
    n_days <- as.integer(n_days)
  }
  if (!isTRUE(n_days %in% c(365L, 366L))) {
    stop("`n_days` must be 365 or 366.")
  }

  if (is.null(n)) n <- n_days
  n <- as.integer(n)
  if (n < 10L) stop("`n` must be >= 10.")

  # Wrap-around trick for circular KDE on DoY
  A_ext <- c(A_doy - n_days, A_doy, A_doy + n_days)

  dens <- stats::density(A_ext, bw = bw, from = 1, to = n_days, n = n)

  # Look up A-derived density at DoY(B)
  w_raw <- stats::approx(dens$x, dens$y, xout = B_doy, rule = 2)$y

  if (scale == "minmax_A") {
    dmin <- min(dens$y)
    dmax <- max(dens$y)
    if (!is.finite(dmin) || !is.finite(dmax) || dmax <= dmin) {
      stop("Density curve has non-finite or degenerate range; cannot min-max scale.")
    }
    w <- (w_raw - dmin) / (dmax - dmin)
    w <- pmin(pmax(w, 0), 1)
  } else if (scale == "percentile_A") {
    d_all <- dens$y
    w <- vapply(w_raw, function(v) mean(d_all <= v), numeric(1))
    if (!is.numeric(eps) || length(eps) != 1L || eps < 0 || eps >= 0.5) {
      stop("`eps` must be a single numeric value in [0, 0.5).")
    }
    if (eps > 0) w <- pmin(pmax(w, eps), 1 - eps)
  } else {
    stop("Unknown `scale` method.")
  }

  list(
    weights_raw = w_raw,
    weights     = w,
    B_doy       = B_doy,
    density     = list(x = dens$x, y = dens$y),
    scale       = scale,
    n_days      = n_days
  )
}
