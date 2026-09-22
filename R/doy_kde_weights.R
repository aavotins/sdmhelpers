#' Calculate Day-of-Year Phenology Weights Using Circular Approximation to Kernel Density Estimation
#'
#' Estimates a seasonal density curve from a set of reference dates
#' and assigns phenological weights to another set of dates. Seasonality is
#' represented by day of year and estimated using a circular approximation to
#' kernel density estimation.
#'
#' @description
#' `doy_kde_weights()` converts `A_dates` and `B_dates` to day of year. A kernel
#' density curve is then estimated from `A_dates`, which are interpreted as
#' reference activity, occurrence, or phenology dates. The estimated density of this
#' reference distribution is evaluated at each date in `B_dates`.
#'
#' Circularity is approximated by periodically extending the reference day-of-year
#' data, repeating the observations one seasonal cycle before and after the original
#' observations. This approach reduces boundary effects between the end and
#' beginning of the year.
#'
#' The density values associated with `B_dates` can be transformed into weights
#' using either of the following methods:
#'
#'   - `"minmax_A"`: Min-max scaling relative to the complete reference density curve.
#'     The minimum density value along the curve is assigned a weight of zero and
#'     the maximum density value is assigned a weight of one.
#'
#'   - `"percentile_A"`: The empirical percentile rank of each density-at-`B_dates` value among
#'     the density values evaluated over the complete reference-grid. Larger values therefore
#'     indicate dates occurring during relatively high-density portions of the
#'     estimated reference phenology.
#'
#'
#' The resulting weights are relative measures of seasonal correspondence to the reference density.
#' Depending on the interpretation of `A_dates`, they may be used to represent relative seasonal
#' availability, sampling relevance, or phenological correspondence. For
#' example, the weights may be supplied as observation weights in a subsequent
#' spatial kernel-density or sampling-effort analysis.
#'
#' @param A_dates A non-empty vector of reference dates representing the
#'   seasonal activity or phenology from which the seasonal density
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
#'   be analysed using the same seasonal definition. If
#'   `n_days = 365`, neither `A_dates` nor `B_dates` may contain an observation
#'   on day 366; otherwise, the function returns an error.
#'
#' @param bw Bandwidth used for kernel density estimation. May be a positive
#'   numeric value expressed in day-of-year units or one of `"nrd0"`, `"nrd"`,
#'   `"ucv"`, `"bcv"`, or `"SJ"` for automatic bandwidth selection.
#'   The default is `"nrd0"`.
#'
#'   For automatic bandwidth selection, the bandwidth is estimated from the original
#'   reference day-of-year values (`A_doy`) using the corresponding function in `stats`.
#'   This calculation is performed before the reference observations are repeated for
#'   the circular approximation. The resulting numeric bandwidth is then passed
#'   to [stats::density()]
#'
#'
#' @param n Either `NULL` or an integer giving the number of equally spaced
#'   evaluation points used to represent the estimated density curve. When `NULL`,
#'   `n_days` points are used. Values smaller than 10 are not allowed.
#'
#'   Increasing `n` produces a more finely resolved evaluation grid but does not
#'   add information beyond that contained in the input dates or change the
#'   underlying bandwidth.
#'
#' @param eps A single numeric value in the interval
#'   \eqn{[0,\,0.5)}.
#'   Used only when `scale = "percentile_A"`.
#'   If `eps > 0`, percentile weights are clamped between
#'   to the interval \eqn{[0,\,1 - eps]}.
#'   This prevents exact zero and one values before
#'   transformations such as the logit.
#'
#' @return
#' A named list containing:
#'
#'   - `weights_raw`: A numeric vector with one value per element of `B_dates`, containing the
#'     unscaled estimated reference-density value evaluated at the corresponding day of
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
#' specialised circular probability-density estimator. Because each reference
#' day-of-year value is represented three times as periodically shifted copies
#' in the input to [stats::density()], the resulting kernel density is normalised
#' over the full replicated input rather than over a single seasonal cycle.
#' Consequently, the absolute magnitude of the density over the focal day-of-year
#' interval depends on this construction. The scaled weights remain useful for
#' relative comparisons within a result, but `weights_raw` should not be
#' interpreted as a conventional probability density integrating to one over the
#' focal day-of-year interval.
#'
#' @section Min-max scaling:
#' With `"minmax_A"` scaling, weights are calculated as:
#'
#' \deqn{
#' w_i = \frac{f(B_i) - \min(f)}
#'            {\max(f) - \min(f)}
#' }
#'
#' where \eqn{f(B_i)} is the reference-density value at the day of year of the
#' \eqn{i}-th `B_dates` observation, and the minimum and maximum are calculated
#' over the complete estimated reference density curve.
#'
#' @section Percentile scaling:
#' With `"percentile_A"` scaling, each weight represents the proportion of
#' evaluation-grid points at which the estimated reference density is less than
#' or equal to the density at the corresponding `B_dates` observation:
#'
#' \deqn{
#' w_i = \frac{1}{m}\sum_{j=1}^{m} I(f_j \leq f(B_i))
#' }
#'
#' where \eqn{m} is the number of density-grid points and
#' \eqn{\mathbf{I}\{\cdot\}} is the indicator function.
#'
#' This weight is calculated from the density value at the evaluation-grid points
#' rather than over the original observations in `A_dates`. It therefore describes
#' how the density at a given date ranks relative to the estimated density levels
#' across the seasonal cycle,  not the percentile rank of that date among the
#' observed reference dates.
#'
#' Dates from different calendar years are pooled by day of year. Consequently,
#' the function estimates a seasonal pattern pooled across years and does not retain
#' interannual differences. The function uses literal calendar day of year.
#' Therefore, for dates from March 1 onward, the same month-and-day date has a
#' day-of-year value one greater in a leap year than in a non-leap year. Users
#' requiring a leap-day-free or biologically standardised seasonal axis should
#' preprocess their dates before calling this function.
#'
#' @section Input validation:
#' Both date vectors must be provided and non-empty. Conversion with [as.Date()]
#' must not produce missing values. The function also requires `n_days` to be
#' either 365 or 366 and `n` to be at least 10.
#'
#' If n_days = 365, neither date vector may contain an observation on day 366.
#' Setting n_days = 366 is allowed even when no observation falls on day 366;
#' in that case, the density is still estimated on a 366-day seasonal cycle.
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
#' Automatic bandwidth selection is performed using the original reference
#' day-of-year values before observations are repeated for the circular
#' approximation, so the repeated observations used for boundary correction
#' do not influence bandwidth selection.
#'
#' Because these automatic bandwidth selectors are conventional one-dimensional
#' methods, they do not explicitly account for the circular distance between
#' the beginning and end of the year.
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

  if (anyNA(A_dates)) {
    stop("`A_dates` contains NA after coercion to Date.")
  }
  if (anyNA(B_dates)) {
    stop("`B_dates` contains NA after coercion to Date.")
  }

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

  if (
    n_days == 365L &&
    (any(A_doy == 366L) || any(B_doy == 366L))
  ) {
    stop(
      "`n_days = 365` cannot be used when `A_dates` or `B_dates` contains day 366."
    )
  }

  if (is.null(n)) {
    n <- n_days
  }

  n <- as.integer(n)

  if (n < 10L) {
    stop("`n` must be >= 10.")
  }

  # Unwrap circular DoY data for automatic bandwidth estimation.
  #
  # Standard bandwidth selectors treat the data as linear. Therefore,
  # observations around the end and beginning of the year (e.g. December
  # and January) would otherwise appear far apart. The circle is cut at
  # the largest observed gap, and observations before the cut are shifted
  # forward by one year.
  unwrap_doy_for_bw <- function(x, n_days) {
    x <- sort(as.numeric(x))

    if (length(x) < 2L) {
      return(x)
    }

    gaps <- c(
      diff(x),
      x[1L] + n_days - x[length(x)]
    )

    k <- which.max(gaps)

    # If the largest gap is already across the year boundary,
    # no unwrapping is necessary.
    if (k == length(x)) {
      return(x)
    }

    cut <- x[k + 1L]

    x[x < cut] <- x[x < cut] + n_days

    sort(x)
  }

  # Data used only for automatic bandwidth estimation
  A_bw <- unwrap_doy_for_bw(A_doy, n_days)

  bw_used <- if (is.character(bw)) {
    switch(
      bw,
      nrd0 = stats::bw.nrd0(A_bw),
      nrd  = stats::bw.nrd(A_bw),
      ucv  = stats::bw.ucv(A_bw),
      bcv  = stats::bw.bcv(A_bw),
      SJ   = stats::bw.SJ(A_bw),
      stop("Unknown bandwidth selection method.")
    )
  } else {
    bw
  }

  # Wrap-around trick for circular KDE on DoY.
  #
  # The replicated observations ensure that density from observations
  # near one end of the year contributes to density at the other end.
  A_ext <- c(
    A_doy - n_days,
    A_doy,
    A_doy + n_days
  )

  dens <- stats::density(
    A_ext,
    bw = bw_used,
    from = 1,
    to = n_days,
    n = n
  )

  # Look up A-derived density at DoY(B)
  w_raw <- stats::approx(
    dens$x,
    dens$y,
    xout = B_doy,
    rule = 2
  )$y

  if (scale == "minmax_A") {

    dmin <- min(dens$y)
    dmax <- max(dens$y)

    if (!is.finite(dmin) || !is.finite(dmax) || dmax <= dmin) {
      stop(
        "Density curve has non-finite or degenerate range; cannot min-max scale."
      )
    }

    w <- (w_raw - dmin) / (dmax - dmin)
    w <- pmin(pmax(w, 0), 1)

  } else if (scale == "percentile_A") {

    d_all <- dens$y

    w <- vapply(
      w_raw,
      function(v) mean(d_all <= v),
      numeric(1)
    )

    if (
      !is.numeric(eps) ||
      length(eps) != 1L ||
      eps < 0 ||
      eps >= 0.5
    ) {
      stop("`eps` must be a single numeric value in [0, 0.5).")
    }

    if (eps > 0) {
      w <- pmin(pmax(w, eps), 1 - eps)
    }

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
