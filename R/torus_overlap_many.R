#' Estimate Seasonal and Daily Activity Overlap on a Torus
#'
#' @description
#' Estimates the joint seasonal and daily (registration) activity density of a
#' reference dataset and compares it with the corresponding densities of one or
#' more groups supplied in a separate comparison dataset.
#'
#' Both calendar date and time of day are treated as circular variables:
#' December is adjacent to January, and the end of one day is adjacent to the
#' beginning of the next. The resulting sample space is therefore a
#' two-dimensional torus.
#'
#' A product-kernel density estimator based on two von Mises kernels is used.
#' Separate bandwidths are specified for the seasonal and daily dimensions.
#' The overlap coefficient between the reference density and each comparison
#' density is calculated by numerical integration of the pointwise minimum of
#' the two densities.
#'
#' @param reference_data A data frame containing the observations defining the
#'   reference activity distribution. It must contain a date column of class
#'   [Date] and a time-of-day column inheriting from class `"hms"`.
#'
#'   A grouping column is not required because all valid rows in
#'   `reference_data` are treated as observations from the same reference
#'   group.
#'
#' @param comparison_data A data frame containing observations for one or more
#'   comparison groups. It must contain:
#'
#'   * a grouping column;
#'   * a date column of class [Date]; and
#'   * a time-of-day column inheriting from class `"hms"`.
#'
#'   The reference distribution is compared separately with every unique,
#'   non-missing value in the grouping column.
#'
#' @param comparison_group_col A character string giving the name of the
#'   grouping column in `comparison_data`.
#'
#' @param reference_date_col A character string giving the name of the
#'   [Date] column in `reference_data`.
#'
#' @param reference_time_col A character string giving the name of the
#'   `"hms"` time-of-day column in `reference_data`.
#'
#' @param comparison_date_col A character string giving the name of the
#'   [Date] column in `comparison_data`. By default, the same column name as
#'   `reference_date_col` is used.
#'
#' @param comparison_time_col A character string giving the name of the
#'   `"hms"` time-of-day column in `comparison_data`. By default, the same
#'   column name as `reference_time_col` is used.
#'
#' @param reference_name A single value used to identify the reference dataset
#'   in the returned tables. The value is converted to character. The default
#'   is `"reference"`.
#'
#' @param season_bw_days A positive numeric value giving the seasonal kernel
#'   bandwidth in days. Larger values produce stronger smoothing across the
#'   annual cycle. The default is `14`.
#'
#' @param daily_bw_hours A positive numeric value giving the daily kernel
#'   bandwidth in hours. Larger values produce stronger smoothing across the
#'   daily cycle. The default is `1`.
#'
#' @param date_resolution_days A positive numeric value giving the requested
#'   approximate spacing, in days, between evaluation points on the seasonal
#'   grid. This controls numerical integration resolution, not kernel
#'   smoothing. The default is `2`.
#'
#' @param time_resolution_minutes A positive numeric value giving the requested
#'   approximate spacing, in minutes, between evaluation points on the daily
#'   grid. This controls numerical integration resolution, not kernel
#'   smoothing. The default is `10`.
#'
#' @param min_n A numeric value giving the minimum number of valid observations
#'   required for density estimation. It is converted to an integer and must
#'   be at least `2`.
#'
#'   If the reference dataset contains fewer than `min_n` valid observations,
#'   the function stops with an error. Comparison groups with fewer than
#'   `min_n` observations are retained in the overlap table, but their overlap
#'   estimate is returned as `NA`.
#'
#' @param return_density Logical. If `FALSE`, the default, only the compact
#'   overlap and settings tables are returned. If `TRUE`, a potentially large
#'   long-format data frame containing grid-level density values is also
#'   returned.
#'
#' @details
#' ## Circular transformation
#'
#' Calendar dates are converted to seasonal angles. For observation \eqn{i},
#' the seasonal angle is
#'
#' \deqn{
#' \theta_i =
#' 2\pi
#' \frac{d_i - 1}{L_i},
#' }
#'
#' where \eqn{d_i} is the day of year and \eqn{L_i} is the length of the
#' observation year: either 365 or 366 days.
#'
#' Consequently, leap-year observations are scaled according to the actual
#' length of their calendar year.
#'
#' Time of day is converted to a daily angle using
#'
#' \deqn{
#' \phi_i =
#' 2\pi
#' \frac{s_i}{86400},
#' }
#'
#' where \eqn{s_i} is the number of seconds since midnight.
#'
#' ## Kernel density estimation
#'
#' The joint density is estimated using a product of two von Mises kernels:
#'
#' \deqn{
#' \widehat{f}(\theta,\phi) =
#' \frac{1}{n}
#' \sum_{i=1}^{n}
#' K_{\kappa_s}(\theta-\theta_i)
#' K_{\kappa_t}(\phi-\phi_i),
#' }
#'
#' where \eqn{\kappa_s} and \eqn{\kappa_t} are the seasonal and daily
#' concentration parameters.
#'
#' User-supplied bandwidths are converted from days and hours to approximate
#' circular standard deviations in radians. Concentration is then approximated
#' as
#'
#' \deqn{
#' \kappa \approx \frac{1}{\sigma^2}.
#' }
#'
#' This approximation is most appropriate for reasonably concentrated
#' kernels. Very large bandwidths correspond to low concentration and should
#' be interpreted cautiously.
#'
#' The von Mises kernel is evaluated using an exponentially scaled modified
#' Bessel function. This avoids numerical overflow for large concentration
#' parameters.
#'
#' ## Overlap coefficient
#'
#' The overlap coefficient between the reference density
#' \eqn{\widehat{f}_r} and comparison density \eqn{\widehat{f}_c} is
#'
#' \deqn{
#' \Delta =
#' \int_0^{2\pi}
#' \int_0^{2\pi}
#' \min\{
#' \widehat{f}_r(\theta,\phi),
#' \widehat{f}_c(\theta,\phi)
#' \}
#' \,d\theta\,d\phi.
#' }
#'
#' It ranges from zero to one:
#'
#' * `0` indicates no estimated overlap;
#' * `1` indicates identical estimated densities.
#'
#' The integral is approximated on a regular toroidal grid. Decreasing
#' `date_resolution_days` or `time_resolution_minutes` increases numerical
#' resolution but also increases computing time and memory use.
#'
#' ## Missing and invalid observations
#'
#' Rows are removed from `reference_data` when their date or time is missing,
#' non-finite, or invalid.
#'
#' Rows are removed from `comparison_data` when:
#'
#' * the comparison group is missing or an empty string;
#' * the date is missing; or
#' * the time is missing, non-finite, or outside the interval from midnight
#'   inclusive to 24:00 exclusive.
#'
#' Counts of removed observations are included in the returned settings table.
#'
#' ## Optional density output
#'
#' When `return_density = TRUE`, density values are returned in a long data
#' frame rather than as matrices. Each valid comparison group contributes
#'
#' \deqn{
#' n_{\mathrm{date}} \times n_{\mathrm{time}}
#' }
#'
#' rows.
#'
#' The reference density is repeated for each comparison group to make each
#' group independently usable for plotting and filtering. This output can
#' therefore become large when many groups or a fine evaluation grid are used.
#'
#' Temporary kernel and density matrices are still created internally to
#' perform the calculations, but they are not returned.
#'
#' @return
#' An object of class `"torus_overlap_many"` and `"list"` containing:
#'
#'   - overlap: A data frame with one row per comparison group and the following
#'     columns:
#'
#'       - reference: Reference identifier supplied through `reference_name`.
#'
#'       - comparison: Comparison-group identifier.
#'
#'       - n_reference: Number of valid observations in the reference dataset.
#'
#'       - n_comparison: Number of valid observations in the comparison group.
#'
#'       - overlap: stimated overlap coefficient. This is `NA` when the comparison
#'         group contains fewer than `min_n` observations.
#'
#'       - status: `"OK"` for successfully estimated comparisons or a message
#'         explaining why an overlap value was not calculated.
#'
#'   - settings: A one-row data frame describing column mappings, bandwidths,
#'     concentration parameters, requested and realised grid resolutions,
#'     grid dimensions, sample-size requirements, and numbers of valid and
#'     removed observations.
#'
#'   - density: Present only when `return_density = TRUE`. A long-format data frame
#'     containing:
#'
#'       - reference: Reference identifier.
#'
#'       - comparison: Comparison-group identifier.
#'
#'       - date_angle: Seasonal grid-cell centre in radians.
#'
#'       - time_angle: Daily grid-cell centre in radians.
#'
#'       - seasonal_position_days: Seasonal grid-cell centre expressed as days from the beginning of a
#'         common 365.2425-day year.
#'
#'       - time_minutes: Daily grid-cell centre expressed as minutes after midnight.
#'
#'       - decimal_hour: Daily grid-cell centre expressed as decimal hours after midnight.
#'
#'       - reference_density: Estimated reference density at the grid cell.
#'
#'       - comparison_density: Estimated comparison density at the grid cell.
#'
#'       - shared_density: Pointwise minimum of the reference and comparison densities.
#'
#'     Groups with fewer than `min_n` observations are not represented in this
#'     table.
#'
#' @section Dependencies:
#' The function itself uses only functions from base R.
#'
#' Input time columns must inherit from class `"hms"`. The
#' [hms::as_hms()] function from the `hms` package can be used to construct
#' these columns.
#'
#' @references
#' Mardia, K. V., and Jupp, P. E. (2000).
#' *Directional Statistics*. Wiley.
#'
#' Ridout, M. S., and Linkie, M. (2009).
#' Estimating overlap of daily activity patterns from camera trap data.
#' *Journal of Agricultural, Biological, and Environmental Statistics*,
#' 14, 322--337.
#'
#' @seealso
#' [hms::as_hms()]
#'
#' @examples
#' if (requireNamespace("hms", quietly = TRUE)) {
#'
#'   set.seed(123)
#'
#'   # Reference observations
#'   reference_observations <- data.frame(
#'     observation_date = as.Date("2023-01-01") +
#'       sample(0:364, 80, replace = TRUE),
#'     observation_time = hms::as_hms(
#'       sample(0:(24 * 60 * 60 - 1), 80, replace = TRUE)
#'     )
#'   )
#'
#'   # Observations from three comparison groups
#'   comparison_observations <- data.frame(
#'     species = rep(
#'       c("species_a", "species_b", "species_c"),
#'       each = 60
#'     ),
#'     observation_date = as.Date("2023-01-01") +
#'       c(
#'         sample(20:180, 60, replace = TRUE),
#'         sample(120:300, 60, replace = TRUE),
#'         sample(250:364, 60, replace = TRUE)
#'       ),
#'     observation_time = hms::as_hms(
#'       c(
#'         sample(5:10, 60, replace = TRUE) * 3600,
#'         sample(10:17, 60, replace = TRUE) * 3600,
#'         sample(17:23, 60, replace = TRUE) * 3600
#'       )
#'     )
#'   )
#'
#'   # Compact output
#'   result <- torus_overlap_many(
#'     reference_data = reference_observations,
#'     comparison_data = comparison_observations,
#'     comparison_group_col = "species",
#'     reference_date_col = "observation_date",
#'     reference_time_col = "observation_time",
#'     reference_name = "reference_species",
#'     season_bw_days = 14,
#'     daily_bw_hours = 1,
#'     date_resolution_days = 5,
#'     time_resolution_minutes = 30
#'   )
#'
#'   result$overlap
#'   result$settings
#'
#'   # Optional long-format density output
#'   result_with_density <- torus_overlap_many(
#'     reference_data = reference_observations,
#'     comparison_data = comparison_observations,
#'     comparison_group_col = "species",
#'     reference_date_col = "observation_date",
#'     reference_time_col = "observation_time",
#'     reference_name = "reference_species",
#'     season_bw_days = 14,
#'     daily_bw_hours = 1,
#'     date_resolution_days = 5,
#'     time_resolution_minutes = 30,
#'     return_density = TRUE
#'   )
#'
#'   head(result_with_density$density)
#'
#'   # Different column names in the two input data frames
#'   names(reference_observations) <- c(
#'     "reference_date",
#'     "reference_time"
#'   )
#'
#'   result_different_names <- torus_overlap_many(
#'     reference_data = reference_observations,
#'     comparison_data = comparison_observations,
#'     comparison_group_col = "species",
#'     reference_date_col = "reference_date",
#'     reference_time_col = "reference_time",
#'     comparison_date_col = "observation_date",
#'     comparison_time_col = "observation_time",
#'     reference_name = "reference_species",
#'     date_resolution_days = 5,
#'     time_resolution_minutes = 30
#'   )
#'
#'   result_different_names$overlap
#' }
#'
#' @export
torus_overlap_many <- function(
    reference_data,
    comparison_data,
    comparison_group_col,
    reference_date_col,
    reference_time_col,
    comparison_date_col = reference_date_col,
    comparison_time_col = reference_time_col,
    reference_name = "reference",
    season_bw_days = 14,
    daily_bw_hours = 1,
    date_resolution_days = 2,
    time_resolution_minutes = 10,
    min_n = 5,
    return_density = FALSE
) {

  # ==========================================================================
  # Internal validation helpers
  # ==========================================================================

  check_positive_number <- function(x, argument) {

    if (
      !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0
    ) {
      stop(
        "`", argument,
        "` must be one positive finite number.",
        call. = FALSE
      )
    }
  }

  check_single_column_name <- function(x, argument) {

    if (
      !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
    ) {
      stop(
        "`", argument,
        "` must be one non-empty column name.",
        call. = FALSE
      )
    }
  }

  is_leap_year <- function(year) {

    year %% 400L == 0L |
      (
        year %% 4L == 0L &
          year %% 100L != 0L
      )
  }

  # ==========================================================================
  # Validate data-frame inputs
  # ==========================================================================

  if (!is.data.frame(reference_data)) {
    stop(
      "`reference_data` must be a data frame.",
      call. = FALSE
    )
  }

  if (!is.data.frame(comparison_data)) {
    stop(
      "`comparison_data` must be a data frame.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Validate column-name arguments
  # ==========================================================================

  check_single_column_name(
    comparison_group_col,
    "comparison_group_col"
  )

  check_single_column_name(
    reference_date_col,
    "reference_date_col"
  )

  check_single_column_name(
    reference_time_col,
    "reference_time_col"
  )

  check_single_column_name(
    comparison_date_col,
    "comparison_date_col"
  )

  check_single_column_name(
    comparison_time_col,
    "comparison_time_col"
  )

  missing_reference_columns <- setdiff(
    c(
      reference_date_col,
      reference_time_col
    ),
    names(reference_data)
  )

  if (length(missing_reference_columns) > 0L) {
    stop(
      "The following columns are missing from `reference_data`: ",
      paste(
        missing_reference_columns,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  missing_comparison_columns <- setdiff(
    c(
      comparison_group_col,
      comparison_date_col,
      comparison_time_col
    ),
    names(comparison_data)
  )

  if (length(missing_comparison_columns) > 0L) {
    stop(
      "The following columns are missing from `comparison_data`: ",
      paste(
        missing_comparison_columns,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Validate remaining arguments
  # ==========================================================================

  if (
    length(reference_name) != 1L ||
    is.na(reference_name)
  ) {
    stop(
      "`reference_name` must contain exactly one non-missing value.",
      call. = FALSE
    )
  }

  reference_name <- as.character(
    reference_name
  )

  check_positive_number(
    season_bw_days,
    "season_bw_days"
  )

  check_positive_number(
    daily_bw_hours,
    "daily_bw_hours"
  )

  check_positive_number(
    date_resolution_days,
    "date_resolution_days"
  )

  check_positive_number(
    time_resolution_minutes,
    "time_resolution_minutes"
  )

  if (
    !is.numeric(min_n) ||
    length(min_n) != 1L ||
    is.na(min_n) ||
    !is.finite(min_n) ||
    min_n < 2
  ) {
    stop(
      "`min_n` must be one finite number of at least 2.",
      call. = FALSE
    )
  }

  min_n <- as.integer(
    min_n
  )

  if (
    !is.logical(return_density) ||
    length(return_density) != 1L ||
    is.na(return_density)
  ) {
    stop(
      "`return_density` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Extract required columns
  # ==========================================================================

  reference_dat <- data.frame(
    date = reference_data[[reference_date_col]],
    time = reference_data[[reference_time_col]]
  )

  comparison_dat <- data.frame(
    group = comparison_data[[comparison_group_col]],
    date = comparison_data[[comparison_date_col]],
    time = comparison_data[[comparison_time_col]],
    stringsAsFactors = FALSE
  )

  comparison_dat$group <- as.character(
    comparison_dat$group
  )

  # ==========================================================================
  # Validate date and time classes
  # ==========================================================================

  if (!inherits(reference_dat$date, "Date")) {
    stop(
      "Column `", reference_date_col,
      "` in `reference_data` must have class 'Date'. ",
      "Convert it with `as.Date()` first.",
      call. = FALSE
    )
  }

  if (!inherits(reference_dat$time, "hms")) {
    stop(
      "Column `", reference_time_col,
      "` in `reference_data` must have class 'hms'. ",
      "Convert it with `hms::as_hms()` first.",
      call. = FALSE
    )
  }

  if (!inherits(comparison_dat$date, "Date")) {
    stop(
      "Column `", comparison_date_col,
      "` in `comparison_data` must have class 'Date'. ",
      "Convert it with `as.Date()` first.",
      call. = FALSE
    )
  }

  if (!inherits(comparison_dat$time, "hms")) {
    stop(
      "Column `", comparison_time_col,
      "` in `comparison_data` must have class 'hms'. ",
      "Convert it with `hms::as_hms()` first.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Convert hms values to seconds
  # ==========================================================================

  reference_dat$time_seconds <- as.numeric(
    reference_dat$time,
    units = "secs"
  )

  comparison_dat$time_seconds <- as.numeric(
    comparison_dat$time,
    units = "secs"
  )

  # ==========================================================================
  # Remove invalid reference observations
  # ==========================================================================

  valid_reference_rows <- (
    !is.na(reference_dat$date) &
      is.finite(reference_dat$time_seconds) &
      reference_dat$time_seconds >= 0 &
      reference_dat$time_seconds < 86400
  )

  n_reference_removed <- sum(
    !valid_reference_rows
  )

  reference_dat <- reference_dat[
    valid_reference_rows,
    ,
    drop = FALSE
  ]

  if (nrow(reference_dat) == 0L) {
    stop(
      "No valid observations remain in `reference_data`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Remove invalid comparison observations
  # ==========================================================================

  valid_comparison_rows <- (
    !is.na(comparison_dat$group) &
      nzchar(comparison_dat$group) &
      !is.na(comparison_dat$date) &
      is.finite(comparison_dat$time_seconds) &
      comparison_dat$time_seconds >= 0 &
      comparison_dat$time_seconds < 86400
  )

  n_comparison_removed <- sum(
    !valid_comparison_rows
  )

  comparison_dat <- comparison_dat[
    valid_comparison_rows,
    ,
    drop = FALSE
  ]

  if (nrow(comparison_dat) == 0L) {
    stop(
      "No valid observations remain in `comparison_data`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Convert dates and times to circular angles
  # ==========================================================================

  convert_to_angles <- function(
    date,
    time_seconds
  ) {

    year <- as.integer(
      format(date, "%Y")
    )

    year_length <- ifelse(
      is_leap_year(year),
      366L,
      365L
    )

    day_of_year <- as.integer(
      format(date, "%j")
    )

    # January 1 corresponds to zero radians.
    date_angle <- (
      2 * pi *
        (day_of_year - 1) /
        year_length
    )

    # Midnight corresponds to zero radians.
    time_angle <- (
      2 * pi *
        time_seconds /
        86400
    )

    list(
      date_angle = date_angle,
      time_angle = time_angle
    )
  }

  reference_angles <- convert_to_angles(
    date = reference_dat$date,
    time_seconds = reference_dat$time_seconds
  )

  reference_dat$date_angle <- (
    reference_angles$date_angle
  )

  reference_dat$time_angle <- (
    reference_angles$time_angle
  )

  comparison_angles <- convert_to_angles(
    date = comparison_dat$date,
    time_seconds = comparison_dat$time_seconds
  )

  comparison_dat$date_angle <- (
    comparison_angles$date_angle
  )

  comparison_dat$time_angle <- (
    comparison_angles$time_angle
  )

  # Columns no longer needed internally.
  reference_dat$date <- NULL
  reference_dat$time <- NULL
  reference_dat$time_seconds <- NULL

  comparison_dat$date <- NULL
  comparison_dat$time <- NULL
  comparison_dat$time_seconds <- NULL

  # ==========================================================================
  # Identify comparison groups
  # ==========================================================================

  comparison_groups <- sort(
    unique(comparison_dat$group)
  )

  if (length(comparison_groups) == 0L) {
    stop(
      "No comparison groups remain after filtering.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Construct density-evaluation grid
  # ==========================================================================

  # Used to define a common output grid and translate seasonal bandwidths
  # expressed in days. Individual observations were already transformed using
  # the actual length of their observation year.
  common_year_days <- 365.2425

  n_date <- max(
    2L,
    as.integer(
      ceiling(
        common_year_days /
          date_resolution_days
      )
    )
  )

  n_time <- max(
    2L,
    as.integer(
      ceiling(
        1440 /
          time_resolution_minutes
      )
    )
  )

  actual_date_resolution_days <- (
    common_year_days /
      n_date
  )

  actual_time_resolution_minutes <- (
    1440 /
      n_time
  )

  date_step_angle <- (
    2 * pi /
      n_date
  )

  time_step_angle <- (
    2 * pi /
      n_time
  )

  # Evaluate densities at grid-cell centres.
  date_grid <- (
    (seq_len(n_date) - 0.5) *
      date_step_angle
  )

  time_grid <- (
    (seq_len(n_time) - 0.5) *
      time_step_angle
  )

  cell_area <- (
    date_step_angle *
      time_step_angle
  )

  # Human-readable grid coordinates used for optional long output.
  date_grid_days <- (
    date_grid /
      (2 * pi) *
      common_year_days
  )

  time_grid_minutes <- (
    time_grid /
      (2 * pi) *
      1440
  )

  # ==========================================================================
  # Convert bandwidths to von Mises concentration parameters
  # ==========================================================================

  seasonal_sd_radians <- (
    2 * pi *
      season_bw_days /
      common_year_days
  )

  daily_sd_radians <- (
    2 * pi *
      daily_bw_hours /
      24
  )

  season_kappa <- (
    1 /
      seasonal_sd_radians^2
  )

  daily_kappa <- (
    1 /
      daily_sd_radians^2
  )

  # ==========================================================================
  # Numerically stable von Mises kernel
  # ==========================================================================

  von_mises_kernel <- function(
    evaluation_angle,
    observation_angle,
    kappa
  ) {

    angle_difference <- outer(
      evaluation_angle,
      observation_angle,
      "-"
    )

    scaled_bessel <- besselI(
      x = kappa,
      nu = 0,
      expon.scaled = TRUE
    )

    normalising_constant <- (
      2 * pi *
        scaled_bessel
    )

    exp(
      kappa *
        (
          cos(angle_difference) - 1
        )
    ) / normalising_constant
  }

  # ==========================================================================
  # Product-kernel toroidal density estimator
  # ==========================================================================

  estimate_density <- function(
    date_angle,
    time_angle
  ) {

    n_observations <- length(
      date_angle
    )

    if (
      n_observations !=
      length(time_angle)
    ) {
      stop(
        "Internal error: date and time vectors have unequal lengths.",
        call. = FALSE
      )
    }

    seasonal_kernel <- von_mises_kernel(
      evaluation_angle = date_grid,
      observation_angle = date_angle,
      kappa = season_kappa
    )

    daily_kernel <- von_mises_kernel(
      evaluation_angle = time_grid,
      observation_angle = time_angle,
      kappa = daily_kappa
    )

    density <- (
      seasonal_kernel %*%
        t(daily_kernel)
    ) / n_observations

    density_integral <- (
      sum(density) *
        cell_area
    )

    if (
      !is.finite(density_integral) ||
      density_integral <= 0
    ) {
      stop(
        "Density normalisation failed.",
        call. = FALSE
      )
    }

    # Correct small numerical-integration errors.
    density / density_integral
  }

  # ==========================================================================
  # Estimate reference density once
  # ==========================================================================

  n_reference <- nrow(
    reference_dat
  )

  if (n_reference < min_n) {
    stop(
      "The reference data contain only ",
      n_reference,
      " valid observations; at least ",
      min_n,
      " are required.",
      call. = FALSE
    )
  }

  reference_density <- estimate_density(
    date_angle = reference_dat$date_angle,
    time_angle = reference_dat$time_angle
  )

  # ==========================================================================
  # Prepare output containers
  # ==========================================================================

  overlap_results <- vector(
    mode = "list",
    length = length(comparison_groups)
  )

  if (return_density) {

    density_results <- vector(
      mode = "list",
      length = length(comparison_groups)
    )

    # R matrices are stored column-wise and converted with as.vector().
    long_date_angle <- rep(
      date_grid,
      times = n_time
    )

    long_time_angle <- rep(
      time_grid,
      each = n_date
    )

    long_season_day <- rep(
      date_grid_days,
      times = n_time
    )

    long_time_minutes <- rep(
      time_grid_minutes,
      each = n_date
    )

    long_reference_density <- as.vector(
      reference_density
    )
  }

  # ==========================================================================
  # Compare reference density with every comparison group
  # ==========================================================================

  for (i in seq_along(comparison_groups)) {

    comparison_name <- comparison_groups[i]

    comparison_group_data <- comparison_dat[
      comparison_dat$group == comparison_name,
      ,
      drop = FALSE
    ]

    n_comparison <- nrow(
      comparison_group_data
    )

    if (n_comparison < min_n) {

      overlap_results[[i]] <- data.frame(
        reference = reference_name,
        comparison = comparison_name,
        n_reference = n_reference,
        n_comparison = n_comparison,
        overlap = NA_real_,
        status = paste0(
          "Too few observations: n < ",
          min_n
        ),
        stringsAsFactors = FALSE
      )

      if (return_density) {
        density_results[[i]] <- NULL
      }

      next
    }

    comparison_density <- estimate_density(
      date_angle =
        comparison_group_data$date_angle,
      time_angle =
        comparison_group_data$time_angle
    )

    shared_density <- pmin(
      reference_density,
      comparison_density
    )

    overlap <- (
      sum(shared_density) *
        cell_area
    )

    # Protect against very small floating-point departures from [0, 1].
    overlap <- min(
      max(overlap, 0),
      1
    )

    overlap_results[[i]] <- data.frame(
      reference = reference_name,
      comparison = comparison_name,
      n_reference = n_reference,
      n_comparison = n_comparison,
      overlap = overlap,
      status = "OK",
      stringsAsFactors = FALSE
    )

    # ========================================================================
    # Optional long-format density data
    # ========================================================================

    if (return_density) {

      density_results[[i]] <- data.frame(
        reference = reference_name,
        comparison = comparison_name,
        date_angle = long_date_angle,
        time_angle = long_time_angle,
        seasonal_position_days =
          long_season_day,
        time_minutes =
          long_time_minutes,
        decimal_hour =
          long_time_minutes / 60,
        reference_density =
          long_reference_density,
        comparison_density =
          as.vector(comparison_density),
        shared_density =
          as.vector(shared_density),
        stringsAsFactors = FALSE
      )
    }

    # Do not retain matrices from the current comparison.
    rm(
      comparison_density,
      shared_density
    )
  }

  # ==========================================================================
  # Combine overlap results
  # ==========================================================================

  overlap_table <- do.call(
    rbind,
    overlap_results
  )

  rownames(overlap_table) <- NULL

  # ==========================================================================
  # Compact settings table
  # ==========================================================================

  settings <- data.frame(
    reference =
      reference_name,
    comparison_group_col =
      comparison_group_col,
    reference_date_col =
      reference_date_col,
    reference_time_col =
      reference_time_col,
    comparison_date_col =
      comparison_date_col,
    comparison_time_col =
      comparison_time_col,
    season_bw_days =
      season_bw_days,
    daily_bw_hours =
      daily_bw_hours,
    season_kappa =
      season_kappa,
    daily_kappa =
      daily_kappa,
    requested_date_resolution_days =
      date_resolution_days,
    actual_date_resolution_days =
      actual_date_resolution_days,
    requested_time_resolution_minutes =
      time_resolution_minutes,
    actual_time_resolution_minutes =
      actual_time_resolution_minutes,
    n_date =
      n_date,
    n_time =
      n_time,
    grid_cells =
      n_date * n_time,
    min_n =
      min_n,
    valid_reference_observations =
      n_reference,
    removed_reference_observations =
      n_reference_removed,
    valid_comparison_observations =
      nrow(comparison_dat),
    removed_comparison_observations =
      n_comparison_removed,
    comparison_groups =
      length(comparison_groups),
    stringsAsFactors = FALSE
  )

  # ==========================================================================
  # Construct returned object
  # ==========================================================================

  result <- list(
    overlap = overlap_table,
    settings = settings
  )

  if (return_density) {

    non_empty_density_results <- Filter(
      Negate(is.null),
      density_results
    )

    if (
      length(non_empty_density_results) == 0L
    ) {

      result$density <- data.frame(
        reference = character(0),
        comparison = character(0),
        date_angle = numeric(0),
        time_angle = numeric(0),
        seasonal_position_days = numeric(0),
        time_minutes = numeric(0),
        decimal_hour = numeric(0),
        reference_density = numeric(0),
        comparison_density = numeric(0),
        shared_density = numeric(0),
        stringsAsFactors = FALSE
      )

    } else {

      result$density <- do.call(
        rbind,
        non_empty_density_results
      )

      rownames(result$density) <- NULL
    }
  }

  class(result) <- c(
    "torus_overlap_many",
    class(result)
  )

  result
}
