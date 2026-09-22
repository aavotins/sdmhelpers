# Estimate Seasonal and Diel Activity Overlap on a Torus

Estimates the joint seasonal and diel (registration) activity density of
observations in a reference dataset and compares it with the
corresponding densities of one or more groups supplied in a separate
comparison dataset.

Both calendar date and time of day are treated as circular variables:
December is adjacent to January, and the end of one day is adjacent to
the beginning of the next. The resulting sample space is therefore a
two-dimensional torus.

A product-kernel density estimator based on two von Mises kernels is
used. Separate smoothing bandwidths are specified for the seasonal and
diel dimensions. The overlap coefficient between the reference density
and each comparison density is calculated by numerical integration of
the pointwise minimum of the two densities.

## Usage

``` r
torus_overlap_many(
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
)
```

## Arguments

- reference_data:

  A data frame containing the observations defining the reference
  distribution. It must contain a date column of class
  [Date](https://rdrr.io/r/base/Dates.html) and a time-of-day column
  inheriting from class `"hms"`.

  A grouping column is not required because all valid rows in
  `reference_data` are treated as observations from the same reference
  group.

- comparison_data:

  A data frame containing observations for one or more comparison
  groups. It must contain:

  - a grouping column;

  - a date column of class [Date](https://rdrr.io/r/base/Dates.html);
    and

  - a time-of-day column inheriting from class `"hms"`.

  The reference distribution is compared separately with every unique,
  non-missing, non-empty value in the grouping column.

- comparison_group_col:

  A character string giving the name of the grouping column in
  `comparison_data`.

- reference_date_col:

  A character string giving the name of the
  [Date](https://rdrr.io/r/base/Dates.html) column in `reference_data`.

- reference_time_col:

  A character string giving the name of the `"hms"` time-of-day column
  in `reference_data`.

- comparison_date_col:

  A character string giving the name of the
  [Date](https://rdrr.io/r/base/Dates.html) column in `comparison_data`.
  By default, the same column name as `reference_date_col` is used.

- comparison_time_col:

  A character string giving the name of the `"hms"` time-of-day column
  in `comparison_data`. By default, the same column name as
  `reference_time_col` is used.

- reference_name:

  A single value used to identify the reference dataset in the returned
  tables. The value is converted to character. The default is
  `"reference"`.

- season_bw_days:

  A positive numeric value giving the seasonal kernel bandwidth in days.
  The bandwidth is converted to an angular value in radians and then to
  the concentration parameter of the seasonal von Misses kernel. Larger
  values produce stronger smoothing across the annual cycle. The default
  is `14`.

- daily_bw_hours:

  A positive numeric value giving the daily kernel bandwidth in hours.
  The bandwidth is converted to an angular value in radians and then to
  the concentration parameter of the seasonal von Misses kernel. Larger
  values produce stronger smoothing across the daily cycle. The default
  is `1`.

- date_resolution_days:

  A positive numeric value giving the requested approximate spacing, in
  days, between evaluation points on the seasonal grid. This controls
  numerical integration resolution, not kernel smoothing. The realised
  spacing may differ slightly because the number of grid cells is
  rounded upward to an integer. The default is `2`.

- time_resolution_minutes:

  A positive numeric value giving the requested approximate spacing, in
  minutes, between evaluation points on the daily grid. This controls
  numerical integration resolution, not kernel smoothing. The realised
  spacing may differ slightly because the number of grid cells is
  rounded upward to an integer. The default is `10`.

- min_n:

  An integer giving the minimum number of valid observations required
  for density estimation. It is converted to an integer and must be at
  least `2`. The default is `5`.

  If the reference dataset contains fewer than `min_n` valid
  observations, the function stops with an error. Comparison groups with
  fewer than `min_n` observations are retained in the overlap table, but
  their overlap estimate is returned as `NA`.

- return_density:

  Logical. If `FALSE`, the default, only the compact overlap and
  settings tables are returned. If `TRUE`, a potentially large
  long-format data frame containing grid-level density values is also
  returned.

## Value

An object of class `"torus_overlap_many"` and `"list"` containing:

- overlap: A data frame with one row per comparison group and the
  following columns:

  - reference: Reference identifier supplied through `reference_name`.

  - comparison: Comparison-group identifier.

  - n_reference: Number of valid observations in the reference dataset.

  - n_comparison: Number of valid observations in the comparison group.

  - overlap: Estimated overlap coefficient. This is `NA` when the
    comparison group contains fewer than `min_n` observations.

  - status: `"OK"` for successfully estimated comparisons or a message
    indicating thet the comparison group contained fewer than `min_n`
    observations.

- settings: A one-row data frame describing column mappings, bandwidths,
  concentration parameters, requested and realised grid resolutions,
  grid dimensions, sample-size requirements, and numbers of valid and
  removed observations.

- density: Present only when `return_density = TRUE`. A long-format data
  frame containing:

  - reference: Reference identifier.

  - comparison: Comparison-group identifier.

  - date_angle: Seasonal grid-cell centre in radians.

  - time_angle: Diel grid-cell centre in radians.

  - seasonal_position_days: Seasonal grid-cell centre expressed as days
    from the beginning of a common 365.2425-day year used to construct
    the evaluation grid.

  - time_minutes: Diel grid-cell centre expressed as minutes after
    midnight.

  - decimal_hour: Diel grid-cell centre expressed as decimal hours after
    midnight.

  - reference_density: Normalized estimated reference density at the
    grid cell.

  - comparison_density: Normalized estimated comparison density at the
    grid cell.

  - shared_density: Pointwise minimum of the reference and comparison
    densities.

  Groups with fewer than `min_n` valid observations are not represented
  in this table.

## Details

### Circular transformation

Calendar dates are converted to seasonal angles. For observation \\i\\,
the seasonal angle is

\$\$ \theta_i = 2\pi \frac{d_i - 1}{L_i} \$\$

where \\d_i\\ is the day of year and \\L_i\\ is the length of the
observation year: either 365 or 366 days. January 1 therefore
corresponds to an angle of zero.

Consequently, leap-year observations are scaled according to the actual
length of their calendar year, so that both 365-day and 366-day years
are mapped onto the same \\2\pi\\-periodic seasonal cycle.

Time of day is converted to a diel angle using

\$\$ \phi_i = 2\pi \frac{s_i}{86400} \$\$

where \\s_i\\ is the number of seconds since midnight.

### Kernel density estimation

The joint density is estimated using a product of two von Mises kernels:

\$\$ \widehat{f}(\theta,\phi) = \frac{1}{n} \sum\_{i=1}^{n}
K\_{\kappa_s}(\theta-\theta_i) K\_{\kappa_t}(\phi-\phi_i) \$\$

where \\\kappa_s\\ and \\\kappa_t\\ are the seasonal and diel
concentration parameters.

User-supplied bandwidths are converted from days and hours to angular
bandwidths expressed in radians (approximate circular standard
deviations in radians). Concentration is then approximated as

\$\$ \kappa_s = \left(\frac{365.2425}{2\pi h_s}\right)^2 \$\$

\$\$ \kappa_t = \left(\frac{24}{2\pi h_t}\right)^2 \$\$

where \\h_s\\ is the seasonal bandwidth in days, \\h_t\\ is the diel
bandwidth in hours, and \\\kappa_s\\ and \\\kappa_t\\ are the
corresponding seasonal and diel von Mises concentration parameters.

Larger bandwidths produce smaller concentration parameters and therefore
broader von Mises kernels, resulting in stronger smoothing. Conversely,
smaller bandwidths produce larger concentration parameters and more
concentrated kernels.

The von Mises kernel is evaluated using an exponentially scaled modified
Bessel function. This provides a numerically stable form of the kernel
and reduces the risk of overflow for large concentration parameters.

### Density normalization

Each product-kernel estimate is theoretically a probability density on
the two-dimensional torus and therefore integrates to one. However, the
density is evaluated only at the centres of a finite regular grid, so
numerical integration of the evaluated density may differ slightly from
one.

The integral of the estimated density over the torus is approximated as

\$\$ I = \sum\_{j=1}^{n\_{\mathrm{date}}}
\sum\_{k=1}^{n\_{\mathrm{time}}} \widehat{f}(\theta_j,\phi_k)
\Delta\theta,\Delta\phi \$\$

where \\\Delta\theta = 2\pi/n\_{\mathrm{date}}\\ and \\\Delta\phi =
2\pi/n\_{\mathrm{time}}\\ are the angular dimensions of a grid cell.

To correct for numerical integration error, the evaluated density is
normalized as

\$\$ \widehat{f}{\mathrm{norm}}(\theta_j,\phi_k) =
\frac{\widehat{f}(\theta_j,\phi_k)}{I}. \$\$

Consequently,

\$\$ \sum{j=1}^{n\_{\mathrm{date}}} \sum\_{k=1}^{n\_{\mathrm{time}}}
\widehat{f}\_{\mathrm{norm}}(\theta_j,\phi_k) \Delta\theta,\Delta\phi
= 1. \$\$

The reference density and each comparison density are normalized
independently before their overlap is calculated. This normalization
does not alter the relative shape of either estimated density; it only
ensures that both numerically integrated densities represent probability
distributions on the evaluation grid.

### Overlap coefficient

After numerical normalization, the overlap coefficient between the
normalized reference density \\\widehat{f}\_r\\ and normalized
comparison density \\\widehat{f}\_c\\ is

\$\$ \Delta = \int_0^{2\pi} \int_0^{2\pi} \min\\
\widehat{f}\_r(\theta,\phi), \widehat{f}\_c(\theta,\phi) \\
\\d\theta\\d\phi. \$\$

For nomalized densities, it ranges from zero to one:

- `0` indicates no estimated overlap between the estimated dnsities;

- `1` indicates identical estimated densities.

The integral is approximated on a regular toroidal grid by evaluating
the densities at grid-cell centres and summing the pointwise minima of
the two densities multiplied by the grid-cell area. Decreasing
`date_resolution_days` or `time_resolution_minutes` increases numerical
resolution but also increases computing time and memory use.

### Missing and invalid observations

Rows are removed from `reference_data` when their date is missing or
when the time is missing, non-finite, earlier than 00:00, or equal to or
later than 24:00.

Rows are removed from `comparison_data` when:

- the comparison group is missing or an empty string;

- the date is missing; or

- the time is missing, non-finite, or outside the interval from midnight
  inclusive to 24:00 exclusive.

Counts of removed observations are included in the returned settings
table.

### Optional density output

When `return_density = TRUE`, density values are returned in a
long-format data frame rather than as matrices. Each comparison group
containing at least `\eqn{min_n}` observations contributes

\$\$ n\_{\mathrm{date}} \times n\_{\mathrm{time}} \$\$

rows.

The reference density is repeated for each comparison group with at
least `\eqn{min_n}` observations to make each group independently usable
for plotting and filtering. This output can therefore become large when
many groups or a fine evaluation grid are used.

Temporary kernel and density matrices are still created internally to
perform the calculations, but they are not returned.

## Dependencies

The function itself uses only functions from base R.

Input time columns must inherit from class `"hms"`. The
[`hms::as_hms()`](https://hms.tidyverse.org/reference/hms.html) function
from the `hms` package can be used to construct these columns.

## References

Mardia, K. V., and Jupp, P. E. (2000). *Directional Statistics*. Wiley.

Ridout, M. S., and Linkie, M. (2009). Estimating overlap of daily
activity patterns from camera trap data. *Journal of Agricultural,
Biological, and Environmental Statistics*, 14, 322–337.

## See also

[`hms::as_hms()`](https://hms.tidyverse.org/reference/hms.html)

## Examples

``` r
if (requireNamespace("hms", quietly = TRUE)) {

  set.seed(123)

  # Reference observations
  reference_observations <- data.frame(
    observation_date = as.Date("2023-01-01") +
      sample(0:364, 80, replace = TRUE),
    observation_time = hms::as_hms(
      sample(0:(24 * 60 * 60 - 1), 80, replace = TRUE)
    )
  )

  # Observations from three comparison groups
  comparison_observations <- data.frame(
    species = rep(
      c("species_a", "species_b", "species_c"),
      each = 60
    ),
    observation_date = as.Date("2023-01-01") +
      c(
        sample(20:180, 60, replace = TRUE),
        sample(120:300, 60, replace = TRUE),
        sample(250:364, 60, replace = TRUE)
      ),
    observation_time = hms::as_hms(
      c(
        sample(5:10, 60, replace = TRUE) * 3600,
        sample(10:17, 60, replace = TRUE) * 3600,
        sample(17:23, 60, replace = TRUE) * 3600
      )
    )
  )

  # Compact output
  result <- torus_overlap_many(
    reference_data = reference_observations,
    comparison_data = comparison_observations,
    comparison_group_col = "species",
    reference_date_col = "observation_date",
    reference_time_col = "observation_time",
    reference_name = "reference_species",
    season_bw_days = 14,
    daily_bw_hours = 1,
    date_resolution_days = 5,
    time_resolution_minutes = 30
  )

  result$overlap
  result$settings

  # Optional long-format density output
  result_with_density <- torus_overlap_many(
    reference_data = reference_observations,
    comparison_data = comparison_observations,
    comparison_group_col = "species",
    reference_date_col = "observation_date",
    reference_time_col = "observation_time",
    reference_name = "reference_species",
    season_bw_days = 14,
    daily_bw_hours = 1,
    date_resolution_days = 5,
    time_resolution_minutes = 30,
    return_density = TRUE
  )

  head(result_with_density$density)

  # Different column names in the two input data frames
  names(reference_observations) <- c(
    "reference_date",
    "reference_time"
  )

  result_different_names <- torus_overlap_many(
    reference_data = reference_observations,
    comparison_data = comparison_observations,
    comparison_group_col = "species",
    reference_date_col = "reference_date",
    reference_time_col = "reference_time",
    comparison_date_col = "observation_date",
    comparison_time_col = "observation_time",
    reference_name = "reference_species",
    date_resolution_days = 5,
    time_resolution_minutes = 30
  )

  result_different_names$overlap
}
#>           reference comparison n_reference n_comparison   overlap status
#> 1 reference_species  species_a          80           60 0.2023461     OK
#> 2 reference_species  species_b          80           60 0.2306367     OK
#> 3 reference_species  species_c          80           60 0.1731254     OK
```
