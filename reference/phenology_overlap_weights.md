# Calculate phenological overlap weights among species

Calculates the overlap between the day-of-year distribution of a target
species and the day-of-year distributions of multiple comparison
species. Kernel-density overlap is calculated separately for each
comparison species using
[`overlapping::overlap()`](https://rdrr.io/pkg/overlapping/man/overlap.html).

## Usage

``` r
phenology_overlap_weights(
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
)
```

## Arguments

- target:

  A data frame containing observations of one or more potential target
  species.

- group:

  A data frame containing observations of the comparison species.

- target_code:

  Species code identifying the target species in `target`.

- target_id:

  Character string giving the species-code column in `target`. Defaults
  to `"CODE"`.

- group_id:

  Character string giving the species-code column in `group`. Defaults
  to `"CODE"`.

- target_doy:

  Character string giving the day-of-year column in `target`. Defaults
  to `"DoY"`.

- group_doy:

  Character string giving the day-of-year column in `group`. Defaults to
  `"DoY"`.

- type:

  Character string specifying the overlap measure passed to
  [`overlapping::overlap()`](https://rdrr.io/pkg/overlapping/man/overlap.html).
  Either `"1"` or `"2"`. Defaults to `"2"`.

- nbins:

  Integer giving the number of equally spaced points used for comparing
  the estimated densities. Passed to
  [`overlapping::overlap()`](https://rdrr.io/pkg/overlapping/man/overlap.html).
  Defaults to `1024`.

- boundaries:

  Optional numeric vector of length two defining the lower and upper
  boundaries over which overlap is evaluated. Defaults to `NULL`.

- min_n:

  Minimum number of finite day-of-year observations required for the
  target species and each comparison species. Defaults to `2`.

- na_rm:

  Logical. If `TRUE`, non-finite day-of-year values are removed before
  estimating densities. Defaults to `TRUE`.

- ...:

  Additional arguments passed to
  [`overlapping::overlap()`](https://rdrr.io/pkg/overlapping/man/overlap.html)
  and subsequently to
  [`stats::density()`](https://rdrr.io/r/stats/density.html), for
  example `bw`, `adjust`, or `kernel`.

## Value

A data frame containing the comparison-species code, the estimated
overlap with the target species, and the normalized overlap weight.

## Details

The resulting overlap estimates are additionally normalized across all
comparison species so that their weights sum to 1.

The target species is selected from `target` using `target_code` and the
column specified by `target_id`.

Each species in `group` is then compared independently with the target
species. The normalized weight for species \\i\\ is

\$\$w_i = O_i / \sum_j O_j\$\$

where \\O_i\\ is the estimated overlap between the target species and
comparison species \\i\\.

Missing overlap estimates are excluded from the denominator. If all
valid overlap estimates are zero, weights are returned as `NA`.

Note that
[`overlapping::overlap()`](https://rdrr.io/pkg/overlapping/man/overlap.html)
uses ordinary linear kernel-density estimation. Day of year is circular,
so observations near the beginning and end of the year are not treated
as adjacent.

## References

Pastore, M. (2018). Overlapping: an R package for estimating overlapping
in empirical distributions. Journal of Open Source Software, 3(32),
1023. [doi:10.21105/joss.01023](https://doi.org/10.21105/joss.01023)

Pastore, M., & Calcagni, A. (2019). Measuring distribution similarities
between samples: A distribution-free overlapping index. Frontiers in
Psychology, 10, 1089.
[doi:10.3389/fpsyg.2019.01089](https://doi.org/10.3389/fpsyg.2019.01089)

## See also

[`overlapping::overlap()`](https://rdrr.io/pkg/overlapping/man/overlap.html),
[`stats::density()`](https://rdrr.io/r/stats/density.html)

## Examples

``` r
set.seed(123)

target <- data.frame(
  CODE = rep(c("TARGET1", "TARGET2"), each = 100),
  DoY = c(
    round(rnorm(100, 140, 20)),
    round(rnorm(100, 220, 20))
  )
)

group <- data.frame(
  CODE = rep(c("SP1", "SP2", "SP3"), each = 100),
  DoY = c(
    round(rnorm(100, 145, 20)),
    round(rnorm(100, 190, 20)),
    round(rnorm(100, 250, 20))
  )
)

phenology_overlap_weights(
  target = target,
  group = group,
  target_code = "TARGET1"
)
#>     CODE     overlap      weight
#> SP1  SP1 0.830042695 0.842134033
#> SP2  SP2 0.152366190 0.154585727
#> SP3  SP3 0.003233142 0.003280239
```
