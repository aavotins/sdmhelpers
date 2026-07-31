# Calculate the Maximum True Skill Statistic

Calculates the maximum True Skill Statistic (TSS) across all unique
prediction thresholds for presence and background or absence
predictions.

## Usage

``` r
max_tss(p_pos, p_neg)
```

## Arguments

- p_pos:

  A numeric vector containing predicted suitability values for observed
  presences.

- p_neg:

  A numeric vector containing predicted suitability values for
  background records or observed absences.

## Value

A named numeric vector of length four containing:

- `tss`: Maximum True Skill Statistic.

- `threshold`: Threshold producing the maximum TSS.

- `sensitivity`: Sensitivity at the selected threshold.

- `specificity`: Specificity at the selected threshold.

## Details

The True Skill Statistic is calculated as

\$\$\mathrm{TSS} = \mathrm{sensitivity} + \mathrm{specificity} - 1.\$\$

At each threshold, predictions greater than or equal to the threshold
are classified as positive, whereas predictions below the threshold are
classified as negative.

Thresholds are evaluated at every unique finite prediction value. Two
additional classifications are considered:

- all records classified as negative; and

- all records classified as positive.

Both extreme classifications have a TSS of zero.

Non-finite values, including `NA`, `NaN`, `Inf`, and `-Inf`, are removed
separately from `p_pos` and `p_neg`. If no finite values remain in
either vector, all returned values are `NA_real_`.

When multiple thresholds produce the same maximum TSS, the highest
threshold is returned. This corresponds to the most restrictive
classification among the tied thresholds. One exception occurs when no
evaluated threshold has a TSS greater than zero: the function returns
`Inf`, representing classification of all records as negative.

TSS is independent of prevalence when sensitivity and specificity are
calculated from representative presence and negative samples. However,
results can still depend strongly on the spatial, environmental, and
sampling properties of the background or absence data.

## References

Allouche, O., Tsoar, A., and Kadmon, R. (2006). Assessing the accuracy
of species distribution models: prevalence, kappa and the true skill
statistic (TSS). *Journal of Applied Ecology*, **43**, 1223–1232.
[doi:10.1111/j.1365-2664.2006.01214.x](https://doi.org/10.1111/j.1365-2664.2006.01214.x)

## See also

[`tss_user_eval()`](https://aavotins.github.io/sdmhelpers/reference/tss_user_eval.md)

## Examples

``` r
presence_predictions <- c(0.90, 0.85, 0.75, 0.60, 0.40)
background_predictions <- c(0.70, 0.50, 0.30, 0.20, 0.10)

max_tss(
  p_pos = presence_predictions,
  p_neg = background_predictions
)
#>         tss   threshold sensitivity specificity 
#>        0.60        0.75        0.60        1.00 

# Non-finite values are removed
max_tss(
  p_pos = c(0.9, 0.8, NA, Inf),
  p_neg = c(0.4, 0.2, NA)
)
#>         tss   threshold sensitivity specificity 
#>         1.0         0.8         1.0         1.0 

# Apply the selected threshold manually
result <- max_tss(
  p_pos = presence_predictions,
  p_neg = background_predictions
)

threshold <- unname(result["threshold"])
presence_classes <- presence_predictions >= threshold
background_classes <- background_predictions >= threshold
```
