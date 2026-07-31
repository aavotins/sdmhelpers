# Calculate TSS Statistics During ENMeval Evaluation

A user-defined evaluation function for
[`ENMeval::ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.html)
that calculates the maximum True Skill Statistic for the training and
validation data of each model partition.

## Usage

``` r
tss_user_eval(vars)
```

## Arguments

- vars:

  A named list supplied internally by
  [`ENMeval::ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.html)
  to its `user.eval` function. The function uses the following
  components:

  - `occs.train.pred`: Numeric predictions for training occurrences.

  - `occs.val.pred`: Numeric predictions for validation occurrences.

  - `bg.train.pred`: Numeric predictions for training background
    records.

  - `bg.val.pred`: Numeric predictions for validation background
    records.

  - `other.settings`: A named list containing the settings used by
    [`ENMeval::ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.html),
    including `validation.bg`.

## Value

A one-row data frame containing:

- `tss.val`: Maximum validation TSS.

- `thr.val`: Threshold maximizing validation TSS.

- `sens.val`: Validation sensitivity at `thr.val`.

- `spec.val`: Validation specificity at `thr.val`.

- `tss.train`: Maximum training TSS.

If either the occurrence or background predictions required for a
calculation contain no finite values, the corresponding results are
returned as `NA_real_`.

## Details

This function is intended to be supplied to the `user.eval` argument of
[`ENMeval::ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.html):


    ENMeval::ENMevaluate(
      ...,
      user.eval = tss_user_eval
    )

For each cross-validation fold, the function calculates:

- validation TSS, threshold, sensitivity, and specificity; and

- training TSS.

The validation background follows the setting used by ENMeval:

- when `other.settings$validation.bg = "full"`, validation occurrences
  are compared with the combined training and validation background
  predictions;

- otherwise, validation occurrences are compared only with
  validation-background predictions.

The training statistic always compares training occurrences with
training background predictions.

The maximum TSS threshold is selected separately for the training and
validation data. Consequently, `tss.train` and `tss.val` are optimistic
estimates of performance within their respective datasets because each
is maximized using the same data on which it is evaluated.

For an independent assessment of a threshold-based model, a threshold
should generally be selected using training data and then applied,
without re-optimization, to independent validation or testing data.

ENMeval summarizes fold-level statistics and adds columns such as `.avg`
and `.sd` to its model-level results. The exact naming and placement of
these summaries are controlled by ENMeval.

## References

Allouche, O., Tsoar, A., and Kadmon, R. (2006). Assessing the accuracy
of species distribution models: prevalence, kappa and the true skill
statistic (TSS). *Journal of Applied Ecology*, **43**, 1223–1232.
[doi:10.1111/j.1365-2664.2006.01214.x](https://doi.org/10.1111/j.1365-2664.2006.01214.x)

Kass, J. M., Muscarella, R., Galante, P. J., Bohl, C. L.,
Pinilla-Buitrago, G. E., Boria, R. A., Soley-Guardia, M., and Anderson,
R. P. (2021). ENMeval 2.0: Redesigned for customizable and reproducible
modeling of species' niches and distributions. *Methods in Ecology and
Evolution*, **12**, 1602–1608.
[doi:10.1111/2041-210X.13628](https://doi.org/10.1111/2041-210X.13628)

## See also

[`max_tss()`](https://aavotins.github.io/sdmhelpers/reference/max_tss.md),
[`ENMeval::ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.html)

## Examples

``` r
# A minimal object resembling the list supplied by ENMeval
example_vars <- list(
  occs.train.pred = c(0.95, 0.85, 0.70, 0.55),
  occs.val.pred   = c(0.90, 0.65, 0.45),
  bg.train.pred   = c(0.60, 0.40, 0.25, 0.10),
  bg.val.pred     = c(0.55, 0.35, 0.15),
  other.settings  = list(validation.bg = "partition")
)

tss_user_eval(example_vars)
#>     tss.val thr.val  sens.val spec.val tss.train
#> 1 0.6666667    0.65 0.6666667        1      0.75

# Use the complete background for validation
example_vars$other.settings$validation.bg <- "full"
tss_user_eval(example_vars)
#>     tss.val thr.val sens.val  spec.val tss.train
#> 1 0.7142857    0.45        1 0.7142857      0.75

if (FALSE) { # \dontrun{
# Example integration with ENMeval
evaluation <- ENMeval::ENMevaluate(
  occs = occurrence_data,
  envs = environmental_rasters,
  bg = background_data,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = list(
    fc = c("L", "LQ", "LQH"),
    rm = c(0.5, 1, 2)
  ),
  partition.settings = list(
    orientation = "lat_lon"
  ),
  other.settings = list(
    validation.bg = "partition"
  ),
  user.eval = tss_user_eval
)

ENMeval::eval.results(evaluation)
ENMeval::eval.results.partitions(evaluation)
} # }
```
