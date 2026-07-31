# Select joint presence and background folds for independent testing

Selects a common set of complete fold IDs for an independent testing
dataset. The selected folds are chosen so that the proportion of
selected presence records is close to a requested target, while the
proportions of both presence and background records remain within
user-defined bounds.

## Usage

``` r
select_joint_folds(
  pres_folds,
  bg_folds,
  target_prop = 0.25,
  pres_bounds = c(0.2, 0.3),
  bg_bounds = c(0.2, 0.3),
  max_tries = 3000,
  bg_match_weight = 1,
  include_na = FALSE,
  seed = NULL,
  print_report = TRUE
)
```

## Arguments

- pres_folds:

  A vector containing fold identifiers for presence records. Fold
  identifiers may be numeric, character, factor, or another atomic type
  coercible to character.

- bg_folds:

  A vector containing fold identifiers for background records. Fold
  identifiers may be numeric, character, factor, or another atomic type
  coercible to character.

- target_prop:

  A single numeric value strictly between `0` and `1`. The desired
  proportion of in-scope presence records assigned to the independent
  testing dataset. The default is `0.25`.

- pres_bounds:

  A numeric vector of length two giving the minimum and maximum
  permitted proportions of presence records in the selected folds.
  Values must lie between `0` and `1`, and the first value must not
  exceed the second. The default is `c(0.20, 0.30)`.

- bg_bounds:

  A numeric vector of length two giving the minimum and maximum
  permitted proportions of background records in the selected folds.
  Values must lie between `0` and `1`, and the first value must not
  exceed the second. The default is `c(0.20, 0.30)`.

- max_tries:

  A positive whole number giving the maximum number of randomized fold
  orderings evaluated by the heuristic search. Larger values may improve
  the selected solution when fold sizes are irregular, at the cost of
  additional computation. The default is `3000`.

- bg_match_weight:

  A non-negative numeric value controlling the importance of matching
  the selected background proportion to the selected presence
  proportion. A value of `0` disables proportional matching, although
  `bg_bounds` are still considered. Larger values place progressively
  greater emphasis on making the two selected proportions similar. The
  default is `1`.

- include_na:

  Logical. If `FALSE`, records with missing fold IDs are excluded from
  fold-size calculations and cannot be selected. Their corresponding
  values in the returned row-selection vectors are `FALSE`. If `TRUE`,
  missing fold IDs are treated as one additional fold represented
  internally and in the returned fold IDs by `"<NA>"`. The default is
  `FALSE`.

- seed:

  `NULL` or a single integer-like value passed to
  [`base::set.seed()`](https://rdrr.io/r/base/Random.html). Supplying a
  seed makes the randomized search reproducible. When `NULL`, the
  current random-number generator state is used.

- print_report:

  Logical. If `TRUE`, print a summary of the selected folds, sample
  sizes, proportions, iteration number, and final score. The default is
  `TRUE`.

## Value

Invisibly returns a named list with the following elements:

- `selected_fold_ids`: A character vector containing the selected fold
  IDs.

- `pres_selected_idx`: A logical vector of length `length(pres_folds)`.
  Values are `TRUE` for presence records belonging to selected folds.

- `bg_selected_idx`: A logical vector of length `length(bg_folds)`.
  Values are `TRUE` for background records belonging to selected folds.

- `summary`: A list containing:

  - `n_pres_total`: Total number of elements in `pres_folds`, including
    missing values.

  - `n_bg_total`: Total number of elements in `bg_folds`, including
    missing values.

  - `n_pres_in_scope`: Number of presence records considered during
    selection after applying the requested missing-value handling.

  - `n_bg_in_scope`: Number of background records considered during
    selection after applying the requested missing-value handling.

  - `n_pres_selected`: Number of selected in-scope presence records.

  - `n_bg_selected`: Number of selected in-scope background records.

  - `prop_pres_selected`: Proportion of in-scope presence records
    selected.

  - `prop_bg_selected`: Proportion of in-scope background records
    selected.

  - `target_prop`: Requested presence selection proportion.

  - `bounds_pres`: Requested presence proportion bounds.

  - `bounds_bg`: Requested background proportion bounds.

  - `target_pres_abs`: Rounded target number of presence records.

  - `target_bg_abs_range`: Integer lower and upper bounds for the
    selected background count.

  - `iterations`: Iteration at which the best solution was first
    encountered.

  - `score`: Score of the selected solution. Smaller values indicate
    better agreement with the requested target, bounds, and proportional
    matching criterion.

## Details

The same fold IDs are selected for both datasets. Consequently, all
presence and background records belonging to a selected fold are
assigned to the independent testing dataset.

This function is intended for spatial or otherwise grouped model
evaluation in which complete folds must be withheld from both presence
and background datasets. It does not select individual records
independently.

Fold sizes are first counted separately for the presence and background
data. The union of the presence and background fold IDs defines the
candidate folds. A fold that occurs in only one dataset therefore has
size zero in the other dataset.

The function performs up to `max_tries` randomized greedy searches. For
each randomized ordering, candidate folds are added when doing so
improves the objective score or when additional presence records are
needed to reach the lower presence bound.

The objective score combines:

1.  deviation of the selected presence count from the requested presence
    target;

2.  penalties for selected presence counts outside `pres_bounds`;

3.  penalties for selected background counts outside `bg_bounds`; and

4.  the absolute difference between selected presence and background
    proportions, weighted by `bg_match_weight`.

Because complete folds are indivisible, it may be impossible to satisfy
all requested bounds or to attain `target_prop` exactly. The returned
result is the best solution encountered by the heuristic search; it is
not guaranteed to be the global optimum.

Search quality may be improved by increasing `max_tries`, especially
when:

- the number of folds is large;

- fold sizes vary substantially;

- presence and background fold-size distributions differ; or

- the requested bounds are narrow.

## Note

The function optimizes the number of records assigned to the independent
dataset, not the number of selected folds. A small number of large folds
can therefore be preferred over a larger number of small folds.

The returned logical vectors can be used directly to separate testing
and training records:

`pres_test <- pres_data[result$pres_selected_idx, ]`

`pres_train <- pres_data[!result$pres_selected_idx, ]`

`bg_test <- bg_data[result$bg_selected_idx, ]`

`bg_train <- bg_data[!result$bg_selected_idx, ]`

## Interpretation of missing fold IDs

With `include_na = FALSE`, missing fold IDs are outside the selection
problem. They are included in `n_pres_total` or `n_bg_total`, but
excluded from `n_pres_in_scope` or `n_bg_in_scope`.

With `include_na = TRUE`, all missing fold IDs are treated as belonging
to one common fold named `"<NA>"`. Therefore, selecting that fold
selects every record with a missing fold ID in the corresponding
dataset.

To avoid ambiguity when `include_na = TRUE`, the input fold identifiers
should not contain an actual, non-missing fold ID equal to `"<NA>"`.

## Random-number generation

If `seed` is supplied, the function calls
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) and therefore
changes R's global random-number generator state. The same inputs, seed,
and R version should produce the same result.

## See also

[`base::table()`](https://rdrr.io/r/base/table.html),
[`base::set.seed()`](https://rdrr.io/r/base/Random.html),
[`base::sample.int()`](https://rdrr.io/r/base/sample.html)

## Examples

``` r
## Create fold assignments with unequal fold sizes.
pres_folds <- rep(
  paste0("fold_", 1:8),
  times = c(12, 15, 9, 18, 14, 11, 16, 10)
)

bg_folds <- rep(
  paste0("fold_", 1:8),
  times = c(80, 95, 70, 120, 90, 75, 110, 85)
)

result <- select_joint_folds(
  pres_folds = pres_folds,
  bg_folds = bg_folds,
  target_prop = 0.25,
  pres_bounds = c(0.20, 0.30),
  bg_bounds = c(0.20, 0.30),
  max_tries = 500,
  seed = 123,
  print_report = FALSE
)

result$selected_fold_ids
#> [1] "fold_6" "fold_7"
result$summary
#> $n_pres_total
#> [1] 105
#> 
#> $n_bg_total
#> [1] 725
#> 
#> $n_pres_in_scope
#> [1] 105
#> 
#> $n_bg_in_scope
#> [1] 725
#> 
#> $n_pres_selected
#> fold_6 
#>     27 
#> 
#> $n_bg_selected
#> fold_6 
#>    185 
#> 
#> $prop_pres_selected
#>    fold_6 
#> 0.2571429 
#> 
#> $prop_bg_selected
#>    fold_6 
#> 0.2551724 
#> 
#> $target_prop
#> [1] 0.25
#> 
#> $bounds_pres
#> [1] 0.2 0.3
#> 
#> $bounds_bg
#> [1] 0.2 0.3
#> 
#> $target_pres_abs
#> [1] 26
#> 
#> $target_bg_abs_range
#> [1] 145 217
#> 
#> $iterations
#> [1] 18
#> 
#> $score
#>   fold_6 
#> 2.428571 
#> 

## Extract selected records.
pres_test_folds <- pres_folds[result$pres_selected_idx]
bg_test_folds <- bg_folds[result$bg_selected_idx]

unique(pres_test_folds)
#> [1] "fold_6" "fold_7"
unique(bg_test_folds)
#> [1] "fold_6" "fold_7"

## Verify that both datasets use the same selected fold IDs.
setequal(
  unique(pres_test_folds),
  unique(bg_test_folds)
)
#> [1] TRUE

## Presence and background tables do not need to contain exactly the
## same folds. A missing fold in one dataset receives a size of zero.
pres_folds2 <- rep(letters[1:6], c(7, 13, 10, 15, 9, 11))
bg_folds2 <- rep(letters[c(1:5, 7)], c(60, 90, 75, 100, 70, 50))

result2 <- select_joint_folds(
  pres_folds = pres_folds2,
  bg_folds = bg_folds2,
  target_prop = 0.25,
  pres_bounds = c(0.15, 0.35),
  bg_bounds = c(0.15, 0.35),
  max_tries = 500,
  seed = 42,
  print_report = FALSE
)

result2$selected_fold_ids
#> [1] "d"

## Missing fold IDs can be excluded.
pres_with_na <- c(pres_folds, NA, NA)
bg_with_na <- c(bg_folds, NA)

result3 <- select_joint_folds(
  pres_folds = pres_with_na,
  bg_folds = bg_with_na,
  include_na = FALSE,
  max_tries = 500,
  seed = 7,
  print_report = FALSE
)

result3$pres_selected_idx[is.na(pres_with_na)]
#> [1] FALSE FALSE

## Alternatively, all missing IDs can be treated as one selectable fold.
result4 <- select_joint_folds(
  pres_folds = pres_with_na,
  bg_folds = bg_with_na,
  include_na = TRUE,
  max_tries = 500,
  seed = 7,
  print_report = FALSE
)

result4$selected_fold_ids
#> [1] "fold_1" "fold_7"
```
