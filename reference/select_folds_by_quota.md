# Select Complete Fold Groups to Approximate a Target Sample Proportion

This function is intended for creating an independent testing subset
when observations have already been assigned to spatial, temporal,
environmental, or other grouped folds. It searches for a subset of
complete fold groups whose combined number of observations approximates
`target_prop`.

Candidate subsets are generated using repeated random orderings of the
fold groups. For each ordering, groups are greedily added while their
cumulative size does not exceed the upper bound. If the resulting subset
is smaller than the lower bound, one additional group is selected to
improve the solution.

Because folds are indivisible, an exact match to the target proportion
may not be possible. The function therefore returns the best solution
found during at most `max_tries` randomized attempts. In some cases,
particularly when one or more folds are very large, no combination may
satisfy `bounds`. The returned summary indicates whether the selected
solution lies within the requested bounds.

Missing fold identifiers are excluded by default. When
`include_na = TRUE`, all missing identifiers are treated as one
additional fold group.

## Usage

``` r
select_folds_by_quota(
  folds,
  target_prop = 0.25,
  bounds = c(0.2, 0.3),
  max_tries = 2000L,
  include_na = FALSE,
  seed = NULL,
  print_report = TRUE
)
```

## Arguments

- folds:

  A vector of fold identifiers, with one value per observation. Fold
  identifiers may be numeric, character, factor, logical, or another
  atomic vector type that can be processed by
  [`base::table()`](https://rdrr.io/r/base/table.html).

- target_prop:

  A single numeric value strictly between `0` and `1` giving the desired
  proportion of rows assigned to the selected folds. It must fall within
  the interval specified by `bounds`. The default is `0.25`.

- bounds:

  A numeric vector of length two giving the minimum and maximum
  acceptable proportions of rows assigned to the selected folds. Values
  must satisfy `0 <= bounds[1] <= target_prop <= bounds[2] <= 1`. The
  default is `c(0.20, 0.30)`.

- max_tries:

  A positive integer giving the maximum number of randomized fold
  orderings to evaluate. The search may stop earlier if it finds a
  solution containing exactly the target number of rows and lying within
  `bounds`. The default is `2000`.

- include_na:

  Logical. If `FALSE`, observations with missing fold identifiers are
  excluded from the fold-frequency table and are never selected. If
  `TRUE`, all missing fold identifiers are treated as a single fold
  group represented by `"<NA>"` in `selected_fold_ids`.

- seed:

  An optional single integer used to initialize the random-number
  generator before the search. Use a fixed value to make the selected
  folds reproducible. If `NULL`, the current random-number generator
  state is used.

- print_report:

  Logical. If `TRUE`, print a summary of the selected solution to the
  console. The returned object is invisible regardless of this setting.

## Value

Invisibly returns an object of class `"fold_quota_selection"`, which is
a list with the following components:

- `selected_fold_ids`: A character vector containing the identifiers of
  the selected fold groups. If missing values are included and selected,
  they are represented by `"<NA>"`.

- `selected_idx`: A logical vector with the same length as `folds`.
  Values are `TRUE` for rows belonging to selected fold groups and
  `FALSE` otherwise.

- `summary`: A list containing:

- `n_rows_total`: Number of rows included in the selection problem.
  Missing fold identifiers are omitted when `include_na = FALSE`.

- `n_rows_selected`: Number of rows belonging to the selected folds.

- `prop_selected`: Proportion of included rows belonging to the selected
  folds.

- `bounds`: Requested lower and upper proportional bounds.

- `row_bounds`: Integer lower and upper bounds used during selection.

- `target_prop`: Requested target proportion.

- `target_n`: Target number of selected rows after rounding.

- `within_bounds`: Logical value indicating whether the returned
  solution lies within the requested integer bounds.

- `exact_target`: Logical value indicating whether the selected row
  count equals `target_n`.

- `n_groups`: Total number of fold groups considered.

- `n_groups_selected`: Number of fold groups selected.

- `iterations_attempted`: Number of randomized attempts performed.

- `best_iteration`: Attempt on which the best returned solution was
  first found.

## Details

Selects complete fold groups so that the number of rows belonging to the
selected folds is as close as possible to a requested proportion of all
rows. Fold groups are never split between the selected and unselected
sets.

The row-count constraints are calculated as:


    lower   = ceiling(bounds[1] * N)
    upper   = floor(bounds[2] * N)
    target  = round(target_prop * N)

where `N` is the number of observations included in the selection
problem.

The procedure is a randomized heuristic rather than an exhaustive
combinatorial search. Consequently, it does not guarantee the globally
optimal combination of folds. Increasing `max_tries` can improve the
probability of finding a good combination when many fold groups are
available.

The same fold identifiers and seed produce the same result under the
same R version and random-number generator settings. Calling this
function with a non-`NULL` `seed` changes the global random-number
generator state.

Fold identifiers in `selected_fold_ids` are returned as character values
because [`base::table()`](https://rdrr.io/r/base/table.html) uses group
labels as names. Use `selected_idx` when subsetting the original data,
because it preserves the direct relationship with the input rows and
avoids type-conversion issues.

## Missing fold identifiers

When `include_na = FALSE`, rows with missing fold identifiers:

- are not included in `n_rows_total`;

- do not affect the target or bounds;

- receive `FALSE` in `selected_idx`.

When `include_na = TRUE`, all missing fold identifiers form one complete
group. Therefore, either all or none of the missing-fold rows are
selected.

## Interpretation

This function selects complete fold groups based only on their sizes. It
does not assess whether the selected observations are environmentally,
geographically, temporally, or statistically independent. Independence
must arise from the procedure used to define `folds`.

## See also

[`base::table()`](https://rdrr.io/r/base/table.html),
[`base::sample.int()`](https://rdrr.io/r/base/sample.html),
[`base::set.seed()`](https://rdrr.io/r/base/Random.html)

## Examples

``` r
## Example with equally sized folds
folds <- rep(letters[1:10], each = 10)

selection <- select_folds_by_quota(
  folds,
  target_prop = 0.25,
  bounds = c(0.20, 0.30),
  seed = 123,
  print_report = FALSE
)

selection$selected_fold_ids
#> [1] "b" "c" "j"
selection$summary
#> $n_rows_total
#> [1] 100
#> 
#> $n_rows_selected
#> [1] 30
#> 
#> $prop_selected
#> [1] 0.3
#> 
#> $bounds
#> [1] 0.2 0.3
#> 
#> $row_bounds
#> lower upper 
#>    20    30 
#> 
#> $target_prop
#> [1] 0.25
#> 
#> $target_n
#> [1] 25
#> 
#> $within_bounds
#> [1] TRUE
#> 
#> $exact_target
#> [1] FALSE
#> 
#> $n_groups
#> [1] 10
#> 
#> $n_groups_selected
#> [1] 3
#> 
#> $iterations_attempted
#> [1] 2000
#> 
#> $best_iteration
#> [1] 1
#> 
table(selection$selected_idx)
#> 
#> FALSE  TRUE 
#>    70    30 

## Use selected_idx to divide a data frame
dat <- data.frame(
  occurrence_id = seq_along(folds),
  fold = folds
)

independent_test <- dat[selection$selected_idx, , drop = FALSE]
model_training <- dat[!selection$selected_idx, , drop = FALSE]

## Verify that complete folds were selected
unique(independent_test$fold)
#> [1] "b" "c" "j"
intersect(
  unique(independent_test$fold),
  unique(model_training$fold)
)
#> character(0)

## Unequal fold sizes
unequal_folds <- rep(
  paste0("fold_", 1:8),
  times = c(5, 8, 12, 15, 18, 20, 7, 15)
)

unequal_selection <- select_folds_by_quota(
  unequal_folds,
  target_prop = 0.25,
  bounds = c(0.20, 0.30),
  max_tries = 5000,
  seed = 42,
  print_report = FALSE
)

unequal_selection$summary$prop_selected
#> [1] 0.25
unequal_selection$summary$within_bounds
#> [1] TRUE

## Missing fold identifiers excluded
folds_with_na <- c("A", "A", "B", "B", "C", NA, NA)

without_na <- select_folds_by_quota(
  folds_with_na,
  target_prop = 0.4,
  bounds = c(0.3, 0.5),
  include_na = FALSE,
  seed = 1,
  print_report = FALSE
)

without_na$selected_idx
#> [1]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE

## Treat all missing identifiers as one fold
with_na <- select_folds_by_quota(
  folds_with_na,
  target_prop = 0.4,
  bounds = c(0.3, 0.5),
  include_na = TRUE,
  seed = 1,
  print_report = FALSE
)

with_na$selected_fold_ids
#> [1] "A" "C"
with_na$selected_idx
#> [1]  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE
```
