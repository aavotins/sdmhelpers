# Small (species-wise) example bird-observation dataset

A small bird-observation dataset intended for fast examples. Contains
6539 observations from 5 species.

## Usage

``` r
aves_small
```

## Format

A `data.frame` containing 6359 example observations from 5 species
described by 6 variables:

- `species`: scientific name of observed species.

- `code`: abbreviation of the scientific name.

- `year`: the year of observation.

- `date`: the date (`as.Date` object) of the observation.

- `date_time`: the `POSIXct` object of date-time. Contains missing
  values.

- `id100`: identificator of 100 m grid cell for joining with
  `data(points100)`.

## Source

A reduced subset prepared for package examples.
