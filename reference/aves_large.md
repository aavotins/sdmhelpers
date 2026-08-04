# Large (species-wise) example bird-observation dataset

A reduced bird-observation dataset intended for examples that require
multiple observations, species, dates, or times. Contains 132647
observations from 200 species.

## Usage

``` r
aves_large
```

## Format

A `data.frame` containing 132647 example observations from 200 species
described by 7 variables:

- `species`: pseudonimised species.

- `coordX`: approximate X-coordinate of observation (EPSG: 3059).

- `coordY`: approximate Y-coordinate of observation (EPSG: 3059).

- `year`: the year of observation.

- `date`: the date (`as.Date` object) of the observation. Contains
  missing values.

- `date_time`: the `POSIXct` object of date-time. Contains missing
  values.

- `id100`: identificator of 100 m grid cell for joining with
  `data(points100)`.

## Source

A reduced and pseudonimised subset prepared for package examples.
