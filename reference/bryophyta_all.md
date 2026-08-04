# Example bryophyte observations

A reduced dataset containing 12587 observations of 551 bryophyte taxa.

## Usage

``` r
bryophyta_all
```

## Format

A `data.frame` containing 12587 example observations from 551 species
described by 6 variables:

- `species`: pseudonimised species.

- `coordX`: approximate X-coordinate of observation (EPSG: 3059).

- `coordY`: approximate Y-coordinate of observation (EPSG: 3059).

- `year`: the year of observation.

- `date`: the date (`as.Date` object) of the observation. Contains
  missing values.

- `id100`: identificator of 100 m grid cell for joining with
  `data(points100)`.

## Source

A reduced and pseudonimised subset prepared for package examples.
