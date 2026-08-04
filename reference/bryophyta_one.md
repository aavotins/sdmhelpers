# Example observations for one bryophyte taxon

A small subset containing 170 observations of one bryophyte taxon.

## Usage

``` r
bryophyta_one
```

## Format

A `data.frame` containing 170 example observations of one bryophyte
taxon described by 6 variables:

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
