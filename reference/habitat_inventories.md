# Example habitat inventories

Habitat-inventory observations used to demonstrate seasonal weighting
against species-activity dates.

## Usage

``` r
habitat_inventories
```

## Format

A data frame or spatial object containing inventory dates covering
127700 `data(example_raster)` centroids and 6 associated attributes:

- `habitat_type`: character vector of bvroad habitat class:

- `dune`;

- `forest`;

- `freshwater`;

- `grassland`;

- `heath`;

- `mire`.

- `coordX`: X-coordinate of `data(example_raster)` cell's centre,
  matching `data(points100)` object (EPSG: 3059).

- `coordY`: Y-coordinate of `data(example_raster)` cell's centre,
  matching `data(points100)` object (EPSG: 3059).

- `year`: year of habitat inventory.

- `date`: date ([`as.Date()`](https://rdrr.io/r/base/as.Date.html)) of
  habitat inventory

- `id100`: identificator of 100 m grid cell for joining with
  `data(points100)`.

## Source

A reduced and anonymised subset prepared for package examples.
