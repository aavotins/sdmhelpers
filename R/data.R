#' Example amphibian observations with dates
#'
#' A reduced set of amphibian observations with dates but no locations.
#'
#' @format A `data.frame` containing 224484 example observations described by 4 variables:
#'
#' - `species`: scientific name of observed species.
#'
#' - `code`: abbreviation of the scientific name.
#'
#' - `year`: year of observation.
#'
#' - `date`: date of observation.
#'
#' @source A reduced subset prepared for package examples.
"amphibia_dates"


#' Large (species-wise) example bird-observation dataset
#'
#' A reduced bird-observation dataset intended for examples that require
#' multiple observations, species, dates, or times. Contains 132647 observations
#' from 200 species.
#'
#' @format A `data.frame` containing 132647 example observations from 200 species
#' described by 7 variables:
#'
#' - `species`: pseudonimised species.
#'
#' - `coordX`: approximate X-coordinate of observation (EPSG: 3059).
#'
#' - `coordY`: approximate Y-coordinate of observation (EPSG: 3059).
#'
#' - `year`: the year of observation.
#'
#' - `date`: the date (`as.Date` object) of the observation. Contains missing values.
#'
#' - `date_time`: the `POSIXct` object of date-time. Contains missing values.
#'
#' - `id100`: identificator of 100 m grid cell for joining with `data(points100)`.
#'
#' @source A reduced and pseudonimised subset prepared for package examples.
"aves_large"


#' Small (species-wise) example bird-observation dataset
#'
#' A small bird-observation dataset intended for fast examples. Contains 6539 observations
#' from 5 species.
#'
#' @format A `data.frame` containing 6359 example observations from 5 species
#' described by 6 variables:
#'
#' - `species`: scientific name of observed species.
#'
#' - `code`: abbreviation of the scientific name.
#'
#' - `year`: the year of observation.
#'
#' - `date`: the date (`as.Date` object) of the observation.
#'
#' - `date_time`: the `POSIXct` object of date-time. Contains missing values.
#'
#' - `id100`: identificator of 100 m grid cell for joining with `data(points100)`.
#'
#' @source A reduced subset prepared for package examples.
"aves_small"


#' Example bryophyte observations
#'
#' A reduced dataset containing 12587 observations of 551 bryophyte taxa.
#'
#' @format A `data.frame` containing 12587 example observations from 551 species
#' described by 6 variables:
#'
#' - `species`: pseudonimised species.
#'
#' - `coordX`: approximate X-coordinate of observation (EPSG: 3059).
#'
#' - `coordY`: approximate Y-coordinate of observation (EPSG: 3059).
#'
#' - `year`: the year of observation.
#'
#' - `date`: the date (`as.Date` object) of the observation. Contains missing values.
#'
#' - `id100`: identificator of 100 m grid cell for joining with `data(points100)`.
#'
#' @source A reduced and pseudonimised subset prepared for package examples.
"bryophyta_all"


#' Example observations for one bryophyte taxon
#'
#' A small subset containing 170 observations of one bryophyte taxon.
#'
#' @format A `data.frame` containing 170 example observations of one bryophyte taxon
#' described by 6 variables:
#'
#' - `species`: pseudonimised species.
#'
#' - `coordX`: approximate X-coordinate of observation (EPSG: 3059).
#'
#' - `coordY`: approximate Y-coordinate of observation (EPSG: 3059).
#'
#' - `year`: the year of observation.
#'
#' - `date`: the date (`as.Date` object) of the observation. Contains missing values.
#'
#' - `id100`: identificator of 100 m grid cell for joining with `data(points100)`.
#'
#' @source A reduced and pseudonimised subset prepared for package examples.
"bryophyta_one"


#' Names of example ecogeographical variables
#'
#' Descriptive names and abbreviations corresponding to the ecogeographical
#' variable layers supplied with the package examples.
#'
#' @format A `data.frame` describing 5 ecogeographical (environmental) variables
#' provided in `"inst/extdata/example_egv.tif"` with 3 fields:
#'
#' - `layername`: abbreviated layer name used for `terra::names()`.
#'
#' - `longname_english`: descriptive name of the variable in English.
#'
#' - `longname_latvia`: descriptive name of the variable in Latvian.
#'
"egv_names"


#' Example background locations
#'
#' Spatial background locations used in fold-selection and model-evaluation
#' examples.
#'
#' @format A `data.frame` containing 40000 background location coordinates
#' with 2 fields:
#'
#' - `x`: approximate X-coordinate (EPSG: 3059).
#'
#' - `y`: approximate Y-coordinate (EPSG: 3059).
#'
#' @source Prepared for package examples and practically used in species
#' distribution modelling.
"example_background"


#' Example raster-grid geometry
#'
#' A raster-grid object defining the spatial geometry, pixel alignment and CRS
#' used by the package's raster examples.
#'
#' @format A `SpatRaster` object created with [terra::rast()].
#'
#' @details
#' File contains only presence cells coded with value `1`.
#'
#' @source Derived from a reduced section of the Latvian 100-m raster grid.
"example_grid"


#' Example independent spatial blocks
#'
#' Fold or block identifiers used to demonstrate independent-data selection.
#'
#' @format A `sf` object containing independent spatial-blocks as `POLYGON` and
#' their identifiers.
#'
#' @source A reduced subset prepared for package examples from actual species
#' distribution modelling.
"example_independent_blocks"


#' Example species-presence locations
#'
#' Spatial presence records used in spatial KDE, fold-selection, and model
#' evaluation examples.
#'
#' @format A `data.frame` containing 2778 species presence location coordinates
#' with 2 fields:
#'
#' - `x`: approximate X-coordinate (EPSG: 3059).
#'
#' - `y`: approximate Y-coordinate (EPSG: 3059).
#'
#' @source Prepared for package examples and practically used in species
#' distribution modelling.
"example_presences"


#' Example habitat inventories
#'
#' Habitat-inventory observations used to demonstrate seasonal weighting
#' against species-activity dates.
#'
#' @format A data frame or spatial object containing inventory dates covering
#' 127700 `data(example_raster)` centroids and 6 associated attributes:
#'
#' - `habitat_type`: character vector of bvroad habitat class:
#'
#'  - `dune`;
#'
#'  - `forest`;
#'
#'  - `freshwater`;
#'
#'  - `grassland`;
#'
#'  - `heath`;
#'
#'  - `mire`.
#'
#' - `coordX`: X-coordinate of `data(example_raster)` cell's centre, matching `data(points100)` object (EPSG: 3059).
#'
#' - `coordY`: Y-coordinate of `data(example_raster)` cell's centre, matching `data(points100)` object (EPSG: 3059).
#'
#' - `year`: year of habitat inventory.
#'
#' - `date`: date (`as.Date()`) of habitat inventory
#'
#' - `id100`: identificator of 100 m grid cell for joining with `data(points100)`.
#'
#' @source A reduced and anonymised subset prepared for package examples.
"habitat_inventories"


#' Example 100-m grid points
#'
#' Point centres of a reduced 100-m grid used in spatial examples.
#'
#' @format An `sf` object containing point geometries and grid identifiers:
#'
#' - `id`: identificator of 100 m grid cell.
#'
#' - `yes`: indicator for location in terrestrial Latvia (all the cells contain value `1`).
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `rinda300`: identificator of 300 m grid cell.
#'
#' - `ID1km`: identificator of 1 km grid cell.
#'
#' - `rinda500`: identificator of 500 m grid cell.
#'
#' - `geom`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 100-m analysis grid.
"points100"
