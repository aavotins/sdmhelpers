library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readr)

tks50=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGrids/tks93_50km.parquet")
telpa=tks50 |> filter(NOSAUKUMS=="Jūrmala")
telpa=st_transform(telpa,crs=3059)
telpa_buff=st_buffer(telpa,dist=20000)


habitats=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/00_effortEUhabitats/EUHabitats_points.parquet")
habitats=habitats |>
  sf::st_as_sf(coords=c("x","y"),remove=FALSE,crs=3059)
habitats=habitats[telpa_buff,,]
table(habitats$EUhabitat)
habitat_inventories=data.frame(habitats) |>
  mutate(habitat_type=case_when(EUhabitat=="kapas"~"dune",
                                EUhabitat=="koki"~"forest",
                                EUhabitat=="purvi"~"mire",
                                EUhabitat=="udeni"~"freshwater",
                                EUhabitat=="virsaji"~"heather",
                                EUhabitat=="zalaji"~"grassland")) |>
  mutate(coordX=x) |>
  mutate(coordY=y) |>
  mutate(year=lubridate::year(as.Date(OBS_DATE))) |>
  mutate(date=as.Date(OBS_DATE)) |>
  mutate(id100=id) |>
  dplyr::select(habitat_type,coordX,coordY,year,date,id100)


save(habitat_inventories,
     file = "data/habitat_inventories.rda",
     compress = "xz")
usethis::use_data(habitat_inventories,
                  overwrite = TRUE,
                  compress = "xz")
