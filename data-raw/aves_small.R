library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readr)

plesi=read_rds("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/00_Observations/BirdObs_Mislocs.RDS")
table(plesi$avots)
plesi=plesi |>
  filter(avots=="monitorings_Plēšputnu") |>
  filter(KODS %in% c("BUTBUT","CIRAER","GLAPAS","STRURA","STRALU")) |>
  sf::st_as_sf(coords=c("lksX","lksY"),remove=FALSE,crs=3059)
aves_small=data.frame(plesi) |>
  mutate(species=zinatniski) |>
  mutate(code=KODS) |>
  mutate(coordX=round(lksX/100,0)*100+1) |>
  mutate(coordY=round(lksY/100,0)*100+1) |>
  mutate(year=year) |>
  mutate(date=datums_no) |>
  mutate(date_time=datums_laiks) |>
  dplyr::select(species,code,year,date,date_time,id100)


save(aves_small,
     file = "data/aves_small.rda",
     compress = "xz")
usethis::use_data(aves_small,
                  overwrite = TRUE,
                  compress = "xz")
