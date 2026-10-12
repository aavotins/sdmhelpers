library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readr)

abinieki=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/00_Observations/Observations_v14_RS_tikls.csv")

amphibia_dates=abinieki |>
  filter(Group=="Amphibia")
amphibia_dates=data.frame(amphibia_dates) |>
  mutate(coordX=round(LKS_X/100,0)*100+1) |>
  mutate(coordY=round(LKS_Y/100,0)*100+1) |>
  mutate(year=Year) |>
  mutate(date=Date) |>
  mutate(group=Group) |>
  mutate(species=ScientificName) |>
  mutate(code=CODE) |>
  dplyr::select(species,code,year,date)


save(amphibia_dates,
     file = "data/amphibia_dates.rda",
     compress = "xz")
usethis::use_data(amphibia_dates,
                  overwrite = TRUE,
                  compress = "xz")
