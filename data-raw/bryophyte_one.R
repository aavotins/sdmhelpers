library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readr)

tks50=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGrids/tks93_50km.parquet")
telpa=tks50 |> filter(NOSAUKUMS=="Jūrmala")
telpa=st_transform(telpa,crs=3059)
telpa_buff=st_buffer(telpa,dist=20000)

bryophyte_one=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/00_Observations/Observations_v14_RS_tikls.csv")
bryophyte_one=bryophyte_one |>
  filter(CODE=="ODODEN") |>
  sf::st_as_sf(coords=c("LKS_X","LKS_Y"),remove=FALSE,crs=3059)
bryophyte_one=bryophyte_one[telpa_buff,,]
bryophyta_one=data.frame(bryophyte_one) |>
  mutate(species="Bryophyta_reference") |>
  mutate(coordX=round(LKS_X/100,0)*100+1) |>
  mutate(coordY=round(LKS_Y/100,0)*100+1) |>
  mutate(year=Year) |>
  mutate(date=Date) |>
  dplyr::select(species,coordX,coordY,year,date,id100)


save(bryophyta_one,
     file = "data/bryophyta_one.rda",
     compress = "xz")
usethis::use_data(bryophyta_one,
                  overwrite = TRUE,
                  compress = "xz")
