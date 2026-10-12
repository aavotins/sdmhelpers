library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readr)

tks50=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGrids/tks93_50km.parquet")
telpa=tks50 |> filter(NOSAUKUMS=="Jūrmala")
telpa=st_transform(telpa,crs=3059)
telpa_buff=st_buffer(telpa,dist=20000)


aves=read_rds("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/00_Observations/BirdObs_Mislocs.RDS")
aves=aves |>
  sf::st_as_sf(coords=c("lksX","lksY"),remove=FALSE,crs=3059)
aves=aves[telpa_buff,,]
aves_large=data.frame(aves) |>
  mutate(kodiem=as.numeric(as.factor(KODS))+1,
         species=paste0("Aves_",kodiem)) |>
  mutate(coordX=round(lksX/100,0)*100+1) |>
  mutate(coordY=round(lksY/100,0)*100+1) |>
  mutate(year=year) |>
  mutate(date=datums_no) |>
  mutate(date_time=datums_laiks) |>
  dplyr::select(species,coordX,coordY,year,date,date_time,id100)


save(aves_large,
     file = "data/aves_large.rda",
     compress = "xz")
usethis::use_data(aves_large,
                  overwrite = TRUE,
                  compress = "xz")
