library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readr)

tks50=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGrids/tks93_50km.parquet")
telpa=tks50 |> filter(NOSAUKUMS=="Jūrmala")
telpa=st_transform(telpa,crs=3059)
telpa_buff=st_buffer(telpa,dist=20000)

modelejamas_sunas=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/00_Observations/Observations_v14_RS_tikls.csv")
extra_sunas=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/00_Observations/ObsEffort_Bryophyta.csv")

names(modelejamas_sunas)
modelejamas_sunas=modelejamas_sunas |>
  filter(CODE!="ODODEN") |>
  sf::st_as_sf(coords=c("LKS_X","LKS_Y"),remove=FALSE,crs=3059)
modelejamas_sunas=modelejamas_sunas[telpa_buff,,]
modelejamas_sunas=data.frame(modelejamas_sunas) |>
  mutate(coordX=round(LKS_X/100,0)*100+1) |>
  mutate(coordY=round(LKS_Y/100,0)*100+1) |>
  mutate(year=Year) |>
  mutate(date=Date) |>
  dplyr::select(CODE,coordX,coordY,year,date,id100)

names(extra_sunas)
extra_sunas=extra_sunas |>
  filter(CODE!="ODODEN") |>
  sf::st_as_sf(coords=c("LKS_X","LKS_Y"),remove=FALSE,crs=3059)
extra_sunas=extra_sunas[telpa_buff,,]
extra_sunas=data.frame(extra_sunas) |>
  mutate(coordX=round(LKS_X/100,0)*100+1) |>
  mutate(coordY=round(LKS_Y/100,0)*100+1) |>
  mutate(year=Year) |>
  mutate(date=Date) |>
  dplyr::select(CODE,coordX,coordY,year,date,id100)



sunas=rbind(modelejamas_sunas,extra_sunas)
bryophyta_all=data.frame(sunas) |>
  mutate(kodiem=as.numeric(as.factor(CODE))+1,
         species=paste0("Bryophyta_",kodiem)) |>
  dplyr::select(species,coordX,coordY,year,date,id100)


save(bryophyta_all,
     file = "data/bryophyta_all.rda",
     compress = "xz")
usethis::use_data(bryophyta_all,
                  overwrite = TRUE,
                  compress = "xz")
