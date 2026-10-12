library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readxl)

pres_test=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/Bryophyta/ODODEN_v0c1_CEP/PointsPresTest_ODODEN.csv")
pres_train=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/Bryophyta/ODODEN_v0c1_CEP/PointsPresTrain_ODODEN.csv")

example_presences=rbind(pres_test,pres_train)

save(example_presences,
     file = "data/example_presences.rda",
     compress = "xz")
usethis::use_data(example_presences,
                  overwrite = TRUE,
                  compress = "xz")
