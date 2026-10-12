library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readxl)

bg_test=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/Bryophyta/ODODEN_v0c1_CEP/PointsBgTest_ODODEN.csv")
bg_train=read_csv("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/Bryophyta/ODODEN_v0c1_CEP/PointsBgTrain_ODODEN.csv")

example_background=rbind(bg_test,bg_train)

save(example_presences,
     file = "data/example_background.rda",
     compress = "xz")
usethis::use_data(example_background,
                  overwrite = TRUE,
                  compress = "xz")
