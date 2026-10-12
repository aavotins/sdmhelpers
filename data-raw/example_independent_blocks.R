library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readxl)
library(blockCV)

blocks=read_rds("../00_EGVcreation/HiQBioDiv_CodePrivate/SpeciesModels/Bryophyta/ODODEN_v0c1_CEP/blockCV_ODODEN.RDS")
example_independent_blocks=blocks$plots$map_plot$data
example_independent_blocks=st_transform(example_independent_blocks,crs=3059)

save(example_independent_blocks,
     file = "data/example_independent_blocks.rda",
     compress = "xz")
usethis::use_data(example_independent_blocks,
                  overwrite = TRUE,
                  compress = "xz")
