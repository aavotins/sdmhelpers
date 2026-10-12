library(sf)
library(sfarrow)
library(dplyr)

points100=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGridPoints/pts100_sauzeme.parquet")
tks50=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGrids/tks93_50km.parquet")
telpa=tks50 |> filter(NOSAUKUMS=="Jūrmala")
#telpa=st_transform(telpa,crs=3059)
telpa_buff=st_buffer(telpa,dist=20000)
points100=points100[telpa_buff,,]
points100=st_transform(points100,crs=3059)

sf::st_geometry(points100) <- "geom"

attr(points100, "sf_column")
# should be "geom"

sf::st_geometry_type(points100)
sf::st_crs(points100)

save(points100,
     file = "data/points100.rda",
     compress = "xz")
usethis::use_data(points100,
                  overwrite = TRUE,
                  compress = "xz")
