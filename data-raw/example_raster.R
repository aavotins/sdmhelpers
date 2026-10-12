library(terra)
library(sf)
library(sfarrow)
library(dplyr)

template100=rast("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateRasters/LV100m_10km.tif")
tks50=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGrids/tks93_50km.parquet")
telpa=tks50 |> filter(NOSAUKUMS=="Jūrmala")
telpa=st_transform(telpa,crs=3059)
telpa_buff=st_buffer(telpa,dist=20000)
example_grid=crop(template100,telpa_buff)
plot(example_grid)
example_grid=terra::wrap(example_grid)

save(example_grid,
     file = "data/example_grid.rda",
     compress = "xz")
usethis::use_data(example_grid,
                  overwrite = TRUE,
                  compress = "xz")
