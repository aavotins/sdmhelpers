library(terra)
library(sf)
library(sfarrow)
library(dplyr)
library(readxl)


nosaukumi=read_excel("../00_EGVcreation/HiQBioDiv_CodePrivate/RasterGrids_100m/2024/HiQBioDiv_EGVnosaukumi_20251202.xlsx")
nosaukumi=nosaukumi |>
  filter(layername %in% c("egv_280","egv_293","egv_302","egv_385","egv_400"))
egv_names=nosaukumi |>
  dplyr::select(layername,longname_english,longname_latvian)

save(egv_names,
     file = "data/egv_names.rda",
     compress = "xz")
usethis::use_data(egv_names,
                  overwrite = TRUE,
                  compress = "xz")


slani=terra::rast(paste0("../00_EGVcreation/HiQBioDiv_CodePrivate/RasterGrids_100m/2024/Scaled/",
                         nosaukumi$new_filename))
names(slani)=egv_names$layername
slani
plot(slani)


tks50=st_read_parquet("../00_EGVcreation/HiQBioDiv_CodePrivate/Templates/TemplateGrids/tks93_50km.parquet")
telpa=tks50 |> filter(NOSAUKUMS=="Jūrmala")
telpa=st_transform(telpa,crs=3059)
telpa_buff=st_buffer(telpa,dist=20000)
slani=crop(slani,telpa_buff)
plot(slani)

#example_egv <- round(slani, digits = 4)
example_egv=slani
#save(example_egv,
#     file = "data/example_egv.rda",
#     compress = "xz")
#usethis::use_data(example_egv,
#                  overwrite = TRUE,
#                  compress = "xz")

dir.create("inst/extdata",
           recursive = TRUE,
           showWarnings = FALSE)
terra::writeRaster(example_egv,
                   "inst/extdata/example_egv.tif",
                   overwrite = TRUE,
                   datatype = "FLT4S",
                   gdal = c("COMPRESS=DEFLATE",
                            "PREDICTOR=3",
                            "ZLEVEL=9"))
#egv_file <- system.file("extdata",
#                        "example_egv.tif",
#                        package = "sdmhelpers")
#example_egv <- terra::rast(egv_file)
