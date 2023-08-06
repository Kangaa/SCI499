library(tidyverse)
library(sf)

SA2_shapefile <- read_sf("../../data/ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp")


melb_shp_SA2 <- SA2_shapefile |>
  filter(GCC_CODE21 == "2GMEL")


melb_shp_SA3 |>
  group_by(SA4_CODE21) |>
  summarise(geom = st_union(geometry)) |>
  plot()

