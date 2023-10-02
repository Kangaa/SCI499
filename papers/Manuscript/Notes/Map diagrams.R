#| label: fig-MelbSAexaple
library(sf)
library(tidyverse)
Shapes <- list(SA1 = tibble(),
               SA2 = tibble(),
               SA3 = tibble(),
               SA4 = tibble())

Shapes$SA1 <- read_sf("../../data/ASGS_GDA2020/SA1_2021_AUST_SHP_GDA2020/SA1_2021_AUST_GDA2020.shp")
Shapes$SA2 <- read_sf("../../data/ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp")
Shapes$SA3 <- read_sf("../../data/ASGS_GDA2020/SA3_2021_AUST_SHP_GDA2020/SA3_2021_AUST_GDA2020.shp")
Shapes$SA4 <- read_sf("../../data/ASGS_GDA2020/SA4_2021_AUST_SHP_GDA2020/SA4_2021_AUST_GDA2020.shp")

GMelb_Shapes <- Shapes |>
  map(\(x) filter(x, GCC_NAME21 == "Greater Melbourne"))

GMelb_SA4_plot <- ggplot() + geom_sf(data = GMelb_Shapes$SA4, aes(fill = SA4_CODE21))

GMelb_SA4_colors <- ggplot_build(GMelb_SA4_plot) |>
  pluck("data",1, "fill")


ggplot() +
  geom_sf(data = GMelb_Shapes$SA4, aes(fill = SA4_CODE21), alpha = 0.5, linewidth = NA) +
  geom_sf(data = GMelb_Shapes$SA3 |> filter(SA3_NAME21 == "Maribyrnong"), fill = GMelb_SA4_colors[8], alpha = 0.5) +
  geom_sf(data = GMelb_Shapes$SA2 |> filter(SA2_NAME21 == "Footscray"), fill = GMelb_SA4_colors[8], alpha = 0.5) +
  geom_sf(data = GMelb_Shapes$SA1 |>  filter(SA1_CODE21 == "21303134811"), fill = GMelb_SA4_colors[8], alpha = 0.5)+
  coord_sf(xlim = c(144.8, 145) , ylim = c(-37.75, -37.9))


ggplot() +
  geom_sf(data = GMelb_Shapes$SA4, aes(fill = SA4_CODE21), alpha = 0.5, linewidth = NA) +
  geom_sf(data = GMelb_Shapes$SA3, fill = NA) +
  geom_sf(data = GMelb_Shapes$SA2, lty = "dotted", fill = NA)



