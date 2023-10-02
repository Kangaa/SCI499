
library(tidyverse)
library(sf)
library(glue)

# Chapter 3

#| echo: false
#| warning: false
levels <- list(SA2 = "SA2", SA3 = "SA3", SA4 = "SA4")
Shapes <- map(
  map(
    levels,
    \(x) glue("../../data/ASGS_GDA2020/{x}_2021_AUST_SHP_GDA2020/{x}_2021_AUST_GDA2020.shp")),
  read_sf)


GMelb_Shapes <- Shapes |>
  map(\(x) filter(x, GCC_NAME21 == "Greater Melbourne"))



plotmixmat <- function(mixmat){
  mixmat %>%
    as_tibble(rownames = "i") %>%
    pivot_longer(-1, names_to = "j", values_to = "mixing") %>%
    ggplot() +
    aes(x = j, y = i, fill = log(mixing)) +
    geom_tile(show.legend = FALSE) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())+
    scale_fill_continuous(type = "viridis")}



GMelb_Popns <- map(
  map(
    c("2", "3", "4"),
    \(x) str_replace("../../data/GmelbSA_X_Pop21.csv", "_X_", x)),
  read_csv) |>
  `names<-`(c("SA2", "SA3", "SA4"))



MMpops <- map(c("SA2", "SA3", "SA4"),
              \(x){
                size <- nrow(GMelb_Popns[[x]])

                MMpop <- matrix(NA,
                                size,
                                size,
                                dimnames = list(
                                  GMelb_Popns[[x]][[paste0(x, "_Code")]],
                                  GMelb_Popns[[x]][[paste0(x, "_Code")]]))

                Ntot = sum(GMelb_Popns[[x]][[3]])

                for (i in seq_along(MMpop[1,])){
                  for (j in seq_along(MMpop[,1])){
                    MMpop[i,j] <- GMelb_Popns[[x]][[3]][j]/Ntot
                  }
                }
                MMpop
              })



Hmixmat <- function(data){

  levels <- data %>%
    names() %>%
    subset(grepl("SA._CODE|GCC_CODE21", .))

  xi <- rep(1/length(levels), length(levels))

  npatch = nrow(data)

  level_vec = vector("character", npatch)

  mixmat = matrix(0, npatch, npatch,
                  dimnames = list(data[[levels[1]]],
                                  data[[levels[1]]]))

  for (i in 1:npatch){
    patch_i <- data[i,]
    norm <- rep(0, length(levels))
    for (j in 1:npatch){
      patch_j <- data[j,]
      for (l in seq_along(levels)){ #find lowest level of association (LLA)
        code_l <-  levels[l]
        if (patch_i[[code_l]] == patch_j[[code_l]]){
          level_vec[j] <- levels[l] #place LLA in level_vec(i)
          norm[l:length(norm)] <- norm[l:length(norm)] + 1 #add +1 to the count of patches in the same level L as i
          break
        }
      }
    }
    for (l in seq_along(levels)){
      for (j in 1:npatch){
        #if (level_vec[j] %in% levels[1:l]){
        if (level_vec[j] == levels[l]){
          #mixmat[i,j] <- mixmat[i,j] + (xi[l]/norm[l])
          if (l == 1) mixmat[i,j] <- xi[l]/norm[l]
          else mixmat[i,j] <- xi[l]/(norm[l] - norm[l-1])
        }
      }
    }
  }
  mixmat
}
Hmixmat(GMelb_MM_eg)

HPmixmat <- function(data){
  data <- data %>% tibble()
  levels <- data %>%
    names() %>%
    subset(grepl("SA._CODE|GCC_CODE21", .))

  xi <- rep(1/length(levels), length(levels))

  npatch = nrow(data)

  level_vec = vector("character", npatch)

  mixmat = matrix(0, npatch, npatch,
                  dimnames = list(data[[levels[1]]],
                                  data[[levels[1]]]))
  for (i in 1:npatch){
    patch_i <- data[i,]
    norm <- rep(0, length(levels))
    for (j in 1:npatch){
      patch_j <- data[j,]
      for (l in seq_along(levels)){ #find lowest level of association (LLA)
        code_l <-  levels[l]
        if (patch_i[[code_l]] == patch_j[[code_l]]){
          #mixmat[i,j] <- mixmat[i,j] + (xi[l]/norm[l])
          if (l == 1) mixmat[i,j] <- xi[l]


          else {
            ll_pop <- GMelb_Popns[[1]][GMelb_Popns[[1]][[glue("{substr(levels[1], 1,3)}_Code")]] == patch_i[[levels[1]]],][[3]]
            ul_pop <- GMelb_Popns[[l]][GMelb_Popns[[l]][[glue("{substr(levels[l], 1,3)}_Code")]] == patch_j[[code_l]],][[3]]

            mixmat[i,j] <- xi[l]*(ul_pop/ll_pop)}
        }
      }
    }
  }
  mixmat
}
HPmixmat(GMelb_MM_eg)
