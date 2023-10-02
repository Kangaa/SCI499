library(tidyverse)
library(sf)
library(glue)


# Chapter 4

sim_dirs <- list.files("../../data/sims/new", recursive = TRUE)

tot_dirs <- sim_dirs[!str_detect(sim_dirs, "patchinf|patchsus|summary")] %>%
  paste0("../../data/sims/new/", .)

sims <- read_csv(tot_dirs, col_names = c("t", "S", "I"), skip = 1, id = "sim" )

sims <-  sims |>
  filter(t != "Int64", grepl("HPMM", sim)) %>%
  mutate(
    S = as.numeric(S),
    I = as.numeric(I),
    t = as.numeric(t),
    SA = str_extract(sim, "SA./SA(.)", 1),
    Mixmat = str_extract(sim, "SA./SA._..._..._(HPMM|HMM|OD)_", 1),
    mu = str_extract(sim, "SA./SA._..._..._.{3,4}_(.{0,3})_", 1),
    Intervention = factor(str_extract(sim,"SA./SA._..._..._.{3,4}_.{0,3}_(none|local|travel|total)", 1),levels = c("none", "local", "travel", "total"), ordered = TRUE),
    beta = as.numeric(str_extract(sim, "SA./SA._(...)", group = 1)),
    gamma = as.numeric(str_extract(sim, "SA./SA._..._(...)", group = 1)),
    R0 = beta/gamma,
    sim_num = str_extract(sim, "(_.{1,3}).csv", 1))


## Population proportional SA mixing matrix


sim_stats <- sims %>%
  group_by(sim_num, SA, Mixmat,mu,  Intervention, R0) %>%
  summarise(peaksize = max(I),
            finalsize = max(S) - min(S),
            duration = max(t))
## peak size

sim_stats %>%
  filter(Intervention == "travel", SA == "3") %>%
  ggplot() +
  aes(x = R0, y = peaksize) +
  geom_point(aes(col = mu)) +
  facet_grid(SA ~ Intervention)


## Final Size


## Duration

### SA3


### SA4




