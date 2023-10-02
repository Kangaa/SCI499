### patchplots

library(tidyverse)

sim_dirs <- list.files("../../data/sims", recursive = TRUE)

patch_dirs <- sim_dirs[str_detect(sim_dirs, "patchinf|patchsus")] %>%
  paste0("../../data/sims/", .)

sims <- read_csv(patch_dirs, id = "sim" )

sims <-  sims |>
  filter(t != "Int64", grepl("HPMM", sim)) %>%
  mutate(
    I = as.numeric(I),
    t = as.numeric(t),
    SA = str_extract(sim, "SA./SA(.)", 1),
    Mixmat = str_extract(sim, "SA./SA._..._..._(HPMM|HMM|OD)_", 1),
    mu = str_extract(sim, "SA./SA._..._..._.{3,4}_(.{0,3})_", 1),
    Intervention = factor(str_extract(sim,"SA./SA._..._..._.{3,4}_.{0,3}_(none|local|travel|total)", 1)),
    beta = as.numeric(str_extract(sim, "SA./SA._(...)", group = 1)),
    gamma = as.numeric(str_extract(sim, "SA./SA._..._(...)", group = 1)),
    R0 = beta/gamma,
    sim_num = str_extract(sim, "(_.{1,3}).csv", 1))
