library(tidyverse)

SA2_files <- list.files("C:/Users/reubender/Downloads/SA3", full.names = TRUE)

tot_files <- SA2_files[!str_detect(SA2_files, "patchinf|patchsus|summary")]

sims <- read_csv(tot_files, col_names = c("t", "S", "I"), skip = 2, id = "sim")

sims <- sims |>
  mutate(
    beta = as.numeric(str_extract(sim, "(...)_.{3}_.{3}_.{1,3}.csv", group = 1)),
    gamma = as.numeric(str_extract(sim, "(...)_.{3}_.{1,3}.csv", group = 1)),
    mm = as.numeric(str_extract(sim, "(...)_.{1,3}.csv", group = 1)),
    R0 = beta/gamma)


sims |>
  ggplot() +
  aes(x = t, y = I, col = mm, group = sim) +
  geom_line() +
  facet_wrap(~R0, scales = "free")


## Final infection size
