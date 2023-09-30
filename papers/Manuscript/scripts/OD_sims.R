#Metapop Plots
library(tidyverse)

sim_dirs <- list.files("../../data/sims", recursive = TRUE)

tot_dirs <- sim_dirs[str_detect(sim_dirs, "OD")] %>%
  paste0("../../data/sims/", .)

Summary_dirs <- tot_dirs[!str_detect(tot_dirs, "patchinf|patchsus|summary")]

sims <- read_csv(Summary_dirs, col_names = c("t", "S", "I"), skip = 1, id = "sim")

sims <-  sims |>
  filter() %>%
  mutate(
    S = as.numeric(S),
    I = as.numeric(I),
    t = as.numeric(t),
    SA = str_extract(sim, "SA./SA(.)", 1),
    Mixmat = str_extract(sim, "SA./SA._..._..._(HPMM|HMM|OD)_", 1),
    mu = str_extract(sim, "SA./SA.*_OD_(...)", 1),
    Intervention = factor(str_extract(sim,"/SA.*(none|local|travel|total)", 1),levels = c("none", "local", "travel", "total"), ordered = TRUE),
    beta = as.numeric(str_extract(sim, "SA./SA._(...)", group = 1)),
    gamma = as.numeric(str_extract(sim, "SA./SA._..._(...)", group = 1)),
    R0 = beta/gamma,
    sim_num = str_extract(sim, "_(.{1,3}).csv", 1))%>%
  select(-sim) %>%
  left_join(GMelb_Popns, by = join_by("patch" == "SA3_Name")) %>%
  pivot_wider( names_from = Measure, values_from = Count) %>%
  mutate(cuminf = SA3_pop - sus,
         prop_inf = inf/SA3_pop)

sim_stats <- sims %>%
  group_by(sim_num, SA, Mixmat,mu,  Intervention, R0) %>%
  summarise(peaksize = max(I),
            finalsize = max(S) - min(S),
            duration = max(t))



SA3_OD_Patch_dirs <- list.files("../../data/sims/SA3", pattern = "OD_..._none_(patchinf|patchsus)_([0-9]|10).csv", full.names = TRUE)

SA3_OD_Patch_logs <- read_csv(SA3_OD_Patch_dirs, id = "sim")

GMelb_Popns <- read_csv("../../data/GmelbSA3Pop21.csv")

SA3_patch_sims <- SA3_OD_Patch_logs %>%
  group_by(sim) %>%
  mutate(t = row_number()) %>%
  ungroup() %>%
  pivot_longer(cols = c(-sim, -t), values_to = "Count", names_to = "patch") %>%
  mutate(
         SA = str_extract(sim, "SA./SA(.)", 1),
         Mixmat = str_extract(sim, "SA./SA._..._..._(HPMM|HMM|OD)_", 1),
         mu = as.numeric(str_extract(sim, "SA./SA.*_OD_(...)", 1)),
         Intervention = "none",
         beta = as.numeric(str_extract(sim, "SA./SA._(...)", group = 1)),
         gamma = as.numeric(str_extract(sim, "SA./SA._..._(...)", group = 1)),
         R0 = beta/gamma,
         sim_num = str_extract(sim, "_(.{1,3}).csv", 1),
         Measure = str_extract(sim, "(inf|sus)_.{1,3}.csv", 1),
         ) %>%
  select(-sim) %>%
  left_join(GMelb_Popns, by = join_by("patch" == "SA3_Name")) %>%
  pivot_wider( names_from = Measure, values_from = Count) %>%
  mutate(cuminf = SA3_pop - sus,
         prop_inf = inf/SA3_pop)

patch_summary <- SA3_patch_sims %>%
  group_by(patch, R0, Intervention, mu) %>%
  summarise(
    peak = max(prop_inf),
    final_size = (first(SA3_pop) - min(sus))/first(SA3_pop),
    peak_time = t[which.max(inf)])

## peak size

patch_summary %>%
  filter(R0 == 2) %>%
  ggplot() +
  aes(x = patch, y = peak) +
  geom_point(aes(col = mu)) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0))

##Peaktime

patch_summary %>%
  ggplot() +
  aes(x = patch, y = peak_time) +
  geom_point(aes(col = R0)) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0))

#final_size

patch_summary %>%
  ggplot() +
  aes(x = patch, y = final_size) +
  geom_point(aes(col = R0)) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0))

Patch stats

{r, eval=FALSE}
patch_summary <- SA3_patch_sims %>%
  group_by(patch, R0, Intervention, mu) %>%
  summarise(
    peak = max(prop_inf),
    final_size = (first(SA3_pop) - min(sus))/first(SA3_pop),
    peak_time = t[which.max(inf)])

peak size

{r}
patch_summary %>%
  filter(R0 == 2) %>%
  ggplot() +
  aes(x = patch, y = peak) +
  geom_point(aes(col = mu)) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0))

##Peaktime

{r}
patch_summary %>%
  filter(R0 == 2) %>%
  ggplot() +
  aes(x = patch, y = peak_time) +
  geom_point(aes(col = mu)) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0))

#final_size

{r}
patch_summary %>%
  filter(R0 == 2) %>%
  ggplot() +
  aes(x = patch, y = final_size) +
  geom_point(aes(col = mu)) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0))


