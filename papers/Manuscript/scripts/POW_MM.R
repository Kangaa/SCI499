library(tidyverse)

SA3_OD <- read_csv("../../data/GMelb_UR-POW_2016.csv")
SA3_pop <- read_csv("../../data/GMelbSA3Pop21.csv")


colnames(SA3_OD)[-1] == SA3_OD[1]

SA3_OD <- SA3_OD %>%
  select(-1) %>%
  `colnames<-`(as.character(SA3_pop$SA3_Code[match(SA3_pop$SA3_Name, SA3_OD[[1]])]))


MM <- SA3_OD  %>%
  rowwise() %>%
  mutate(sm = sum(c_across(everything()))) %>%
  mutate(across(-sm, function(x) x/sm)) %>%
  select(-sm)

write.csv(MM, "../../data/Moss_new/GMelb_OD_SA3.csv")

