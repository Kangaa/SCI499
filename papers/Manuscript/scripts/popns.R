#script to generate population summaries from mesh blocks data
library(tidyverse)
library(readxl)
MB_data <- excel_sheets("../../data/ASGS_GDA2020/Mesh Block Counts, 2021.xlsx") |>
  map(\(x) read_excel("../../data/ASGS_GDA2020/Mesh Block Counts, 2021.xlsx", x,skip = 5))


MB_data <- MB_data[2:13] %>%
  tibble() %>%
  unnest()

Hierarchy <- read_excel("../../data/ASGS_GDA2020/MB_2021_AUST.xlsx") %>%
  filter(GCCSA_CODE_2021 == )

MB_data <- full_join(MB_data,Hierarchy, by = "MB_CODE_2021")

levels = list(SA1 = "SA1",
              SA2 = "SA2",
              SA3 = "SA3",
              SA4 = "SA4")


MB_data %>%
  group_by(SA1_CODE_2021) %>%
  summarise(
    Popn = sum(Person)
  ) %>%
  write_csv("../../data/SA1_Popns_AUS.csv")


MB_data %>%
  group_by(SA2_CODE_2021) %>%
  summarise(
    Names = unique(SA2_NAME_2021),
    Popn = sum(Person)
  ) %>%
  write_csv("../../data/SA2_Popns_AUS.csv")

MB_data %>%
  group_by(SA3_CODE_2021) %>%
  summarise(
    Names = unique(SA3_NAME_2021),
    Popn = sum(Person)
  ) %>%
  write_csv("../../data/SA3_Popns_AUS.csv")

MB_data %>%
  group_by(SA4_CODE_2021) %>%
  summarise(
    Names = unique(SA4_NAME_2021),
    Popn = sum(Person)
  ) %>%
  write_csv("../../data/SA4_Popns_AUS.csv")

MB_data %>%
  group_by(GCCSA_CODE_2021) %>%
  summarise(
    Names = unique(GCCSA_NAME_2021),
    Popn = sum(Person)
  ) %>%
  write_csv("../../data/GCCSA_Popns_AUS.csv")

