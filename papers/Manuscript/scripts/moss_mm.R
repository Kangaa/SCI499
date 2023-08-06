## recreating Moss et al. 2018 mixing matrices
library(tidyverse)
moss_mm <- read_csv("../../data/Moss2018/mixing/abs_all.csv", col_select = -1) |>
  as.matrix()

delta_h = c((1:19)/20)
delta_a = 1 - delta_h
mms = replicate(length(delta_h), moss_mm, simplify = FALSE)

for (d in seq_along(delta_h)){
  for (i in 1:nrow(moss_mm)){
    for (j in 1:ncol(moss_mm)){
      ifelse(i == j,
             mms[[d]][i, j] <- delta_h[d],
             mms[[d]][i, j] <-  delta_a[d] * mms[[d]][i, j])
      }
  }
}

for (i in seq_along(mms)){

  write_csv(mms[i] |> as.data.frame(), paste0("../../data/Moss_new/Moss_MM_", i, ".csv"))}

