Hmixmat <- function(data){

  levels <- data %>%
    names() %>%
    subset(grepl("SA._CODE|GCC_CODE21", .))

  xi <- rep(1/length(levels), length(levels))

  npatch = nrow(data)

  level_vec = vector("character", npatch)

  mixmat = matrix(0, npatch, npatch,
                  dimnames = list(data[[glue("{levels[1]}_CODE21")]],
                                  data[[glue("{levels[1]}_CODE21")]]))

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
        if (level_vec[j] %in% levels[1:l]){
          mixmat[i,j] <- mixmat[i,j] + (xi[l]/norm[l])
        }
      }
    }
  }
  mixmat
}
