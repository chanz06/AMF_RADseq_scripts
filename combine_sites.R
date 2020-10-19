## Function for combining sample allel counts and allel frequencies by site
library("purrr")
combine_sites <- function(dflist){
  n <- 0
  df <- NULL
  for (df in dflist) {
    fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(dflist[1])), ' ')[[1]][1])
    n <- n+1
    #  df <- get(i)
    #  df<- df[with(df, order(df[,1],df[,2])),]
    #  df <- df[df$Type=='snp',]
    if (n==1) {
      #fam <- df[,c(1:2,8:9)]
      fam <- df[,c(1:2,6:9)]
    }
    if (n>1) {
      #fam <- list(fam[,1:ncol(fam)], df[,c(1:2,8:9)]) %>% purrr::reduce(full_join, by = c("Sca", "Pos"), all = TRUE)
      fam <- list(fam[,1:ncol(fam)], df[,c(1:2,6:9)]) %>% purrr::reduce(full_join, by = c("Sca", "Pos"), all = TRUE)
      fam <- unique.data.frame(fam)
      fam <- fam[with(fam, order(fam[,1],fam[,2])),]
    }
    else {
      next
    }
  }
  return(assign(paste0(fam_name, "_fam"),fam))
}