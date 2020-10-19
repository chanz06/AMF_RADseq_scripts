## Function to collapse dataframes by uniqe site, select SNPs and compute allele frequencies
new_df <- function(df,samp){
  samp <- samp
  df %>%
    group_by(Sca,Pos,Ref,Type) %>%
    summarise_all(sum) %>%
    data.frame() -> df_new
  df_new <- df_new[which(df_new$Type %in% "snp"),]
  df_new[,8] <- c(df_new[,6]/df_new[,5])
  df_new[,9] <- c(df_new[,7]/df_new[,5])
  names(df_new)[8] <- paste0(samp, "_RF")
  names(df_new)[9] <- paste0(samp, "_AF")
  return(df_new)
}