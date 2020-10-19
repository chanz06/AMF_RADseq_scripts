## Function to filter based on 66% (two-thirds) consistency among replicates and replace possible NA values
library("matrixStats")
filter_df <- function(df,fam_name){
  fam_name <- fam_name
  #samp <- deparse(substitute(df))
  n <- ncol(df)
  s <- ((ncol(df)-2)/2)
  odds <- seq(5,ncol(df)-1,4)
  oddC <- seq(3,ncol(df)-1,4)
  evens <- seq(6,ncol(df),4)
  evenC <- seq(4,ncol(df),4)
  #odds <- seq(3,ncol(df)-1,2)
  #evens <- seq(4,ncol(df),2)
  df$NA_count <- apply(df[,3:n], 1, function(x) sum(is.na(x)))
  df <- df[which(df$NA_count <=round(((ncol(df)-2)*0.25),digits = 0)),-ncol(df)]
  df[,odds] <- (is.na(df[,odds]))*rowMedians(as.matrix(df[,odds]), na.rm=TRUE)[row(df[,odds])] + replace(df[,odds], is.na(df[,odds]), 0)
  df <- df[which(round(rowMedians(as.matrix(df[,odds])),3) >= 0.1),]
  df$RC_Total <- rowSums(df[,oddC],na.rm = TRUE)
  df$RC_Median <- round(rowMedians(as.matrix(df[,oddC]),na.rm = TRUE),digits = 0)
  df$RF_Median <- round(rowMedians(as.matrix(df[,odds])),3)
  df$RF_SD <- round(apply(df[,odds],1,sd),3)
  df$RF_SE <- round(df$RF_SD/s,3)
  df[,evens] <- (is.na(df[,evens]))*rowMedians(as.matrix(df[,evens]), na.rm=TRUE)[row(df[,evens])] + replace(df[,evens], is.na(df[,evens]), 0)
  df <- df[which(round(rowMedians(as.matrix(df[,evens])),3) >= 0.1),]
  df$AC_Total <- rowSums(df[,evenC],na.rm = TRUE)
  df$AC_Median <- round(rowMedians(as.matrix(df[,evenC]),na.rm = TRUE),digits = 0)
  df$AF_Median <- round(rowMedians(as.matrix(df[,evens])),3)
  df$AF_SD <- round(apply(df[,evens],1,sd),3)
  df$AF_SE <- round(df$AF_SD/s,3)
  df <- df[!duplicated(df[,1:2]),]
  #df <- df[which(df$RF_SD <= 0.15),]
  rownames(df) <- paste0(df$Sca,"_" ,df$Pos)
  return(df)
}