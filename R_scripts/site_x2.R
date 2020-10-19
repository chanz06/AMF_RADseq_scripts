#### Function to do chisq.test on allele counts of Ref and Alt alleles in at a given alpha
site_x2 <- function(x,y,alpha=0.05){
  i=0
  sig_list=NULL
  sig_temp=NULL
  pvals=NULL
  #j=0
  alpha=alpha
  dat1 <- x[,c(ncol(x)-8,ncol(x)-3)]
  dat2 <- y[,c(ncol(y)-8,ncol(y)-3)]
  for(i in 1:nrow(dat1)){
    t1 <- dat1[i,]
    t2 <- dat2[i,]
    rownames(t1) <- paste0("X_",rownames(t1))
    rownames(t2) <- paste0("Y_", rownames(t2))
    tab <- rbind(t1,t2)
    X2 <- chisq.test(tab, correct = TRUE)
    p <- X2$p.value
    sig_temp <- data.frame(rows=rownames(x[i,]), pvalue=X2$p.value)
    sig_list <- rbind(sig_list,sig_temp)
  }
  rownames(sig_list) <- sig_list[,1]
  sig_list <- sig_list[which(sig_list$pvalue < alpha),]
  return(sig_list)
}
