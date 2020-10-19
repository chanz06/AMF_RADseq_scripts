####function for one.sided t.test compare y frequencies to x mean frequency for Ref Allele
site_test <- function(x,y,alpha=0.05){
  alpha=alpha
  refx <- x[,seq(5,ncol(x)-11,4)]
  altx <- x[,seq(6,ncol(x)-10,4)]
  musx <- x[,ncol(x)-7]
  refy <- y[,seq(5,ncol(y)-11,4)]
  alty <- y[,seq(6,ncol(y)-10,4)]
  musy <- y[,ncol(y)-7]
  test <- vector("list", nrow(x))
  results <- data.frame()
  results_tmp <- data.frame()
  for (j in 1:nrow(x)){
    refxvals <- refx[j,]
    refyvals <- refy[j,]
    mux <- musx[j]
    #### for t-test
    #### test[[j]] <- t.test(refyvals, refxvals)
    test[[j]] <- wilcox.test(as.matrix(refyvals),as.matrix(refxvals),correct = TRUE)
    results_tmp <- data.frame(df=test[[j]]$statistic, pvalue=test[[j]]$p.value)
    #### for t-test
    #### results_tmp <- data.frame(df=test[[j]]$parameter, pvalue=test[[j]]$p.value)
    results <- rbind(results,results_tmp)
  }
  rownames(results) <- rownames(x)
  sig_results <- results[which(results$pvalue < alpha),]
  return(sig_results)
}