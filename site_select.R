#### Function for selecting all sites to plotting figures
site_select <- function(x,y,sites){
  xname <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', colnames(x)[3]), ' ')[[1]][1])
  yname <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', colnames(y)[3]), ' ')[[1]][1])
  xsites <- x[which(rownames(x) %in% rownames(sites)),]
  xfreq <- xsites[,c(ncol(xsites)-7,ncol(xsites)-2)]
  colnames(xfreq) <- paste0(xname,"_",colnames(xfreq))
  ysites <- y[which(rownames(y) %in% rownames(sites)),]
  yfreq <- ysites[,c(ncol(ysites)-7,ncol(ysites)-2)]
  colnames(yfreq) <- paste0(yname,"_",colnames(yfreq))
  results <- cbind(xfreq,yfreq)
  return(results)
}