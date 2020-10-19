library("dplyr")
library("purrr")
#library(phyloseq)
library(ape)
library(ggplot2)
library(EnvStats)
library(miscTools)
library(devtools)


############################ GLOBAL ANALYSIS MAPPED TO DAOM ############################

#set output directories

setwd("/Users/wrobbins/Desktop/Phd/RADseq/common_sites/new_added_lines/")
dir <- "/Users/wrobbins/Desktop/Phd/RADseq/common_sites/new_added_lines/"
figs <- "/Users/wrobbins/Desktop/Phd/RADseq/common_sites/new_added_lines/"


## Read in binary data of RAD sites
RAD_files <- list.files(pattern = '*_binary.txt')
for (i in RAD_files) {
  x <- read.table(i, header = TRUE, sep = "")
  #x <- as.data.table(x)
  assign(i,x)
}

## Add file names to list
## B12_common_RAD_SNP_binary.txt not used because of 
RAD_dfs <- list(A1_common_RAD_SNP_binary.txt, A1ssp1_common_RAD_SNP_binary.txt, A1ssp2_common_RAD_SNP_binary.txt,A5_common_RAD_SNP_binary.txt,A5ssp1_common_RAD_SNP_binary.txt,A5ssp2_common_RAD_SNP_binary.txt,A5ssp3_common_RAD_SNP_binary.txt, B12ssp1_common_RAD_SNP_binary.txt, B12ssp3_common_RAD_SNP_binary.txt, B12ssp4_common_RAD_SNP_binary.txt, B12ssp5_common_RAD_SNP_binary.txt, C2_common_RAD_SNP_binary.txt,C2ssp1_common_RAD_SNP_binary.txt, C2ssp2_common_RAD_SNP_binary.txt, C2ssp3_common_RAD_SNP_binary.txt, C3_common_RAD_SNP_binary.txt, C3ssp3_common_RAD_SNP_binary.txt, C5_common_RAD_SNP_binary.txt, Sc1_common_RAD_SNP_binary.txt, Sc1a_common_RAD_SNP_binary.txt, Sc1b_common_RAD_SNP_binary.txt, Sc1c_common_RAD_SNP_binary.txt, Sc1d_common_RAD_SNP_binary.txt, Sc1e_common_RAD_SNP_binary.txt, Sc1f_common_RAD_SNP_binary.txt, Sc2_common_RAD_SNP_binary.txt, Sc2a_common_RAD_SNP_binary.txt, Sc2b_common_RAD_SNP_binary.txt,Sc2c_common_RAD_SNP_binary.txt, Sc2d_common_RAD_SNP_binary.txt, Sc2e_common_RAD_SNP_binary.txt, X21_common_RAD_SNP_binary.txt, X27_common_RAD_SNP_binary.txt, X31_common_RAD_SNP_binary.txt, X35_common_RAD_SNP_binary.txt, X36_common_RAD_SNP_binary.txt, X38_common_RAD_SNP_binary.txt, X39_common_RAD_SNP_binary.txt, Z12_common_RAD_SNP_binary.txt, Z15_common_RAD_SNP_binary.txt, Z16_common_RAD_SNP_binary.txt,Z17_common_RAD_SNP_binary.txt, Z18_common_RAD_SNP_binary.txt, Z19_common_RAD_SNP_binary.txt, Z25_common_RAD_SNP_binary.txt, Z31_common_RAD_SNP_binary.txt, Z34_common_RAD_SNP_binary.txt, Z34ssp2_common_RAD_SNP_binary.txt, Z41_common_RAD_SNP_binary.txt, Z43_common_RAD_SNP_binary.txt, Z64_common_RAD_SNP_binary.txt, Z65_common_RAD_SNP_binary.txt, Z68_common_RAD_SNP_binary.txt)

## Add a column of ROWNAMES as identifies
for (i in 1:length(RAD_dfs)){
  RAD_dfs[[i]]$ROWNAMES <- rownames(RAD_dfs[[i]])
}

## Load plyr and merge all sample dataframs keeping all data and rename the rownames
library("plyr")
all_RAD_sites <- join_all(RAD_dfs, by="ROWNAMES", type = "full")
rownames(all_RAD_sites) <- all_RAD_sites$ROWNAMES
all_RAD_sites$ROWNAMES <- NULL
all_RAD_sites <- all_RAD_sites[,which(colSums(all_RAD_sites,na.rm = TRUE) > 7500)]

## Load gtools and reorder the dataframe in ascending order of CHR_POS remove ROWNAMES column and replace all NA data with 0
library("gtools")
all_RAD_sites <- all_RAD_sites[mixedorder(rownames(all_RAD_sites)),]
all_RAD_sites[is.na(all_RAD_sites)] <- 0
all_RAD_sites_dataframe <- all_RAD_sites
## Transform the dataframe to a transposed matrix so it reads rowXcolumn = sampleXsite and compute binary distances with heirarchical clustering
all_RAD_sites <- t(as.matrix(all_RAD_sites))
all_RAD_dist <- dist(all_RAD_sites, method = "binary")
hclust_all_RAD_sites <- hclust(all_RAD_dist, method = "average")

## Load ape, dendextend and circlize to generate the rounded dendrogram
library("ape")
library("dendextend")
colors <- c("red","purple","blue")
line_groups <- cutree(hclust_all_RAD_sites, k = 3)
plot(as.phylo(hclust_all_RAD_sites), type="fan", tip.color = colors[line_groups], label.offset = 0.005, cex = 0.25)

RAD_dend <- all_RAD_sites %>% 
  dist(method = "binary") %>%
  hclust(method = "average")  %>%
  as.dendrogram %>%
  set("labels_col", k=3) %>% 
  set("branches_k_color", k=3) %>%
  set("labels_cex", 0.5) %>%
  set("leaves_pch", c(1)) %>% 
  set("leaves_cex", 1) %>%
  set("leaves_col", "black")

library("circlize")
circlize_dendrogram(RAD_dend)


##Calculate eigenvalues for the dataframe of sites PCA -> quantitative, CA -> tqualitative
library("tidyverse")
library("magrittr")
library("FactoMineR")
library("factoextra")
all_RAD_sites_mat <- t(all_RAD_sites_dataframe)
#all_RAD_dist <- dist(all_RAD_sites_mat,method = "binary")
RAD_sites_pca <- PCA(all_RAD_sites_mat, graph = FALSE)
#RAD_sites_ca <- CA(all_RAD_sites_mat)
RAD_sites_pca_eigen <- get_eigenvalue(RAD_sites_pca)
RAD_sites_pca_eigen

## Load sample data and generate figure
RAD_samples <- "RAD_sample_list.txt"
RAD_sampleIDs <- read.table(file=RAD_samples, sep = "", header = FALSE)
colnames(RAD_sampleIDs) <- c("SampleID", "Parent_Group", "New_Name", "Data_Source","Genetic_Group", "Nuclear_Genotype")
rownames(RAD_sampleIDs) <- RAD_sampleIDs$New_Name
sample_names <- rownames(all_RAD_sites_mat)
RAD_samples <- RAD_sampleIDs[which(rownames(RAD_sampleIDs) %in% sample_names),]


RAD_sampleIDs$Genetic_Group <- as.factor(RAD_sampleIDs$Genetic_Group)

RAD_sampleIDs$Nuclear_Genotype <- as.factor(RAD_sampleIDs$Nuclear_Genotype)

RAD_samples <- RAD_sampleIDs[which(RAD_sampleIDs$New_Name %in% all_RAD_sites_mat)]

fviz_pca_ind(RAD_sites_pca,
            geom.ind = "point",
            #select.ind = list(name = "Z17_R1", "Z17_R2", "Z17_R3", "Z17_R4"),
            pointsize = 4,
            pointshape = 19,
           # col.var = "black",
            #repel = TRUE,
            alpha.ind = 0.25,
            col.ind = as.factor(RAD_samples$Genetic_Group),
            palette = c("#0082CE", "#228B00", "#CC476B"),
           #"#00FF66", "#FF0000", "#FF9933", "#00CC33", "#33FF33", "#FF66FF", "#993333", "#00CC99", "#CC0000", "#0000FF"),
            #addEllipses = TRUE,
            legend.title = "Families") +
 theme(axis.line=element_blank(),
       panel.background=element_blank(),
       panel.border=element_blank(),
       #panel.grid.major=element_blank(),
       panel.grid.minor=element_blank(),
       plot.background=element_blank())


fviz_pca_var(RAD_sites_pca,
             #geom.ind = "point",
             #pointsize = 2,
             #pointshape = 19,
             #col.var = "black",
             select = list(cos2=15),
             repel = TRUE,
             #alpha.ind = 0.5,
             #col.ind = RAD_sampleIDs$Family_Group,
             #palette = c("#660066", "#33CC33", "#0066CC","#FF0000","#FF9900", "#CCFF33", "#999900", "#FF9933", "#3399FF", "#33FF33", "#00CC99", "#FF66FF"),
             #addEllipses = TRUE,
             legend.title = "Families") +
  
  theme(axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())


################# ISOLATE SPECIFIC ANALYSIS MAPPED TO INDIVIDUAL PARENT #################

## Reset directory locations

setwd("/Users/wrobbins/Desktop/Phd/RADseq/")
dir <- "/Users/wrobbins/Desktop/Phd/RADseq/"
figs <- "/Users/wrobbins/Desktop/Phd/RADseq/"

##Read in in silico RAD digest files
file1 <- "A1_insidi"
file2 <- "A5_insidi"
file3 <- "B12_insidi"
file4 <- "C2_insidi"
file5 <- "C3_insidi"
#file6 <- "C5_insidi"
A1 <- read.table(file=file1,  sep= "", header = TRUE)
A5 <- read.table(file=file2,  sep= "", header = TRUE)
B12 <- read.table(file=file3,  sep= "", header = TRUE)
C2 <- read.table(file=file4,  sep= "", header = TRUE)
C3 <- read.table(file=file5,  sep= "", header = TRUE)
#C5 <- read.table(file=file6,  sep= "", header = TRUE)

##Filter targets
inds <- which(A1$Enzyme_name == "EcoRI")
inds2 <- as.integer(inds-1)
ind <- sort(c(inds,inds2))
#rows <- lapply(inds, function(x) (x-1):(x))
A1_RAD_insidi <- A1[which(rownames(A1) %in% ind),]
A1_RAD_targets <- A1_RAD_insidi[which(A1_RAD_insidi$X3frag > 150 & A1_RAD_insidi$X3frag < 450),]
inds <- which(A5$Enzyme_name == "EcoRI")
inds2 <- as.integer(inds-1)
ind <- sort(c(inds,inds2))
#rows <- lapply(inds, function(x) (x-1):(x))
A5_RAD_insidi <- A5[which(rownames(A5) %in% ind),]
A5_RAD_targets <- A5_RAD_insidi[which(A5_RAD_insidi$X3frag > 150 & A5_RAD_insidi$X3frag < 450),]
inds <- which(B12$Enzyme_name == "EcoRI")
inds2 <- as.integer(inds-1)
ind <- sort(c(inds,inds2))
#rows <- lapply(inds, function(x) (x-1):(x))
B12_RAD_insidi <- B12[which(rownames(B12) %in% ind),]
B12_RAD_targets <- B12_RAD_insidi[which(B12_RAD_insidi$X3frag > 150 & B12_RAD_insidi$X3frag < 450),]
inds <- which(C2$Enzyme_name == "EcoRI")
inds2 <- as.integer(inds-1)
ind <- sort(c(inds,inds2))
#rows <- lapply(inds, function(x) (x-1):(x))
C2_RAD_insidi <- C2[which(rownames(C2) %in% ind),]
C2_RAD_targets <- C2_RAD_insidi[which(C2_RAD_insidi$X3frag > 150 & C2_RAD_insidi$X3frag < 450),]
inds <- which(C3$Enzyme_name == "EcoRI")
inds2 <- as.integer(inds-1)
ind <- sort(c(inds,inds2))
#rows <- lapply(inds, function(x) (x-1):(x))
C3_RAD_insidi <- C3[which(rownames(C3) %in% ind),]
C3_RAD_targets <- C3_RAD_insidi[which(C3_RAD_insidi$X3frag > 150 & C3_RAD_insidi$X3frag < 450),]

## Generage overlapping graph of all RAD targets
p <- ggplot(C3_RAD_targets, aes(x=X3frag)) +
  geom_density(alpha=0.15, fill="darkblue")+
  geom_density(data = A5_RAD_targets, aes(x=X3frag), fill="darkgreen", alpha=0.15)+
  geom_density(data = B12_RAD_targets, aes(x=X3frag), fill="yellow", alpha=0.15)+
  geom_density(data = C2_RAD_targets, aes(x=X3frag),fill="darkred", alpha=0.15)+
  geom_density(data = C3_RAD_targets, aes(x=X3frag), fill="purple", alpha=0.15)+
  scale_x_continuous(limits= c(0,650))+
  theme(panel.border= element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        # legend.position = "none",
        text = element_text(size=15),
        axis.line= element_line(colour="black"))

p

## Generate BED files for filtering and viewing in IGV
A1_RAD_insidi <- A1_RAD_insidi[,c(1:3,5)]
A1_sub <- "A1_RAD_insidi.bed"
write.table(A1_RAD_insidi, A1_sub, sep="\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
A5_RAD_insidi <- A5_RAD_insidi[,c(1:3,5)]
A5_sub <- "A5_RAD_insidi.bed"
write.table(A5_RAD_insidi, A5_sub, sep="\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
B12_RAD_insidi <- B12_RAD_insidi[,c(1:3,5)]
B12_sub <- "B12_RAD_insidi.bed"
write.table(B12_RAD_insidi, B12_sub, sep="\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
C2_RAD_insidi <- C2_RAD_insidi[,c(1:3,5)]
C2_sub <- "C2_RAD_insidi.bed"
write.table(C2_RAD_insidi, C2_sub, sep="\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
C3_RAD_insidi <- C3_RAD_insidi[,c(1:3,5)]
C3_sub <- "C3_RAD_insidi.bed"
write.table(C3_RAD_insidi, C3_sub, sep="\t", row.names = FALSE, col.names = FALSE,quote = FALSE)






################ ISOLATE SPECIFIC ALLELE FREQUENCIES COMPARED TO PARENT ################

###########################   A1   ############################
setwd("/Users/wrobbins/Desktop/AF/A1")
library("dplyr")
## Initiate a list of all files with allele frequency data and import to the environment
tmp <- list.files(pattern = "*.allfreq.txt")

dflist <- NULL
for (i in 1:length(tmp)) { 
  id <- strsplit(tmp[i], "\\." [[1]][1])
  assign(x = id[[1]][1], as.data.frame(read.table(tmp[i], sep = "\t", col.names = c("Sca", "Pos", "Ref", "Type",paste0(id[[1]][1],"_D"), paste0(id[[1]][1], "_RC"), paste0(id[[1]][1], "_AC"), paste0(id[[1]][1], "_C3")))[,1:7]))
}

dflist <- list(A1_R1=A1_R1, A1_R2 = A1_R2, A1_R3=A1_R3,A1_R4=A1_R4,
               A1ssp1_R1=A1ssp1_R1,A1ssp1_R2=A1ssp1_R2,A1ssp1_R3=A1ssp1_R3,A1ssp1_R4=A1ssp1_R4,A1ssp1_R5=A1ssp1_R5,A1ssp1_R6=A1ssp1_R6,
               A1ssp2_R1=A1ssp2_R1,A1ssp2_R2=A1ssp2_R2,A1ssp2_R3=A1ssp2_R3,A1ssp2_R4=A1ssp2_R4,A1ssp2_R5=A1ssp2_R5,A1ssp2_R6=A1ssp2_R6,
               Z25_R1=Z25_R1,Z25_R2=Z25_R2,Z25_R3=Z25_R3,Z25_R4=Z25_R4,Z25_R5=Z25_R5,Z25_R6=Z25_R6)
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

## For loop for all dataframes 
for(i in 1:length(dflist)){
  new <- new_df(dflist[[i]],noquote(names(dflist[i])))
  assign(names(dflist[i]),new)
}

## Sample specific lists
A1list <- list(A1_R1=A1_R1,A1_R2 = A1_R2, A1_R3=A1_R3,A1_R4=A1_R4)
A1ssp1list <- list(A1ssp1_R1=A1ssp1_R1,A1ssp1_R2=A1ssp1_R2,A1ssp1_R3=A1ssp1_R3,A1ssp1_R4=A1ssp1_R4,A1ssp1_R5=A1ssp1_R5,A1ssp1_R6=A1ssp1_R6)
A1ssp2list <- list(A1ssp2_R1=A1ssp2_R1,A1ssp2_R2=A1ssp2_R2,A1ssp2_R3=A1ssp2_R3,A1ssp2_R4=A1ssp2_R4,A1ssp2_R5=A1ssp2_R5,A1ssp2_R6=A1ssp2_R6)
Z25list <- list(Z25_R1=Z25_R1,Z25_R2=Z25_R2,Z25_R3=Z25_R3,Z25_R4=Z25_R4,Z25_R5=Z25_R5,Z25_R6=Z25_R6)

###########################   A5   ############################

setwd("/Users/wrobbins/Desktop/AF/A5")

## Initiate a list of all files with allele frequency data and import to the environment
tmp <- list.files(pattern = "*.allfreq.txt")

dflist <- NULL
for (i in 1:length(tmp)) { 
  id <- strsplit(tmp[i], "\\." [[1]][1])
  assign(x = id[[1]][1], as.data.frame(read.table(tmp[i], sep = "\t", col.names = c("Sca", "Pos", "Ref", "Type",paste0(id[[1]][1],"_D"), paste0(id[[1]][1], "_RC"), paste0(id[[1]][1], "_AC"), paste0(id[[1]][1], "_C3")))[,1:7]))
}

dflist <- list(A5_R1 = A5_R1, A5_R2 = A5_R2, A5_R3=A5_R3,A5_R4=A5_R4,A5_R5=A5_R5,
               A5ssp1_R1=A5ssp1_R1,A5ssp1_R2=A5ssp1_R2,A5ssp1_R3=A5ssp1_R3,A5ssp1_R4=A5ssp1_R4,A5ssp1_R5=A5ssp1_R5,
               A5ssp2_R1=A5ssp2_R1,A5ssp2_R2=A5ssp2_R2,A5ssp2_R3=A5ssp2_R3,A5ssp2_R4=A5ssp2_R4,A5ssp2_R5=A5ssp2_R5,A5ssp2_R6=A5ssp2_R6,
               A5ssp3_R1=A5ssp3_R1,A5ssp3_R2=A5ssp3_R2,A5ssp3_R3=A5ssp3_R3,A5ssp3_R4=A5ssp3_R4,A5ssp3_R5=A5ssp3_R5,A5ssp3_R6=A5ssp3_R6,
               X21_R1=X21_R1,X21_R2=X21_R2,X21_R3=X21_R3,X21_R4=X21_R4,X21_R5=X21_R5,
               X27_R1=X27_R1,X27_R2=X27_R2,X27_R3=X27_R3,X27_R4=X27_R4,X27_R5=X27_R5,X27_R6=X27_R6,
               Z64_R1=Z64_R1,Z64_R2=Z64_R2,Z64_R3=Z64_R3,Z64_R4=Z64_R4,Z64_R5=Z64_R5,
               Z65_R1=Z65_R1,Z65_R2=Z65_R2,Z65_R3=Z65_R3,Z65_R4=Z65_R4,Z65_R5=Z65_R5,Z65_R6=Z65_R6,
               Z68_R1=Z68_R1,Z68_R2=Z68_R2,Z68_R3=Z68_R3,Z68_R4=Z68_R4,Z68_R5=Z68_R5)


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

## For loop for all dataframes 
for(i in 1:length(dflist)){
  new <- new_df(dflist[[i]],noquote(names(dflist[i])))
  assign(names(dflist[i]),new)
}

## Sample specific lists
A5list <- list(A5_R1 = A5_R1, A5_R2 = A5_R2, A5_R3=A5_R3,A5_R4=A5_R4,A5_R5=A5_R5)
A5ssp1list <- list(A5ssp1_R1=A5ssp1_R1,A5ssp1_R2=A5ssp1_R2,A5ssp1_R3=A5ssp1_R3,A5ssp1_R4=A5ssp1_R4,A5ssp1_R5=A5ssp1_R5)
A5ssp2list <- list(A5ssp2_R1=A5ssp2_R1,A5ssp2_R2=A5ssp2_R2,A5ssp2_R3=A5ssp2_R3,A5ssp2_R4=A5ssp2_R4,A5ssp2_R5=A5ssp2_R5,A5ssp2_R6=A5ssp2_R6)
A5ssp3list <- list(A5ssp3_R1=A5ssp3_R1,A5ssp3_R2=A5ssp3_R2,A5ssp3_R3=A5ssp3_R3,A5ssp3_R4=A5ssp3_R4,A5ssp3_R5=A5ssp3_R5,A5ssp3_R6=A5ssp3_R6)
X21list <- list(X21_R1=X21_R1,X21_R2=X21_R2,X21_R3=X21_R3,X21_R4=X21_R4,X21_R5=X21_R5)
X27list <- list(X27_R1=X27_R1,X27_R2=X27_R2,X27_R3=X27_R3,X27_R4=X27_R4,X27_R5=X27_R5,X27_R6=X27_R6)
Z64list <- list(Z64_R1=Z64_R1,Z64_R2=Z64_R2,Z64_R3=Z64_R3,Z64_R4=Z64_R4,Z64_R5=Z64_R5)
Z65list <- list(Z65_R1=Z65_R1,Z65_R2=Z65_R2,Z65_R3=Z65_R3,Z65_R4=Z65_R4,Z65_R5=Z65_R5,Z65_R6=Z65_R6)
Z68list <- list(Z68_R1=Z68_R1,Z68_R2=Z68_R2,Z68_R3=Z68_R3,Z68_R4=Z68_R4,Z68_R5=Z68_R5)

###########################   B12   ############################
########### B12 parent samples came from Tania Wyss ISME paper
##B12-E1 --> B12_R1
##B12-I --> B12_R2
##B12-J --> B12_R3
setwd("/Users/wrobbins/Desktop/AF/B12")

## Initiate a list of all files with allele frequency data and import to the environment
tmp <- list.files(pattern = "*.allfreq.txt")

dflist <- NULL
for (i in 1:length(tmp)) { 
  id <- strsplit(tmp[i], "\\." [[1]][1])
  assign(x = id[[1]][1], as.data.frame(read.table(tmp[i], sep = "\t", col.names = c("Sca", "Pos", "Ref", "Type",paste0(id[[1]][1],"_D"), paste0(id[[1]][1], "_RC"), paste0(id[[1]][1], "_AC"), paste0(id[[1]][1], "_C3")))[,1:7]))
}

dflist <- list(B12_R1 = B12_R1, B12_R2 = B12_R2, B12_R3 = B12_R3,
               B12ssp1_R1 = B12ssp1_R1, B12ssp1_R2 = B12ssp1_R2, B12ssp1_R3=B12ssp1_R3,B12ssp1_R4=B12ssp1_R4,B12ssp1_R5=B12ssp1_R5,
               B12ssp3_R1=B12ssp3_R1,B12ssp3_R2=B12ssp3_R2,B12ssp3_R3=B12ssp3_R3,B12ssp3_R4=B12ssp3_R4,B12ssp3_R5=B12ssp3_R5,B12ssp3_R6=B12ssp3_R6,
               B12ssp4_R1=B12ssp4_R1,B12ssp4_R2=B12ssp4_R2,B12ssp4_R3=B12ssp4_R3,B12ssp4_R4=B12ssp4_R4,B12ssp4_R5=B12ssp4_R5,
               B12ssp5_R1=B12ssp5_R1,B12ssp5_R2=B12ssp5_R2,B12ssp5_R3=B12ssp5_R3,B12ssp5_R4=B12ssp5_R4,B12ssp5_R5=B12ssp5_R5,
               Z31_R1=Z31_R1,Z31_R2=Z31_R2,Z31_R3=Z31_R3,Z31_R4=Z31_R4,Z31_R5=Z31_R5,Z31_R6=Z31_R6,
               Z34_R1=Z34_R1,Z34_R2=Z34_R2,Z34_R3=Z34_R3,Z34_R4=Z34_R4,Z34_R5=Z34_R5,Z34_R6=Z34_R6,
               Z34ssp2_R1=Z34ssp2_R1,Z34ssp2_R2=Z34ssp2_R2,Z34ssp2_R3=Z34ssp2_R3)

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

## For loop for all dataframes 
for(i in 1:length(dflist)){
  new <- new_df(dflist[[i]],noquote(names(dflist[i])))
  assign(names(dflist[i]),new)
}

## Sample specific lists
B12list <- list(B12_R1 = B12_R1, B12_R2 = B12_R2, B12_R3 = B12_R3)
B12ssp1list <- list(B12ssp1_R1 = B12ssp1_R1, B12ssp1_R2 = B12ssp1_R2, B12ssp1_R3=B12ssp1_R3,B12ssp1_R4=B12ssp1_R4,B12ssp1_R5=B12ssp1_R5)
B12ssp3list <- list(B12ssp3_R1=B12ssp3_R1,B12ssp3_R2=B12ssp3_R2,B12ssp3_R3=B12ssp3_R3,B12ssp3_R4=B12ssp3_R4,B12ssp3_R5=B12ssp3_R5,B12ssp3_R6=B12ssp3_R6)
B12ssp4list <- list(B12ssp4_R1=B12ssp4_R1,B12ssp4_R2=B12ssp4_R2,B12ssp4_R3=B12ssp4_R3,B12ssp4_R4=B12ssp4_R4,B12ssp4_R5=B12ssp4_R5)
B12ssp5list <- list(B12ssp5_R1=B12ssp5_R1,B12ssp5_R2=B12ssp5_R2,B12ssp5_R3=B12ssp5_R3,B12ssp5_R4=B12ssp5_R4,B12ssp5_R5=B12ssp5_R5)
Z31list <- list(Z31_R1=Z31_R1,Z31_R2=Z31_R2,Z31_R3=Z31_R3,Z31_R4=Z31_R4,Z31_R5=Z31_R5,Z31_R6=Z31_R6)
Z34list <- list(Z34_R1=Z34_R1,Z34_R2=Z34_R2,Z34_R3=Z34_R3,Z34_R4=Z34_R4,Z34_R5=Z34_R5,Z34_R6=Z34_R6)
Z34ssp2list <- list(Z34ssp2_R1=Z34ssp2_R1,Z34ssp2_R2=Z34ssp2_R2,Z34ssp2_R3=Z34ssp2_R3)

###########################   C2/C5   ############################

setwd("/Users/wrobbins/Desktop/AF/C2")

## Initiate a list of all files with allele frequency data and import to the environment
tmp <- list.files(pattern = "*.allfreq.txt")

dflist <- NULL
for (i in 1:length(tmp)) { 
  id <- strsplit(tmp[i], "\\." [[1]][1])
  assign(x = id[[1]][1], as.data.frame(read.table(tmp[i], sep = "\t", col.names = c("Sca", "Pos", "Ref", "Type",paste0(id[[1]][1],"_D"), paste0(id[[1]][1], "_RC"), paste0(id[[1]][1], "_AC"), paste0(id[[1]][1], "_C3")))[,1:7]))
}

dflist <- list(#C5_R1 = C5_R1, C5_R2 = C5_R2, C5_R3 = C5_R3, C5_R4 = C5_R4, C5_R5 = C5_R5, C5_R6 = C5_R6,
               C2_R1 = C2_R1, C2_R2 = C2_R2, C2_R3=C2_R3,C2_R4=C2_R4,C2_R5=C2_R5,C2_R6=C2_R6,
               C2ssp1_R1=C2ssp1_R1,C2ssp1_R2=C2ssp1_R2,C2ssp1_R3=C2ssp1_R3,C2ssp1_R4=C2ssp1_R4,C2ssp1_R5=C2ssp1_R5,
               C2ssp2_R1=C2ssp2_R1,C2ssp2_R2=C2ssp2_R2,C2ssp2_R3=C2ssp2_R3,C2ssp2_R4=C2ssp2_R4,C2ssp2_R5=C2ssp2_R5,C2ssp2_R6=C2ssp2_R6,
               C2ssp3_R1=C2ssp3_R1,C2ssp3_R2=C2ssp3_R2,C2ssp3_R3=C2ssp3_R3,C2ssp3_R4=C2ssp3_R4,C2ssp3_R5=C2ssp3_R5,C2ssp3_R6=C2ssp3_R6,
               X35_R1=X35_R1,X35_R2=X35_R2,X35_R3=X35_R3,X35_R4=X35_R4,X35_R5=X35_R5,
               X36_R1=X36_R1,X36_R2=X36_R2,X36_R3=X36_R3,X36_R4=X36_R4,X36_R5=X36_R5,X36_R6=X36_R6,
               X38_R1=X38_R1,X38_R2=X38_R2,X38_R3=X38_R3,X38_R4=X38_R4,X38_R5=X38_R5,X38_R6=X38_R6,
               X39_R1=X39_R1,X39_R2=X39_R2,X39_R3=X39_R3,X39_R4=X39_R4)#,X39_R5=X39_R5,X39_R6=X39_R6,

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

## For loop for all dataframes 
for(i in 1:length(dflist)){
  new <- new_df(dflist[[i]],noquote(names(dflist[i])))
  assign(names(dflist[i]),new)
}

## Sample specific lists
#C5list <- list (C5_R1=C5_R1, C5_R2=C5_R2, C5_R3=C5_R3, C5_R4=C5_R4, C5_R5=C5_R5, C5_R6=C5_R6)
C2list <- list(C2_R1=C2_R1, C2_R2=C2_R2, C2_R3=C2_R3, C2_R4=C2_R4, C2_R5=C2_R5, C2_R6=C2_R6)
C2ssp1list <- list(C2ssp1_R1=C2ssp1_R1, C2ssp1_R2=C2ssp1_R2, C2ssp1_R3=C2ssp1_R3, C2ssp1_R4=C2ssp1_R4, C2ssp1_R5=C2ssp1_R5)
C2ssp2list <- list(C2ssp2_R1=C2ssp2_R1, C2ssp2_R2=C2ssp2_R2, C2ssp2_R3=C2ssp2_R3, C2ssp2_R4=C2ssp2_R4, C2ssp2_R5=C2ssp2_R5, C2ssp2_R6=C2ssp2_R6)
C2ssp3list <- list(C2ssp3_R1=C2ssp3_R1, C2ssp3_R2=C2ssp3_R2, C2ssp3_R3=C2ssp3_R3, C2ssp3_R4=C2ssp3_R4, C2ssp3_R5=C2ssp3_R5, C2ssp3_R6=C2ssp3_R6)
X35list <- list(X35_R1=X35_R1, X35_R2=X35_R2, X35_R3=X35_R3, X35_R4=X35_R4, X35_R5=X35_R5)
X36list <- list(X36_R1=X36_R1, X36_R2=X36_R2, X36_R3=X36_R3, X36_R4=X36_R4, X36_R5=X36_R5, X36_R6=X36_R6)
X38list <- list(X38_R1=X38_R1, X38_R2=X38_R2, X38_R3=X38_R3, X38_R4=X38_R4, X38_R5=X38_R5, X38_R6=X38_R6)
X39list <- list(X39_R1=X39_R1, X39_R2=X39_R2, X39_R3=X39_R3, X39_R4=X39_R4)#, X39_R5, X39_R6

###########################   C3   ############################
########### C3 parent samples came from Tania Wyss and Pawel Rosikiewicz
##Sc2-1 --> Sc2_R6
##Sc2-2 --> Sc2_R7
##Sc2-3 --> Sc2_R8
##Sc2c-1 --> Sc2c_R4
##Sc2c-2 --> Sc2c_R2
##Sc2c-3 --> Sc2c_R3
##Sc1-1 --> Sc1_R1
##Sc1-3 --> Sc1_R3
##Sc1-4 --> Sc1_R4
setwd("/Users/wrobbins/Desktop/AF/C3")

## Initiate a list of all files with allele frequency data and import to the environment
tmp <- list.files(pattern = "*.allfreq.txt")

dflist <- NULL
for (i in 1:length(tmp)) { 
  id <- strsplit(tmp[i], "\\." [[1]][1])
  assign(x = id[[1]][1], as.data.frame(read.table(tmp[i], sep = "\t", col.names = c("Sca", "Pos", "Ref", "Type",paste0(id[[1]][1],"_D"), paste0(id[[1]][1], "_RC"), paste0(id[[1]][1], "_AC"), paste0(id[[1]][1], "_C3")))[,1:7]))
}

dflist <- list(C3_R1=C3_R1,C3_R2=C3_R2,C3_R3=C3_R3,C3_R4=C3_R4,C3_R5=C3_R5,
               C3ssp3_R1=C3ssp3_R1,C3ssp3_R2=C3ssp3_R2,C3ssp3_R3=C3ssp3_R3,C3ssp3_R4=C3ssp3_R4,C3ssp3_R5=C3ssp3_R5,C3ssp3_R6=C3ssp3_R6,
               Sc1_R1=Sc1_R1,Sc1_R3=Sc1_R3,Sc1_R4=Sc1_R4,
               Sc1a_R1=Sc1a_R1,Sc1a_R2=Sc1a_R2,Sc1a_R3=Sc1a_R3,Sc1a_R4=Sc1a_R4,Sc1a_R5=Sc1a_R5,
               Sc1b_R1=Sc1b_R1,Sc1b_R2=Sc1b_R2,Sc1b_R3=Sc1b_R3,Sc1b_R4=Sc1b_R4,Sc1b_R5=Sc1b_R5,
               Sc1c_R1=Sc1c_R1,Sc1c_R2=Sc1c_R2,Sc1c_R3=Sc1c_R3,Sc1c_R4=Sc1c_R4,
               Sc1d_R1=Sc1d_R1,Sc1d_R2=Sc1d_R2,Sc1d_R3=Sc2d_R3,Sc1d_R4=Sc1d_R4,Sc1d_R5=Sc1d_R5,
               Sc1e_R1=Sc1e_R1,Sc1e_R2=Sc1e_R2,Sc1e_R3=Sc1e_R3,Sc1e_R4=Sc1e_R4,Sc1e_R5=Sc1e_R5, Sc1e_R6=Sc1e_R6,
               Sc1f_R1=Sc1f_R1,Sc1f_R2=Sc1f_R2,Sc1f_R3=Sc1f_R3,Sc1f_R4=Sc1f_R4,Sc1f_R5=Sc1f_R5,
               Sc2_R1=Sc2_R1,Sc2_R2=Sc2_R2,Sc2_R3=Sc2_R3,Sc2_R4=Sc2_R4,Sc2_R5=Sc2_R5,Sc2_R6=Sc2_R6,Sc2_R7=Sc2_R7,Sc2_R8=Sc2_R8,
               Sc2a_R1=Sc2a_R1,Sc2a_R2=Sc2a_R2,Sc2a_R3=Sc2a_R3,Sc2a_R4=Sc2a_R4,
               Sc2c_R1=Sc2c_R1,Sc2c_R2=Sc2c_R2,Sc2c_R3=Sc2c_R3,Sc2c_R4=Sc2c_R4,
               Sc2b_R1=Sc2b_R1,Sc2b_R2=Sc2b_R2,Sc2b_R3=Sc2b_R3,Sc2b_R4=Sc2b_R4,Sc2b_R5=Sc2b_R5,
               Sc2d_R1=Sc2d_R1,Sc2d_R2=Sc2d_R2,Sc2d_R3=Sc2d_R3,Sc2d_R4=Sc2d_R4,Sc2d_R5=Sc2d_R5,
               Sc2e_R1=Sc2e_R1,Sc2e_R2=Sc2e_R2,Sc2e_R3=Sc2e_R3,Sc2e_R4=Sc2e_R4,
               X31_R1=X31_R1,X31_R2=X31_R2,X31_R3=X31_R3,X31_R4=X31_R4,X31_R5=X31_R5,X31_R6=X31_R6,
               Z12_R1=Z12_R1,Z12_R2=Z12_R2,Z12_R3=Z12_R3,Z12_R4=Z12_R4,Z12_R5=Z12_R5,
               Z15_R1=Z15_R1,Z15_R2=Z15_R2,Z15_R3=Z15_R3,Z15_R4=Z15_R4,Z15_R5=Z15_R5,Z15_R6=Z15_R6,
               Z16_R1=Z16_R1,Z16_R2=Z16_R2,Z16_R3=Z16_R3,Z16_R4=Z16_R4,Z16_R5=Z16_R5,Z16_R6=Z16_R6,
               Z17_R1=Z17_R1,Z17_R2=Z17_R2,Z17_R3=Z17_R3,Z17_R4=Z17_R4,
               Z18_R1=Z18_R1,Z18_R2=Z18_R2,Z18_R3=Z18_R3,Z18_R4=Z18_R4,Z18_R5=Z18_R5,Z18_R6=Z18_R6,
               Z19_R1=Z19_R1,Z19_R2=Z19_R2,Z19_R3=Z19_R3,Z19_R4=Z19_R4,Z19_R5=Z19_R5,Z19_R6=Z19_R6,
               Z41_R1=Z41_R1,Z41_R2=Z41_R2,Z41_R3=Z41_R3,Z41_R4=Z41_R4,Z41_R5=Z41_R5,
               Z43_R1=Z43_R1,Z43_R2=Z43_R2,Z43_R3=Z43_R3,Z43_R4=Z43_R4,Z43_R5=Z43_R5,Z43_R6=Z43_R6)


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

## For loop for all dataframes 
for(i in 1:length(dflist)){
  new <- new_df(dflist[[i]],noquote(names(dflist[i])))
  assign(names(dflist[i]),new)
}

## Sample specific lists
C3list <- list(C3_R1=C3_R1, C3_R2=C3_R2, C3_R3=C3_R3,C3_R4=C3_R4,C3_R5=C3_R5)
C3ssp3list <- list(C3ssp3_R1=C3ssp3_R1,C3ssp3_R2=C3ssp3_R2,C3ssp3_R3=C3ssp3_R3,C3ssp3_R4=C3ssp3_R4,C3ssp3_R5=C3ssp3_R5,C3ssp3_R6=C3ssp3_R6)
Sc1list <- list(Sc1_R1=Sc1_R1,Sc1_R3=Sc1_R3,Sc1_R4=Sc1_R4)
Sc1alist <- list(Sc1a_R1=Sc1a_R1,Sc1a_R2=Sc1a_R2,Sc1a_R3=Sc1a_R3,Sc1a_R4=Sc1a_R4,Sc1a_R5=Sc1a_R5)
Sc1blist <- list(Sc1b_R1=Sc1b_R1,Sc1b_R2=Sc1b_R2,Sc1b_R3=Sc1b_R3,Sc1b_R4=Sc1b_R4,Sc1b_R5=Sc1b_R5)
Sc1clist <- list(Sc1c_R1=Sc1c_R1,Sc1c_R2=Sc1c_R2,Sc1c_R3=Sc1c_R3,Sc1c_R4=Sc1c_R4)
Sc1dlist <- list(Sc1d_R1=Sc1d_R1,Sc1d_R2=Sc1d_R2,Sc1d_R3=Sc1d_R3,Sc1d_R4=Sc1d_R4,Sc1d_R5=Sc1d_R5)
Sc1elist <- list(Sc1e_R1=Sc1e_R1,Sc1e_R2=Sc1e_R2,Sc1e_R3=Sc1e_R3,Sc1e_R4=Sc1e_R4,Sc1e_R5=Sc1e_R5, Sc1e_R6=Sc1e_R6)
Sc1flist <- list(Sc1f_R1=Sc1f_R1,Sc1f_R2=Sc1f_R2,Sc1f_R3=Sc1f_R3,Sc1f_R4=Sc1f_R4,Sc1f_R5=Sc1f_R5)
Sc2list <- list(Sc2_R1=Sc2_R1,Sc2_R2=Sc2_R2,Sc2_R3=Sc2_R3,Sc2_R4=Sc2_R4,Sc2_R5=Sc2_R5,Sc2_R6=Sc2_R6,Sc2_R7=Sc2_R7,Sc2_R8=Sc2_R8)
Sc2clist <- list(Sc2c_R1=Sc2c_R1,Sc2c_R2=Sc2c_R2,Sc2c_R3=Sc2c_R3,Sc2c_R4=Sc2c_R4)
Sc2alist <- list(Sc2a_R1=Sc2a_R1,Sc2a_R2=Sc2a_R2,Sc2a_R3=Sc2a_R3,Sc2a_R4=Sc2a_R4)
Sc2blist <- list(Sc2b_R1=Sc2b_R1,Sc2b_R2=Sc2b_R2,Sc2b_R3=Sc2b_R3,Sc2b_R4=Sc2b_R4,Sc2b_R5=Sc2b_R5)
Sc2dlist <- list(Sc2d_R1=Sc2d_R1,Sc2d_R2=Sc2d_R2,Sc2d_R3=Sc2d_R3,Sc2d_R4=Sc2d_R4,Sc2d_R5=Sc2d_R5)
Sc2elist <- list(Sc2e_R1=Sc2e_R1,Sc2e_R2=Sc2e_R2,Sc2e_R3=Sc2e_R3,Sc2e_R4=Sc2e_R4)
X31list <- list(X31_R1=X31_R1,X31_R2=X31_R2,X31_R3=X31_R3,X31_R4=X31_R4,X31_R5=X31_R5,X31_R6=X31_R6)
Z12list <- list(Z12_R1=Z12_R1,Z12_R2=Z12_R2,Z12_R3=Z12_R3,Z12_R4=Z12_R4,Z12_R5=Z12_R5)
Z15list <- list(Z15_R1=Z15_R1,Z15_R2=Z15_R2,Z15_R3=Z15_R3,Z15_R4=Z15_R4,Z15_R5=Z15_R5,Z15_R6=Z15_R6)
Z16list <- list(Z16_R1=Z16_R1,Z16_R2=Z16_R2,Z16_R3=Z16_R3,Z16_R4=Z16_R4,Z16_R5=Z16_R5,Z16_R6=Z16_R6)
Z17list <- list(Z17_R1=Z17_R1,Z17_R2=Z17_R2,Z17_R3=Z17_R3,Z17_R4=Z17_R4)
Z18list <- list(Z18_R1=Z18_R1,Z18_R2=Z18_R2,Z18_R3=Z18_R3,Z18_R4=Z18_R4,Z18_R5=Z18_R5,Z18_R6=Z18_R6)
Z19list <- list(Z19_R1=Z19_R1,Z19_R2=Z19_R2,Z19_R3=Z19_R3,Z19_R4=Z19_R4,Z19_R5=Z19_R5,Z19_R6=Z19_R6)
Z41list <- list(Z41_R1=Z41_R1,Z41_R2=Z41_R2,Z41_R3=Z41_R3,Z41_R4=Z41_R4,Z41_R5=Z41_R5)
Z43list <- list(Z43_R1=Z43_R1,Z43_R2=Z43_R2,Z43_R3=Z43_R3,Z43_R4=Z43_R4,Z43_R5=Z43_R5,Z43_R6=Z43_R6)


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

###########################   A1   ############################
## Generate combined dataframes with all RAD sites from all samples
A1_fam <- combine_sites(A1list) 
A1ssp1_fam <- combine_sites(A1ssp1list)
A1ssp2_fam <- combine_sites(A1ssp2list)
Z25_fam <- combine_sites(Z25list) 

###########################   A5   ############################
## Generate combined dataframes with all RAD sites from all samples
A5_fam <- combine_sites(A5list)
A5ssp1_fam <- combine_sites(A5ssp1list) 
A5ssp2_fam <- combine_sites(A5ssp2list) 
A5ssp3_fam <- combine_sites(A5ssp3list)
X21_fam <- combine_sites(X21list)
X27_fam <- combine_sites(X27list)
Z64_fam <- combine_sites(Z64list)
Z65_fam <- combine_sites(Z65list)
Z68_fam <- combine_sites(Z68list) 

###########################   B12   ############################
## Generate combined dataframes with all RAD sites from all samples
B12_fam <- combine_sites(B12list)
B12ssp1_fam <- combine_sites(B12ssp1list)
B12ssp3_fam <- combine_sites(B12ssp3list)
B12ssp4_fam <- combine_sites(B12ssp4list)
B12ssp5_fam <- combine_sites(B12ssp5list)
Z31_fam <- combine_sites(Z31list)
Z34_fam <- combine_sites(Z34list)
Z34ssp2_fam <- combine_sites(Z34ssp2list) 

###########################   C2   ############################
## Generate combined dataframes with all RAD sites from all samples
#C5_fam <- combine_sites(C5list)
C2_fam <- combine_sites(C2list) 
C2ssp1_fam <- combine_sites(C2ssp1list)
C2ssp2_fam <- combine_sites(C2ssp2list)
C2ssp3_fam <- combine_sites(C2ssp3list)
X35_fam <- combine_sites(X35list)
X36_fam <- combine_sites(X36list)
X38_fam <- combine_sites(X38list)
X39_fam <- combine_sites(X39list)

###########################   C3   ############################
## Generate combined dataframes with all RAD sites from all samples
C3_fam <- combine_sites(C3list)
C3ssp3_fam <- combine_sites(C3ssp3list)
Sc1_fam <- combine_sites(Sc1list)
Sc1a_fam <- combine_sites(Sc1alist)
Sc1b_fam <- combine_sites(Sc1blist)
Sc1c_fam <- combine_sites(Sc1clist)
Sc1d_fam <- combine_sites(Sc1dlist)
Sc1e_fam <- combine_sites(Sc1elist)
Sc1f_fam <- combine_sites(Sc1flist)
Sc2_fam <- combine_sites(Sc2list)
Sc2a_fam <- combine_sites(Sc2alist)
Sc2b_fam <- combine_sites(Sc2blist)
Sc2c_fam <- combine_sites(Sc2clist)
Sc2d_fam <- combine_sites(Sc2dlist)
Sc2e_fam <- combine_sites(Sc2elist)
X31_fam <- combine_sites(X31list)
Z12_fam <- combine_sites(Z12list)
Z15_fam <- combine_sites(Z15list)
Z16_fam <- combine_sites(Z16list)
Z17_fam <- combine_sites(Z17list)
Z18_fam <- combine_sites(Z18list)
Z19_fam <- combine_sites(Z19list)
Z41_fam <- combine_sites(Z41list)
Z43_fam <- combine_sites(Z43list)

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


###########################   A1   ############################
## Combine family groups list
famlist <- list(A1_fam =A1_fam,A1ssp1_fam = A1ssp1_fam, A1ssp2_fam =A1ssp2_fam, Z25_fam =Z25_fam) 

## For loop to filter sites in all dataframes 
for(i in 1:length(famlist)){
  new <- filter_df(famlist[[i]],noquote(names(famlist[i])))
  assign(names(famlist[i]),new)
}

## Find the intersect of common RAD loci by family groups
common_A1_sites <- Reduce(intersect, list(rownames(A1_fam), rownames(A1ssp1_fam), rownames(A1ssp2_fam), rownames(Z25_fam)))


###########################   A5   ############################
## Combine family groups list
famlist <- list(A5_fam=A5_fam,A5ssp1_fam=A5ssp1_fam, A5ssp2_fam=A5ssp2_fam,A5ssp3_fam=A5ssp3_fam,X21_fam=X21_fam,X27_fam=X27_fam,Z64_fam=Z64_fam,Z65_fam=Z65_fam, Z68_fam =Z68_fam)

## For loop to filter sites in all dataframes 
for(i in 1:length(famlist)){
  new <- filter_df(famlist[[i]],noquote(names(famlist[i])))
  assign(names(famlist[i]),new)
}

## Find the intersect of common RAD loci by family groups
common_A5_sites <- Reduce(intersect, list(rownames(A5_fam), rownames(A5ssp1_fam), rownames(A5ssp3_fam), rownames(A5ssp2_fam), rownames(X21_fam), rownames(X27_fam), rownames(Z64_fam), rownames(Z65_fam), rownames(Z68_fam)))


###########################   B12   ############################
## Combine family groups list
famlist <- list(B12_fam=B12_fam,B12ssp1_fam = B12ssp1_fam, B12ssp3_fam=B12ssp3_fam, B12ssp4_fam = B12ssp4_fam, B12ssp5_fam=B12ssp5_fam,Z31_fam =Z31_fam,Z34_fam =Z34_fam, Z34ssp2_fam =Z34ssp2_fam)

## For loop to filter sites in all dataframes 
for(i in 1:length(famlist)){
  new <- filter_df(famlist[[i]],noquote(names(famlist[i])))
  assign(names(famlist[i]),new)
}

## Find the intersect of common RAD loci by family groups
common_B12_sites <- Reduce(intersect, list(rownames(B12_fam),rownames(B12ssp1_fam), rownames(B12ssp3_fam), rownames(B12ssp4_fam), rownames(B12ssp5_fam), rownames(Z31_fam), rownames(Z34_fam), rownames(Z34ssp2_fam)))

###########################   C2   ############################
## Combine family groups list
#C5_fam =C5_fam,
famlist <- list( C2_fam =C2_fam,C2ssp1_fam =C2ssp1_fam, C2ssp2_fam =C2ssp2_fam, C2ssp3_fam=C2ssp3_fam, X35_fam =X35_fam, X36_fam =X36_fam, X38_fam=X38_fam, X39_fam= X39_fam)

## For loop to filter sites in all dataframes 
for(i in 1:length(famlist)){
  new <- filter_df(famlist[[i]],noquote(names(famlist[i])))
  assign(names(famlist[i]),new)
}

## Find the intersect of common RAD loci by family groups
#rownames(C5_fam),
common_C2_sites <- Reduce(intersect, list( rownames(C2_fam), rownames(C2ssp1_fam), rownames(C2ssp2_fam), rownames(C2ssp3_fam), rownames(X35_fam), rownames(X36_fam), rownames(X38_fam), rownames(X39_fam)))

###########################   C3   ############################
## Combine family groups list
famlist <- list(C3_fam=C3_fam, C3ssp3_fam =C3ssp3_fam, Sc1_fam=Sc1_fam, Sc1a_fam =Sc1a_fam, Sc1b_fam =Sc1b_fam, Sc1c_fam =Sc1c_fam, Sc1d_fam=Sc1d_fam, Sc1e_fam =Sc1e_fam, Sc1f_fam= Sc1f_fam,Sc2_fam=Sc2_fam, Sc2a_fam =Sc2a_fam, Sc2b_fam =Sc2b_fam, Sc2c_fam=Sc2c_fam, Sc2d_fam =Sc2d_fam, Sc2e_fam =Sc2e_fam, X31_fam = X31_fam, Z12_fam =Z12_fam, Z15_fam =Z15_fam,Z16_fam =Z16_fam, Z17_fam =Z17_fam, Z18_fam =Z18_fam, Z19_fam =Z19_fam,Z41_fam =Z41_fam, Z43_fam =Z43_fam)


## For loop to filter sites in all dataframes 
for(i in 1:length(famlist)){
  new <- filter_df(famlist[[i]],noquote(names(famlist[i])))
  assign(names(famlist[i]),new)
}

## Find the intersect of common RAD loci by family groups
common_C3_sites <- Reduce(intersect, list(rownames(C3_fam), rownames(C3ssp3_fam),  rownames(Sc1_fam), rownames(Sc1a_fam), rownames(Sc1b_fam), rownames(Sc1c_fam), rownames(Sc1d_fam), rownames(Sc1e_fam), rownames(Sc1f_fam),rownames(Sc2_fam), rownames(Sc2a_fam), rownames(Sc2b_fam), rownames(Sc2c_fam), rownames(Sc2d_fam), rownames(Sc2e_fam), rownames(X31_fam), rownames(Z12_fam),  rownames(Z15_fam), rownames(Z16_fam), rownames(Z17_fam), rownames(Z18_fam), rownames(Z19_fam),rownames(Z41_fam), rownames(Z43_fam)))


###########################   A1   ############################
## Refresh family groups list
famlist <- list(A1_fam =A1_fam,A1ssp1_fam = A1ssp1_fam, A1ssp2_fam =A1ssp2_fam, Z25_fam =Z25_fam) 

## Filter all family groups for commonly occuring RAD SNP loci ####CHANGE sites
for(i in 1:length(famlist)){
  sites <- common_A1_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
}

# For loop to get dataframe for producing PCA
A1pca_common <- NULL
for(i in 1:length(famlist)){
  sites <- common_A1_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
  if(i<2){
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    A1pca_common <- pcasites
  }else{
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    A1pca_common <- cbind(A1pca_common,pcasites)
  }
  names(A1pca_common) <- gsub("_RF","", names(A1pca_common))
}

###########################   A5   ############################
## Refresh family groups list
famlist <- list(A5_fam=A5_fam,A5ssp1_fam=A5ssp1_fam, A5ssp2_fam=A5ssp2_fam,A5ssp3_fam=A5ssp3_fam,X21_fam=X21_fam,X27_fam=X27_fam,Z64_fam=Z64_fam,Z65_fam=Z65_fam, Z68_fam =Z68_fam)

## Filter all family groups for commonly occuring RAD SNP loci ####CHANGE sites
for(i in 1:length(famlist)){
  sites <- common_A5_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
}

## For loop to get dataframe for producing PCA
A5pca_common <- NULL
for(i in 1:length(famlist)){
  sites <- common_A5_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
  if(i<2){
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    A5pca_common <- pcasites
  }else{
  pcasites <- new[,seq(5,c(ncol(new)-11),4)]
  A5pca_common <- cbind(A5pca_common,pcasites)
  }
  names(A5pca_common) <- gsub("_RF","", names(A5pca_common))
}

###########################   B12   ############################
## Refresh family groups list
famlist <- list(B12_fam=B12_fam, B12ssp1_fam = B12ssp1_fam, B12ssp3_fam = B12ssp3_fam, B12ssp4_fam = B12ssp4_fam, B12ssp5_fam=B12ssp5_fam,Z31_fam =Z31_fam,Z34_fam =Z34_fam, Z34ssp2_fam =Z34ssp2_fam)

## Filter all family groups for commonly occuring RAD SNP loci ####CHANGE sites
for(i in 1:length(famlist)){
  sites <- common_B12_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
}

## For loop to get datafram for producing PCA
B12pca_common <- NULL
for(i in 1:length(famlist)){
  sites <- common_B12_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
  if(i<2){
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    B12pca_common <- pcasites
  }else{
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    B12pca_common <- cbind(B12pca_common,pcasites)
  }
  names(B12pca_common) <- gsub("_RF","", names(B12pca_common))
}

###########################   C2   ############################
## Refresh family groups list
#C5_fam =C5_fam,
famlist <- list( C2_fam =C2_fam,C2ssp1_fam =C2ssp1_fam, C2ssp2_fam =C2ssp2_fam, C2ssp3_fam=C2ssp3_fam, X35_fam =X35_fam, X36_fam =X36_fam, X38_fam=X38_fam, X39_fam= X39_fam)

## Filter all family groups for commonly occuring RAD SNP loci ####CHANGE sites
for(i in 1:length(famlist)){
  sites <- common_C2_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
}

## For loop to get datafram for producing PCA
C2pca_common <- NULL
for(i in 1:length(famlist)){
  sites <- common_C2_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
  if(i<2){
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    C2pca_common <- pcasites
  }else{
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    C2pca_common <- cbind(C2pca_common,pcasites)
  }
  names(C2pca_common) <- gsub("_RF","", names(C2pca_common))
}

###########################   C3   ############################
## Refresh family groups list
famlist <- list(C3_fam=C3_fam, C3ssp3_fam =C3ssp3_fam, Sc1_fam=Sc1_fam, Sc1a_fam =Sc1a_fam, Sc1b_fam =Sc1b_fam, Sc1c_fam =Sc1c_fam, Sc1d_fam=Sc1d_fam, Sc1e_fam =Sc1e_fam, Sc1f_fam= Sc1f_fam,Sc2_fam=Sc2_fam, Sc2a_fam =Sc2a_fam, Sc2b_fam =Sc2b_fam, Sc2c_fam=Sc2c_fam,Sc2d_fam =Sc2d_fam, Sc2e_fam =Sc2e_fam, X31_fam=X31_fam, Z12_fam =Z12_fam, Z15_fam =Z15_fam, Z16_fam =Z16_fam, Z17_fam =Z17_fam, Z18_fam =Z18_fam, Z19_fam =Z19_fam,Z41_fam =Z41_fam, Z43_fam =Z43_fam)
## 

## Filter all family groups for commonly occuring RAD SNP loci ####CHANGE sites
for(i in 1:length(famlist)){
  sites <- common_C3_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
}

## For loop to get datafram for producing PCA
C3pca_common <- NULL
for(i in 1:length(famlist)){
  sites <- common_C3_sites
  fam_name <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(famlist[i])), ' ')[[1]][1])
  new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(fam_name, "_common_sites"),new)
  if(i<2){
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    C3pca_common <- pcasites
  }else{
    pcasites <- new[,seq(5,c(ncol(new)-11),4)]
    C3pca_common <- cbind(C3pca_common,pcasites)
  }
  names(C3pca_common) <- gsub("_RF","", names(C3pca_common))
}


##################### Family group PCA ###################
library("tidyverse")
library("magrittr")
library("FactoMineR")
library("factoextra")
library("ggfortify")
A1pca_mat <- t(A1pca_common)
A5pca_mat <- t(A5pca_common)
B12pca_mat <- t(B12pca_common)
C2pca_mat <- t(C2pca_common)
C3pca_mat <- t(C3pca_common)
#A5pca_dist <- dist(A5pca_mat,method = "euclidian")
A1fam_pca <- PCA(A1pca_mat, scale.unit = TRUE, ncp = 20, graph = FALSE)
A5fam_pca <- PCA(A5pca_mat, scale.unit = TRUE, ncp = 20, graph = FALSE)
B12fam_pca <- PCA(B12pca_mat, scale.unit = TRUE, ncp = 20, graph = FALSE)
C2fam_pca <- PCA(C2pca_mat, scale.unit = TRUE, ncp = 20, graph = FALSE)
C3fam_pca <- PCA(C3pca_mat, scale.unit = TRUE, ncp = 20, graph = FALSE)

A1fam_pca <- prcomp(A1pca_mat, center=TRUE, scale. = TRUE)
A5fam_pca <- prcomp(A5pca_mat, center=TRUE, scale. = TRUE)
B12fam_pca <- prcomp(B12pca_mat, center=TRUE, scale. = TRUE)
C2fam_pca <- prcomp(C2pca_mat, center=TRUE, scale. = TRUE)
C3fam_pca <- prcomp(C3pca_mat, center=TRUE, scale. = TRUE)

ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE, 
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = muted('red'))
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  # Label the variable axes
  if(var.axes) {
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar, y = yvar, 
                    angle = angle, hjust = hjust), 
                color = 'darkred', size = varname.size)
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}


## Load sample data and generate figure
RAD_samples <- "../../Phd/RADseq/common_sites/new_added_lines/RAD_sample_list.txt"
RAD_sampleIDs <- read.table(file=RAD_samples, sep = "", header = FALSE)
colnames(RAD_sampleIDs) <- c("SampleID","Old_Name", "Parent_Group", "New_Name", "Isolate","Data_Source","Genetic_Group", "Nuclear_Genotype")
rownames(RAD_sampleIDs) <- RAD_sampleIDs$SampleID
RAD_sampleIDs$Genetic_Group <- as.factor(RAD_sampleIDs$Genetic_Group)
RAD_sampleIDs$Nuclear_Genotype <- as.factor(RAD_sampleIDs$Nuclear_Genotype)
A1RAD_names <- RAD_sampleIDs[which(rownames(RAD_sampleIDs) %in% colnames(A1pca_common)),]
A5RAD_names <- RAD_sampleIDs[which(rownames(RAD_sampleIDs) %in% colnames(A5pca_common)),]
B12RAD_names <- RAD_sampleIDs[which(rownames(RAD_sampleIDs) %in% colnames(B12pca_common)),]
C2RAD_names <- RAD_sampleIDs[which(rownames(RAD_sampleIDs) %in% colnames(C2pca_common)),]
C3RAD_names <- RAD_sampleIDs[which(rownames(RAD_sampleIDs) %in% colnames(C3pca_common)),]


ggbiplot(A1fam_pca,var.axes = FALSE, groups = A1RAD_names$Isolate,
         )+
  theme(axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color = "grey"),
        #panel.grid.minor=element_line(color = "black"),
        panel.background=element_blank())

ggbiplot(A5fam_pca,var.axes = FALSE, groups = A5RAD_names$Isolate,
)+
  theme(axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color = "grey"),
        #panel.grid.minor=element_line(color = "black"),
        panel.background=element_blank())

ggbiplot(C3fam_pca,var.axes = FALSE, groups = C3RAD_names$Isolate,
)+
  theme(axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color = "grey"),
        #panel.grid.minor=element_line(color = "black"),
        panel.background=element_blank())


#RAD_sites_ca <- CA(all_RAD_sites_mat)
A1fam_eigen <- get_eigenvalue(A1fam_pca)
A5fam_eigen <- get_eigenvalue(A5fam_pca)
B12fam_eigen <- get_eigenvalue(B12fam_pca)
C2fam_eigen <- get_eigenvalue(C2fam_pca)
C3fam_eigen <- get_eigenvalue(C3fam_pca)


fviz_pca_ind(A1fam_pca,
             geom.ind = "point",
             pointsize = 2,
             pointshape = 19,
             #col.var = "black",
             #repel = TRUE,
             alpha.ind = 0.5,
             col.ind = A1RAD_names$Isolate,
             palette = c("#CCCC00", "#33CC33", "#0066CC", "#FF0000"),
             #addEllipses = TRUE,
             legend.title = "Families") +
  theme(axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

fviz_pca_ind(A5fam_pca,
             geom.ind = "point",
             pointsize = 2,
             pointshape = 19,
             col.var = "black",
             #repel = TRUE,
             alpha.ind = 0.5,
             col.ind = A5RAD_names$Isolate,
             palette = c("#CCCC00", "#33CC33", "#0066CC", "#FF0000", "#FF9933", "#FF66FF", "#993333", "#00CCFF", "#6600CC"),
             addEllipses = TRUE,
             legend.title = "Families") +
  theme(axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

fviz_pca_ind(B12fam_pca,
             geom.ind = "point",
             pointsize = 2,
             pointshape = 19,
             col.var = "black",
             #repel = TRUE,
             alpha.ind = 0.5,
             col.ind = B12RAD_names$Isolate,
             palette = c("#CCCC00", "#33CC33", "#0066CC", "#FF0000", "#FF9933", "#FF66FF", "#993333", "#00CCFF", "#6600CC"),
             #addEllipses = TRUE,
             legend.title = "Families") +
  theme(axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

fviz_pca_ind(C2fam_pca,
             geom.ind = "point",
             pointsize = 2,
             pointshape = 19,
             #col.var = "black",
             #repel = TRUE,
             alpha.ind = 0.5,
             col.ind = C2RAD_names$Isolate,
             palette = c("#CCCC00", "#33CC33", "#0066CC", "#FF0000", "#FF9933", "#FF66FF", "#993333", "#00CCFF", "#6600CC"),
             #addEllipses = TRUE,
             legend.title = "Families") +
  theme(axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

fviz_pca_ind(C3fam_pca,
             geom.ind = "point",
             #pointsize = 2,
             pointshape = 19,
             #col.var = "black",
             repel = TRUE,
             alpha.ind = 0.5,
             col.ind = C3RAD_names$Isolate,
             #palette = c("#CCCC00", "#33CC33", "#0066CC", "#FF0000", "#FF9933", "#FF66FF", "#993333", "#00CCFF", "#6600CC"),
             #addEllipses = TRUE,
             legend.title = "Families") +
  theme(axis.line=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())


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

### Run tests on paired sets of family groups

### For loop to do all chi2 and ttest on sites between parent and progeny
### example:
### A5_A5p1_chi2 <- site_x2(A5_common_sites,A5ssp1_common_sites,alpha = 0.05)
### A5_A5p1_ttest <- site_test(A5_common_sites,A5ssp1_common_sites,alpha = 0.05)


###########################   A1   ############################
# A1_family -> A1ssp1, A1ssp2, Z25
A1fam <- list(A1ssp1_common_sites, A1ssp2_common_sites, Z25_common_sites)

for(i in 1:length(A1fam)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(A1_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(A1fam[[i]][3])), ' ')[[1]][1])
  sites <- A1_common_sites
  x2test <- site_x2(sites,A1fam[[i]],alpha = 0.05)
  tttest  <- site_test(sites,A1fam[[i]],alpha = 0.05)
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_chi2test"),x2test)
  assign(paste0(par,"_",prog,"_tttest"),tttest)
}

###########################   A5   ############################
# A5_family -> A5ssp1, A5ssp2, A5ssp3, X21, X27, Z64, Z65, Z68
A5fam <- list(A5ssp1_common_sites, A5ssp2_common_sites, A5ssp3_common_sites, X21_common_sites, X27_common_sites, Z64_common_sites, Z65_common_sites, Z68_common_sites)

for(i in 1:length(A5fam)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(A5_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(A5fam[[i]][3])), ' ')[[1]][1])
  sites <- A5_common_sites
  x2test <- site_x2(sites,A5fam[[i]],alpha = 0.05)
  tttest  <- site_test(sites,A5fam[[i]],alpha = 0.05)
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_chi2test"),x2test)
  assign(paste0(par,"_",prog,"_tttest"),tttest)
}

###########################   B12   ############################
# B12_family -> B12ssp1, B12ssp3, B12ssp4, B12ssp5, Z31, Z34, Z34ssp2

##using B12ssp1_common_sites for testing as parent
B12fam <- list(B12ssp1_common_sites, B12ssp3_common_sites, B12ssp4_common_sites, B12ssp5_common_sites, Z31_common_sites, Z34_common_sites, Z34ssp2_common_sites)

for(i in 1:length(B12fam)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(B12_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(B12fam[[i]][3])), ' ')[[1]][1])
  sites <- B12_common_sites
  x2test <- site_x2(sites,B12fam[[i]],alpha = 0.05)
  tttest  <- site_test(sites,B12fam[[i]],alpha = 0.05)
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_chi2test"),x2test)
  assign(paste0(par,"_",prog,"_tttest"),tttest)
}

###########################   C2   ############################
# C2_C5_family -> C2ssp1, C2ssp2, C2ssp3, X35, X36, X38, X39
# C5_common_sites,
C2fam <- list (C2ssp1_common_sites, C2ssp2_common_sites, C2ssp3_common_sites, X35_common_sites, X36_common_sites, X38_common_sites, X39_common_sites)

for(i in 1:length(C2fam)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(C2_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(C2fam[[i]][3])), ' ')[[1]][1])
  sites <- C2_common_sites
  x2test <- site_x2(sites,C2fam[[i]],alpha = 0.05)
  tttest  <- site_test(sites,C2fam[[i]],alpha = 0.05)
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_chi2test"),x2test)
  assign(paste0(par,"_",prog,"_tttest"),tttest)
}

###########################   C3   ############################
# C3_family -> C3ssp3, Sc1a, Sc1b, Sc1c, Sc1d, Sc1e, Sc1f, Sc2, Sc2a, Sc2b, Sc2d, Sc2e, X31, Z12, Z15, Z17, Z18, Z19, Z41, Z43
C3fam <- list(C3ssp3_common_sites, Sc1_common_sites, Sc1a_common_sites, Sc1b_common_sites, Sc1c_common_sites, Sc1d_common_sites, Sc1e_common_sites, Sc1f_common_sites, Sc2_common_sites, Sc2a_common_sites, Sc2b_common_sites, Sc2c_common_sites, Sc2d_common_sites, Sc2e_common_sites, Z12_common_sites, X31_common_sites, Z15_common_sites, Z16_common_sites, Z17_common_sites, Z18_common_sites, Z19_common_sites,
Z41_common_sites, Z43_common_sites)
##X31_common_sites, 
for(i in 1:length(C3fam)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(C3_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(C3fam[[i]][3])), ' ')[[1]][1])
  sites <- C3_common_sites
  x2test <- site_x2(sites,C3fam[[i]],alpha = 0.05)
  tttest  <- site_test(sites,C3fam[[i]],alpha = 0.05)
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_chi2test"),x2test)
  assign(paste0(par,"_",prog,"_tttest"),tttest)
}


## For checking COMMON SIGNIFICANT SITES (not really necessary for the moment)
#common_A5_chisq_sites <- Reduce(intersect, list(rownames(A5_A5ssp1_chi2test), rownames(A5_A5ssp2_chi2test), rownames(A5_A5ssp3_chi2test), rownames(A5_X21_chi2test), rownames(A5_X27_chi2test), rownames(A5_Z64_chi2test), rownames(A5_Z65_chi2test), rownames(A5_Z68_chi2test)))

#common_A5_ttest_sites <- Reduce(intersect, list(rownames(A5_A5p1_ttest), rownames(A5_A5p2_ttest), rownames(A5_A5p3_ttest), rownames(A5_X21_ttest), rownames(A5_X27_ttest), rownames(A5_Z64_ttest), rownames(A5_Z65_ttest), rownames(A5_Z68_ttest)))


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

###########################   A1   ############################
### Selecting significant sites in comparisons
# A1_family -> A1ssp1, A1ssp2, Z25
A1prog <- list(A1ssp1_common_sites, A1ssp2_common_sites, Z25_common_sites)
# chi2 test results
A1x2prog <- list(A1_A1ssp1_chi2test, A1_A1ssp2_chi2test, A1_Z25_chi2test)
# ttest results
A1ttprog <- list(A1_A1ssp1_tttest, A1_A1ssp2_tttest, A1_Z25_tttest)

### For loop for selecting significant sites for plotting and making heatmap matrix
### example:
### A1_A1ssp1_x2sites <- site_select(A1_common_sites,A1ssp1_common_sites,A1_A1ssp1_chi2test)
### A1_A1ssp1_ttsites <- site_select(A1_common_sites,A1ssp1_common_sites,A1_A1ssp1_tttest)
refx2_list <- NULL
reftt_list <- NULL
altx2_list <- NULL
alttt_list <- NULL
progeny <- NULL
for(i in 1:length(A1prog)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(A1_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(A1prog[[i]][3])), ' ')[[1]][1])
  allchi2 <- A1x2prog
  alltt <- A1ttprog
  sites <- A1_common_sites
  x2sites <- site_select(sites,A1prog[[i]],allchi2[[i]])
  ttsites <- site_select(sites,A1prog[[i]],alltt[[i]])
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_x2sites"),x2sites)
  assign(paste0(par,"_",prog,"_ttsites"),ttsites)
  if(i<2){
    parrefx2 <- median(x2sites[,1])-median(x2sites[,1])
    parreftt <- median(ttsites[,1])-median(ttsites[,1])
    paraltx2 <- median(x2sites[,2])-median(x2sites[,2])
    paralttt <- median(ttsites[,2])-median(ttsites[,2])
    
    progeny <- append(progeny,par)
    refx2_list <- append(refx2_list,parrefx2)
    reftt_list <- append(reftt_list,parreftt)
    altx2_list <- append(altx2_list,paraltx2)
    alttt_list <- append(alttt_list,paralttt)
    
    refx2 <- median(x2sites[,3])-median(x2sites[,1])
    reftt <- median(ttsites[,3])-median(ttsites[,1])
    altx2 <- median(x2sites[,4])-median(x2sites[,2])
    alttt <- median(ttsites[,4])-median(ttsites[,2])
     
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
  }else{
    refx2 <- median(x2sites[,3])-median(x2sites[,1])
    reftt <- median(ttsites[,3])-median(ttsites[,1])
    altx2 <- median(x2sites[,4])-median(x2sites[,2])
    alttt <- median(ttsites[,4])-median(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
    
    A1freqs <- data.frame(refx2_list,reftt_list,altx2_list,alttt_list)
    rownames(A1freqs) <- c(progeny)
    
    A1_freqmat <- t(A1freqs)
  }
}

## GENERATE HEATMAP of change in progeny frequencies compared to parent
library("pheatmap")
library("gtools")
my_paletteLength <- 101
my_palette <- colorRampPalette(c("red", "white", "green"))(n = my_paletteLength)
colnames(A1_freqmat) <- unique(RAD_sampleIDs[which(RAD_sampleIDs$Old_Name %in% colnames(A1_freqmat)),5])
A1_freqmat <- A1_freqmat[,mixedorder(colnames(A1_freqmat))]


pheatmap(A1_freqmat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         cellwidth = 70,
         cellheight = 70,
         gaps_col = 1,
         gaps_row = 2,
         color = my_palette)

## t.test whether allele frequency means from progeny are different from parent
t.test(A1_A1ssp1_x2sites[,3],mu=mean(A1_A1ssp1_x2sites[,1]))
t.test(A1_A1ssp2_x2sites[,3],mu=mean(A1_A1ssp2_x2sites[,1]))
t.test(A1_Z25_x2sites[,3],mu=mean(A1_Z25_x2sites[,1]))

t.test(A1_A1ssp1_ttsites[,3],mu=mean(A1_A1ssp1_ttsites[,1]))
t.test(A1_A1ssp2_ttsites[,3],mu=mean(A1_A1ssp2_ttsites[,1]))
t.test(A1_Z25_ttsites[,3],mu=mean(A1_Z25_ttsites[,1]))


###########################   A5   ############################
# A5_family -> A5ssp1, A5ssp2, A5ssp3, X21, X27, Z64, Z65, Z68
 A5prog <- list(A5ssp1_common_sites, A5ssp2_common_sites, A5ssp3_common_sites, X21_common_sites, X27_common_sites, Z64_common_sites, Z65_common_sites, Z68_common_sites)
# chi2 test results
 A5x2prog <- list(A5_A5ssp1_chi2test, A5_A5ssp2_chi2test, A5_A5ssp3_chi2test, A5_X21_chi2test, A5_X27_chi2test, A5_Z64_chi2test, A5_Z65_chi2test, A5_Z68_chi2test)
# ttest results
 A5ttprog <-  list(A5_A5ssp1_tttest, A5_A5ssp2_tttest, A5_A5ssp3_tttest, A5_X21_tttest, A5_X27_tttest, A5_Z64_tttest, A5_Z65_tttest,A5_Z68_tttest)

### For loop for selecting significant sites for plotting and making heatmap matrix
### example:
### A5_A5p1_x2sites <- site_select(A5_common_sites,A5ssp1_common_sites,A5_A5p1_chi2test)
### A5_A5p1_ttsites <- site_select(A5_common_sites,A5ssp1_common_sites,A5_A5p1_tttest)

 refx2_list <- NULL
 reftt_list <- NULL
 altx2_list <- NULL
 alttt_list <- NULL
 progeny <- NULL
for(i in 1:length(A5prog)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(A5_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(A5prog[[i]][3])), ' ')[[1]][1])
  allchi2 <- A5x2prog
  alltt <- A5ttprog
  sites <- A5_common_sites
  x2sites <- site_select(sites,A5prog[[i]],allchi2[[i]])
  ttsites <- site_select(sites,A5prog[[i]],alltt[[i]])
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_x2sites"),x2sites)
  assign(paste0(par,"_",prog,"_ttsites"),ttsites)
  if(i<2){
    parrefx2 <- mean(x2sites[,1])-mean(x2sites[,1])
    parreftt <- mean(ttsites[,1])-mean(ttsites[,1])
    paraltx2 <- mean(x2sites[,2])-mean(x2sites[,2])
    paralttt <- mean(ttsites[,2])-mean(ttsites[,2])
    
    progeny <- append(progeny,par)
    refx2_list <- append(refx2_list,parrefx2)
    reftt_list <- append(reftt_list,parreftt)
    altx2_list <- append(altx2_list,paraltx2)
    alttt_list <- append(alttt_list,paralttt)
    
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
  }else{
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
    
    A5freqs <- data.frame(refx2_list,reftt_list,altx2_list,alttt_list)
    rownames(A5freqs) <- c(progeny)
    
    A5_freqmat <- t(A5freqs)
  }
}

 ## GENERATE HEATMAP of change in progeny frequencies compared to parent
 library("pheatmap")
 my_paletteLength <- 101
 my_palette <- colorRampPalette(c("red", "white", "green"))(n = my_paletteLength)
 colnames(A5_freqmat) <- unique(RAD_sampleIDs[which(RAD_sampleIDs$Old_Name %in% colnames(A5_freqmat)),5])
 A5_freqmat <- A5_freqmat[,mixedorder(colnames(A5_freqmat))]
 
 
 pheatmap(A5_freqmat, 
          cluster_rows = FALSE, 
          cluster_cols = FALSE,
          cellwidth = 70,
          cellheight = 70,
          gaps_col = 1,
          gaps_row = 2,
          color = my_palette)
 
 ## Simple density plot for allele frequencies
 plot(density(A5_A5ssp1_x2sites$A5_AF_Median,bw=0.1), col="orange3",main="Significantly different alleles between A5 and A5ssp1", xlab="Frequency", ylab="Density", xlim=c(0, 1), ylim=c(0, 4))
 lines(density(A5_A5ssp1_x2sites$A5ssp1_AF_Median, bw=0.1), col="orange")
 abline(v=0.5, col="red",lty=1)
 abline(v=mean(A5_A5ssp1_x2sites$A5_AF_Median), col="orange3",lty=3)
 abline(v=mean(A5_A5ssp1_x2sites$A5ssp1_AF_Median), col="orange",lty=3)

## t.test whether allele frequency medians from progeny are different from parent
t.test(Z65_common_sites[,ncol(Z65_common_sites)-7],mu=mean(A5_common_sites[,ncol(A5_common_sites)-7]))

t.test(A5_A5ssp1_x2sites[,3],mu=mean(A5_A5ssp1_x2sites[,1]))
t.test(A5_A5ssp2_x2sites[,3],mu=mean(A5_A5ssp2_x2sites[,1]))
t.test(A5_A5ssp3_x2sites[,3],mu=mean(A5_A5ssp3_x2sites[,1]))
t.test(A5_X21_x2sites[,3],mu=mean(A5_X21_x2sites[,1]))
t.test(A5_X27_x2sites[,3],mu=mean(A5_X27_x2sites[,1]))
t.test(A5_Z64_x2sites[,3],mu=mean(A5_Z64_x2sites[,1]))
t.test(A5_Z65_x2sites[,3],mu=mean(A5_Z65_x2sites[,1]))
t.test(A5_Z68_x2sites[,3],mu=mean(A5_Z68_x2sites[,1]))

t.test(A5_A5ssp1_ttsites[,3],mu=mean(A5_A5ssp1_ttsites[,1]))
t.test(A5_A5ssp2_ttsites[,3],mu=mean(A5_A5ssp2_ttsites[,1]))
t.test(A5_A5ssp3_ttsites[,3],mu=mean(A5_A5ssp3_ttsites[,1]))
t.test(A5_X21_ttsites[,3],mu=mean(A5_X21_ttsites[,1]))
t.test(A5_X27_ttsites[,3],mu=mean(A5_X27_ttsites[,1]))
t.test(A5_Z64_ttsites[,3],mu=mean(A5_Z64_ttsites[,1]))
t.test(A5_Z65_ttsites[,3],mu=mean(A5_Z65_ttsites[,1]))
t.test(A5_Z68_ttsites[,3],mu=mean(A5_Z68_ttsites[,1]))


###########################   B12   ############################
# B12_family -> B12ssp1, B12ssp3, B12ssp4, B12ssp5, Z31, Z34, Z34ssp2
###using B12ssp1_common_sites for comparisons
B12prog <- list( B12ssp1_common_sites,B12ssp3_common_sites, B12ssp4_common_sites, B12ssp5_common_sites, Z31_common_sites, Z34_common_sites, Z34ssp2_common_sites)
# chi2 test results
##B12_B12ssp1_chi2test,
B12x2prog <- list(B12_B12ssp1_chi2test, B12_B12ssp3_chi2test, B12_B12ssp4_chi2test, B12_B12ssp5_chi2test, B12_Z31_chi2test, B12_Z34_chi2test, B12_Z34ssp2_chi2test)
# ttest results
##B12_B12ssp1_tttest, 
B12ttprog <-  list(B12_B12ssp1_tttest, B12_B12ssp3_tttest, B12_B12ssp4_tttest, B12_B12ssp5_tttest, B12_Z31_tttest, B12_Z34_tttest, B12_Z34ssp2_tttest)

### For loop for selecting significant sites for plotting and making heatmap matrix
### example:
### B12_B12ssp1_x2sites <- site_select(B12_common_sites,B12ssp1_common_sites,B12_B12ssp1_chi2test)
### B12_B12ssp1_ttsites <- site_select(B12_common_sites,B12ssp1_common_sites,B12_B12ssp1_tttest)

refx2_list <- NULL
reftt_list <- NULL
altx2_list <- NULL
alttt_list <- NULL
progeny <- NULL
###using B12ssp1 as comparison
for(i in 1:length(B12prog)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(B12_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(B12prog[[i]][3])), ' ')[[1]][1])
  allchi2 <- B12x2prog
  alltt <- B12ttprog
  sites <- B12_common_sites
  x2sites <- site_select(sites,B12prog[[i]],allchi2[[i]])
  ttsites <- site_select(sites,B12prog[[i]],alltt[[i]])
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_x2sites"),x2sites)
  assign(paste0(par,"_",prog,"_ttsites"),ttsites)
  if(i<2){
    parrefx2 <- mean(x2sites[,1])-mean(x2sites[,1])
    parreftt <- mean(ttsites[,1])-mean(ttsites[,1])
    paraltx2 <- mean(x2sites[,2])-mean(x2sites[,2])
    paralttt <- mean(ttsites[,2])-mean(ttsites[,2])
    
    progeny <- append(progeny,par)
    refx2_list <- append(refx2_list,parrefx2)
    reftt_list <- append(reftt_list,parreftt)
    altx2_list <- append(altx2_list,paraltx2)
    alttt_list <- append(alttt_list,paralttt)
    
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
  }else{
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
    
    B12freqs <- data.frame(refx2_list,reftt_list,altx2_list,alttt_list)
    rownames(B12freqs) <- c(progeny)
    
    B12_freqmat <- t(B12freqs)
  }
}

## GENERATE HEATMAP of change in progeny frequencies compared to parent
library("pheatmap")
my_paletteLength <- 101
my_palette <- colorRampPalette(c("red", "white", "green"))(n = my_paletteLength)
colnames(B12_freqmat) <- unique(RAD_sampleIDs[which(RAD_sampleIDs$Old_Name %in% colnames(B12_freqmat)),5])
B12_freqmat <- B12_freqmat[,mixedorder(colnames(B12_freqmat))]

pheatmap(B12_freqmat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         cellwidth = 70,
         cellheight = 70,
         gaps_col = 1,
         gaps_row = 2,
         color = my_palette)

## t.test whether allele frequency medians from progeny are different from parent
t.test(B12_B12ssp1_x2sites[,3],mu=mean(B12_B12ssp1_x2sites[,1]))
t.test(B12_B12ssp3_x2sites[,3],mu=mean(B12_B12ssp3_x2sites[,1]))
t.test(B12_B12ssp4_x2sites[,3],mu=mean(B12_B12ssp4_x2sites[,1]))
t.test(B12_B12ssp5_x2sites[,3],mu=mean(B12_B12ssp5_x2sites[,1]))
t.test(B12_Z31_x2sites[,3],mu=mean(B12_Z31_x2sites[,1]))
t.test(B12_Z34_x2sites[,3],mu=mean(B12_Z34_x2sites[,1]))
t.test(B12_Z34ssp2_x2sites[,3],mu=mean(B12_Z34ssp2_x2sites[,1]))

t.test(B12_B12ssp1_ttsites[,3],mu=mean(B12_B12ssp1_ttsites[,1]))
t.test(B12_B12ssp3_ttsites[,3],mu=mean(B12_B12ssp3_ttsites[,1]))
t.test(B12_B12ssp4_ttsites[,3],mu=mean(B12_B12ssp4_ttsites[,1]))
t.test(B12_B12ssp5_ttsites[,3],mu=mean(B12_B12ssp5_ttsites[,1]))
t.test(B12_Z31_ttsites[,3],mu=mean(B12_Z31_ttsites[,1]))
t.test(B12_Z34_ttsites[,3],mu=mean(B12_Z34_ttsites[,1]))
t.test(B12_Z34ssp2_ttsites[,3],mu=mean(B12_Z34ssp2_ttsites[,1]))

###########################   C2   ############################
# C2_C5_family -> C2ssp1, C2ssp2, C2ssp3, X35, X36, X38, X39
#C5_common_sites, 
C2prog <- list(C2ssp1_common_sites, C2ssp2_common_sites, C2ssp3_common_sites, X35_common_sites, X36_common_sites, X38_common_sites, X39_common_sites)
# chi2 test results
#C2_C5_chi2test,
C2x2prog <- list( C2_C2ssp1_chi2test, C2_C2ssp2_chi2test, C2_C2ssp3_chi2test, C2_X35_chi2test, C2_X36_chi2test, C2_X38_chi2test, C2_X39_chi2test)
# ttest results
#C2_C5_tttest,
C2ttprog <-  list( C2_C2ssp1_tttest, C2_C2ssp2_tttest, C2_C2ssp3_tttest, C2_X35_tttest, C2_X36_tttest, C2_X38_tttest, C2_X39_tttest)

### For loop for selecting significant sites for plotting and making heatmap matrix
### example:
### C2_C5_x2sites <- site_select(C2_common_sites,C5_common_sites,C2_C5_chi2test)
### C2_C5_ttsites <- site_select(C2_common_sites,C5_common_sites,C2_C5_tttest)

refx2_list <- NULL
reftt_list <- NULL
altx2_list <- NULL
alttt_list <- NULL
progeny <- NULL
for(i in 1:length(C2prog)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(C2_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(C2prog[[i]][3])), ' ')[[1]][1])
  allchi2 <- C2x2prog
  alltt <- C2ttprog
  sites <- C2_common_sites
  x2sites <- site_select(sites,C2prog[[i]],allchi2[[i]])
  ttsites <- site_select(sites,C2prog[[i]],alltt[[i]])
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_x2sites"),x2sites)
  assign(paste0(par,"_",prog,"_ttsites"),ttsites)
  if(i<2){
    parrefx2 <- mean(x2sites[,1])-mean(x2sites[,1])
    parreftt <- mean(ttsites[,1])-mean(ttsites[,1])
    paraltx2 <- mean(x2sites[,2])-mean(x2sites[,2])
    paralttt <- mean(ttsites[,2])-mean(ttsites[,2])
    
    progeny <- append(progeny,par)
    refx2_list <- append(refx2_list,parrefx2)
    reftt_list <- append(reftt_list,parreftt)
    altx2_list <- append(altx2_list,paraltx2)
    alttt_list <- append(alttt_list,paralttt)
    
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
  }else{
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
    
    C2freqs <- data.frame(refx2_list,reftt_list,altx2_list,alttt_list)
    rownames(C2freqs) <- c(progeny)
    
    C2_freqmat <- t(C2freqs)
  }
}

## GENERATE HEATMAP of change in progeny frequencies compared to parent
library("pheatmap")
my_paletteLength <- 101
my_palette <- colorRampPalette(c("red", "white", "green"))(n = my_paletteLength)
colnames(C2_freqmat) <- unique(RAD_sampleIDs[which(RAD_sampleIDs$Old_Name %in% colnames(C2_freqmat)),5])
C2_freqmat <- C2_freqmat[,mixedorder(colnames(C2_freqmat))]

pheatmap(C2_freqmat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         cellwidth = 70,
         cellheight = 70,
         gaps_col = 1,
         gaps_row = 2,
         color = my_palette)

## t.test whether allele frequency medians from progeny are different from parent
# C2_C5_family -> C2ssp1, C2ssp2, C2ssp3, X35, X36, X38, X39

t.test(X36_common_sites[,ncol(X36_common_sites)-7],C2_common_sites[,ncol(C2_common_sites)-7])

#t.test(C2_C5_x2sites[,3],mu=mean(C2_C5_x2sites[,1]))
t.test(C2_C2ssp1_x2sites[,3],mu=mean(C2_C2ssp1_x2sites[,1]))
t.test(C2_C2ssp2_x2sites[,3],mu=mean(C2_C2ssp2_x2sites[,1]))
t.test(C2_C2ssp3_x2sites[,3],mu=mean(C2_C2ssp3_x2sites[,1]))
t.test(C2_X35_x2sites[,3],mu=mean(C2_X35_x2sites[,1]))
t.test(C2_X36_x2sites[,3],mu=mean(C2_X36_x2sites[,1]))
t.test(C2_X38_x2sites[,3],mu=mean(C2_X38_x2sites[,1]))
t.test(C2_X39_x2sites[,3],mu=mean(C2_X39_x2sites[,1]))

t.test(C2_C5_ttsites[,3],mu=mean(C2_C5_ttsites[,1]))
t.test(C2_C2ssp1_ttsites[,3],mu=mean(C2_C2ssp1_ttsites[,1]))
t.test(C2_C2ssp2_ttsites[,3],mu=mean(C2_C2ssp2_ttsites[,1]))
t.test(C2_C2ssp3_ttsites[,3],mu=mean(C2_C2ssp3_ttsites[,1]))
t.test(C2_X35_ttsites[,3],mu=mean(C2_X35_ttsites[,1]))
t.test(C2_X36_ttsites[,3],mu=mean(C2_X36_ttsites[,1]))
t.test(C2_X38_ttsites[,3],mu=mean(C2_X38_ttsites[,1]))
t.test(C2_X39_ttsites[,3],mu=mean(C2_X39_ttsites[,1]))


###########################   C3   ############################
# C3_family -> C3ssp3, Sc1a, Sc1b, Sc1c, Sc1d, Sc1e, Sc1f, Sc2, Sc2a, Sc2b, Sc2d, Sc2e, X31, Z12, Z15, Z17, Z18, Z19, Z41, Z43
C3prog <- list(C3ssp3_common_sites, Sc1_common_sites, Sc1a_common_sites, Sc1b_common_sites, Sc1c_common_sites, Sc1d_common_sites, Sc1e_common_sites, Sc1f_common_sites, Sc2_common_sites, Sc2a_common_sites, Sc2b_common_sites, Sc2c_common_sites, Sc2d_common_sites, Sc2e_common_sites, X31_common_sites, Z12_common_sites, Z15_common_sites, Z16_common_sites, Z17_common_sites, Z18_common_sites, Z19_common_sites, Z41_common_sites, Z43_common_sites)
##X31_common_sites, 
# chi2 test results
C3x2prog <- list(C3_C3ssp3_chi2test, C3_Sc1_chi2test, C3_Sc1a_chi2test, C3_Sc1b_chi2test, C3_Sc1c_chi2test, C3_Sc1d_chi2test, C3_Sc1e_chi2test, C3_Sc1f_chi2test, C3_Sc2_chi2test, C3_Sc2a_chi2test, C3_Sc2b_chi2test, C3_Sc2c_chi2test, C3_Sc2d_chi2test, C3_Sc2e_chi2test, C3_X31_chi2test, C3_Z12_chi2test, C3_Z15_chi2test, C3_Z16_chi2test, C3_Z17_chi2test, C3_Z18_chi2test, C3_Z19_chi2test, C3_Z41_chi2test, C3_Z43_chi2test)
## C3_X31_chi2test, 
# ttest results
C3ttprog <-  list(C3_C3ssp3_tttest, C3_Sc1_tttest, C3_Sc1a_tttest, C3_Sc1b_tttest, C3_Sc1c_tttest, C3_Sc1d_tttest, C3_Sc1e_tttest, C3_Sc1f_tttest, C3_Sc2_tttest, C3_Sc2a_tttest, C3_Sc2b_tttest, C3_Sc2c_tttest, C3_Sc2d_tttest, C3_Sc2e_tttest, C3_X31_tttest, C3_Z12_tttest, C3_Z15_tttest,  C3_Z16_tttest, C3_Z17_tttest, C3_Z18_tttest, C3_Z19_tttest,C3_Z41_tttest, C3_Z43_tttest)
##  C3_X31_tttest, 

### For loop for selecting significant sites for plotting and making heatmap matrix
### example:
### C3_C3ssp3_x2sites <- site_select(C3_common_sites,C3ssp3_common_sites,C3_C3ssp3_chi2test)
### C3_C3ssp3_ttsites <- site_select(C3_common_sites,C3ssp3_common_sites,C3_C3ssp3_tttest)

refx2_list <- NULL
reftt_list <- NULL
altx2_list <- NULL
alttt_list <- NULL
progeny <- NULL
for(i in 1:length(C3prog)){
  par <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', deparse(substitute(C3_common_sites))), ' ')[[1]][1])
  prog <- noquote(strsplit(sub('(^[^_]+)_(.*)$', '\\1 \\2', names(C3prog[[i]][3])), ' ')[[1]][1])
  allchi2 <- C3x2prog
  alltt <- C3ttprog
  sites <- C3_common_sites
  x2sites <- site_select(sites,C3prog[[i]],allchi2[[i]])
  ttsites <- site_select(sites,C3prog[[i]],alltt[[i]])
  #new <- famlist[[i]][which(rownames(famlist[[i]]) %in% sites),]
  assign(paste0(par,"_",prog,"_x2sites"),x2sites)
  assign(paste0(par,"_",prog,"_ttsites"),ttsites)
  if(i<2){
    parrefx2 <- mean(x2sites[,1])-mean(x2sites[,1])
    parreftt <- mean(ttsites[,1])-mean(ttsites[,1])
    paraltx2 <- mean(x2sites[,2])-mean(x2sites[,2])
    paralttt <- mean(ttsites[,2])-mean(ttsites[,2])
    
    progeny <- append(progeny,par)
    refx2_list <- append(refx2_list,parrefx2)
    reftt_list <- append(reftt_list,parreftt)
    altx2_list <- append(altx2_list,paraltx2)
    alttt_list <- append(alttt_list,paralttt)
    
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
  }else{
    refx2 <- mean(x2sites[,3])-mean(x2sites[,1])
    reftt <- mean(ttsites[,3])-mean(ttsites[,1])
    altx2 <- mean(x2sites[,4])-mean(x2sites[,2])
    alttt <- mean(ttsites[,4])-mean(ttsites[,2])
    
    progeny <- append(progeny,prog)
    refx2_list <- append(refx2_list,refx2)
    reftt_list <- append(reftt_list,reftt)
    altx2_list <- append(altx2_list,altx2)
    alttt_list <- append(alttt_list,alttt)
    
    C3freqs <- data.frame(refx2_list,reftt_list,altx2_list,alttt_list)
    rownames(C3freqs) <- c(progeny)
    
    C3_freqmat <- t(C3freqs)
  }
}

## GENERATE HEATMAP of change in progeny frequencies compared to parent
library("pheatmap")
library("gtools")
my_paletteLength <- 101
my_palette <- colorRampPalette(c("red", "white", "green"))(n = my_paletteLength)
colnames(C3_freqmat) <- unique(RAD_sampleIDs[which(RAD_sampleIDs$Old_Name %in% colnames(C3_freqmat)),5])
C3_freqmat <- C3_freqmat[,mixedorder(colnames(C3_freqmat))]
c3_freqmat <- mixedorder(colnames(C3_freqmat))
c3_freqnew <- C3_freqmat[,c3_freqmat]

pheatmap(C3_freqmat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         cellwidth = 70,
         cellheight = 70,
         gaps_col = 1,
         gaps_row = 2,
         color = my_palette)

##save as 

## t.test whether allele frequency medians from progeny are different from parent
# C3_family -> C3ssp3, Sc1a, Sc1b, Sc1c, Sc1d, Sc1e, Sc1f, Sc2, Sc2a, Sc2b, Sc2d, Sc2e, X31, Z12, Z15, Z17, Z18, Z19, Z41, Z43
t.test(Sc2d_common_sites[,ncol(Sc2d_common_sites)-7], mu= mean(C3_common_sites[,ncol(C3_common_sites)-7]))

t.test(C3_C3ssp3_x2sites[,3],mu=mean(C3_C3ssp3_x2sites[,1]))
t.test(C3_Sc1_x2sites[,3], mu=mean(C3_Sc1_x2sites[,1]))
t.test(C3_Sc1a_x2sites[,3],mu=mean(C3_Sc1a_x2sites[,1]))
t.test(C3_Sc1b_x2sites[,3],mu=mean(C3_Sc1b_x2sites[,1]))
t.test(C3_Sc1c_x2sites[,3],mu=mean(C3_Sc1c_x2sites[,1]))
t.test(C3_Sc1d_x2sites[,3],mu=mean(C3_Sc1d_x2sites[,1]))
t.test(C3_Sc1e_x2sites[,3],mu=mean(C3_Sc1e_x2sites[,1]))
t.test(C3_Sc1f_x2sites[,3],mu=mean(C3_Sc1f_x2sites[,1]))
t.test(C3_Sc2_x2sites[,3],mu=mean(C3_Sc2_x2sites[,1]))
t.test(C3_Sc2a_x2sites[,3],mu=mean(C3_Sc2a_x2sites[,1]))
t.test(C3_Sc2b_x2sites[,3],mu=mean(C3_Sc2b_x2sites[,1]))
t.test(C3_Sc2c_x2sites[,3],mu=mean(C3_Sc2c_x2sites[,1]))
t.test(C3_Sc2d_x2sites[,3],mu=mean(C3_Sc2d_x2sites[,1]))
t.test(C3_Sc2e_x2sites[,3],mu=mean(C3_Sc2e_x2sites[,1]))
t.test(C3_X31_x2sites[,3],mu=mean(C3_X31_x2sites[,1]))
t.test(C3_Z12_x2sites[,3],mu=mean(C3_Z12_x2sites[,1]))
t.test(C3_Z15_x2sites[,3],mu=mean(C3_Z15_x2sites[,1]))
t.test(C3_Z16_x2sites[,3],mu=mean(C3_Z16_x2sites[,1]))
t.test(C3_Z17_x2sites[,3],mu=mean(C3_Z17_x2sites[,1]))
t.test(C3_Z18_x2sites[,3],mu=mean(C3_Z18_x2sites[,1]))
t.test(C3_Z19_x2sites[,3],mu=mean(C3_Z19_x2sites[,1]))
t.test(C3_Z41_x2sites[,3],mu=mean(C3_Z41_x2sites[,1]))
t.test(C3_Z43_x2sites[,3],mu=mean(C3_Z43_x2sites[,1]))

t.test(C3_C3ssp3_ttsites[,3],mu=mean(C3_C3ssp3_ttsites[,1]))
t.test(C3_Sc1_ttsites[,3],mu=mean(C3_Sc1_ttsites[,1]))
t.test(C3_Sc1a_ttsites[,3],mu=mean(C3_Sc1a_ttsites[,1]))
t.test(C3_Sc1b_ttsites[,3],mu=mean(C3_Sc1b_ttsites[,1]))
t.test(C3_Sc1c_ttsites[,3],mu=mean(C3_Sc1c_ttsites[,1]))
t.test(C3_Sc1d_ttsites[,3],mu=mean(C3_Sc1d_ttsites[,1]))
t.test(C3_Sc1e_ttsites[,3],mu=mean(C3_Sc1e_ttsites[,1]))
t.test(C3_Sc1f_ttsites[,3],mu=mean(C3_Sc1f_ttsites[,1]))
t.test(C3_Sc2_ttsites[,3],mu=mean(C3_Sc2_ttsites[,1]))
t.test(C3_Sc2a_ttsites[,3],mu=mean(C3_Sc2a_ttsites[,1]))
t.test(C3_Sc2b_ttsites[,3],mu=mean(C3_Sc2b_ttsites[,1]))
t.test(C3_Sc2c_ttsites[,3],mu=mean(C3_Sc2c_ttsites[,1]))
t.test(C3_Sc2d_ttsites[,3],mu=mean(C3_Sc2d_ttsites[,1]))
t.test(C3_Sc2e_ttsites[,3],mu=mean(C3_Sc2e_ttsites[,1]))
t.test(C3_X31_ttsites[,3],mu=mean(C3_X31_ttsites[,1]))
t.test(C3_Z12_ttsites[,3],mu=mean(C3_Z12_ttsites[,1]))
t.test(C3_Z15_ttsites[,3],mu=mean(C3_Z15_ttsites[,1]))
t.test(C3_Z16_ttsites[,3],mu=mean(C3_Z16_ttsites[,1]))
t.test(C3_Z17_ttsites[,3],mu=mean(C3_Z17_ttsites[,1]))
t.test(C3_Z18_ttsites[,3],mu=mean(C3_Z18_ttsites[,1]))
t.test(C3_Z19_ttsites[,3],mu=mean(C3_Z19_ttsites[,1]))
t.test(C3_Z41_ttsites[,3],mu=mean(C3_Z41_ttsites[,1]))
t.test(C3_Z43_ttsites[,3],mu=mean(C3_Z43_ttsites[,1]))


##########################Box-whisker plot for site count##############################
setwd("/Users/wrobbins/Desktop/Phd/RADseq/common_sites/")
dir <- "/Users/wrobbins/Desktop/Phd/RADseq/common_sites/"

RAD_sig_loci <- "significant_sites_counts.txt"
Significant_loci <- read.table(file=RAD_sig_loci, sep = "", header = FALSE)
colnames(Significant_loci) <- c("Isolate", "Parent", "Nuclear_Genotype", "Significant_loci","Common_loci")
Significant_loci$Nuclear_Genotype <- factor(Significant_loci$Nuclear_Genotype, levels = c("mono", "di"))

p <-ggplot(Significant_loci, aes(x=Nuclear_Genotype, y=log(Significant_loci))) +
  geom_boxplot( aes(Nuclear_Genotype)) +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2), alpha=0.5)+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


t.test(log(Significant_loci[which(Significant_loci$Nuclear_Genotype == c("mono")),4]),mu=mean(log(Significant_loci[which(Significant_loci$Nuclear_Genotype == c("di")),4])))


############Create randomized sample of common sites for regression#################
#### A1 family
names(sample(famlist, 4))
common2famsites_A1 <- Reduce(intersect, list(rownames(A1_fam), rownames(A1ssp1_fam)))
common4famsites_A1 <- Reduce(intersect, list( rownames(Z25_fam),rownames(A1ssp1_fam) ,rownames(A1ssp1_fam),  rownames(A1_fam)))

A1_2famsites <- Reduce(intersect, list(rownames(A1_fam), rownames(A1ssp1_fam), rownames(A1ssp2_fam), rownames(Z25_fam)))

#### A5 family
names(sample(famlist, 8))
common8famsites_A5 <- Reduce(intersect,  list(rownames(A5_fam), rownames(X21_fam), rownames(A5ssp3_fam), rownames(A5ssp2_fam), rownames(X27_fam), rownames(Z64_fam), rownames(A5ssp1_fam), rownames(Z68_fam)))

common_A5_sites <- Reduce(intersect, list(rownames(A5_fam), rownames(A5ssp1_fam), rownames(A5ssp3_fam), rownames(A5ssp2_fam), rownames(X21_fam), rownames(X27_fam), rownames(Z64_fam), rownames(Z65_fam), rownames(Z68_fam)))

#### B12 family
names(sample(famlist, 6))
common6famsites_B12 <- Reduce(intersect, list(rownames(B12ssp1_fam), rownames(B12ssp3_fam),rownames(B12ssp5_fam), rownames(B12ssp4_fam), rownames(Z34_fam), rownames(Z31_fam)))


common_B12_sites <- Reduce(intersect, list(rownames(B12ssp1_fam), rownames(B12ssp3_fam), rownames(B12ssp4_fam), rownames(B12ssp5_fam), rownames(Z31_fam), rownames(Z34_fam), rownames(Z34ssp2_fam)))

#### C2 family
names(sample(famlist, 6))
common6famsites_C2 <- Reduce(intersect, list(  rownames(C2ssp1_fam), rownames(C2ssp3_fam),rownames(C2ssp2_fam), rownames(X36_fam), rownames(X38_fam), rownames(X39_fam)))

common_C2_sites <- Reduce(intersect, list( rownames(C2_fam), rownames(C2ssp1_fam), rownames(C2ssp2_fam), rownames(C2ssp3_fam), rownames(X35_fam), rownames(X36_fam), rownames(X38_fam), rownames(X39_fam)))


#### C3 family
names(sample(famlist, 20))
common8famsites_C3 <- Reduce(intersect, list(rownames(C3_fam), rownames(C3ssp3_fam),  rownames(Sc1a_fam), rownames(Sc1b_fam), rownames(Z18_fam), rownames(Sc1d_fam), rownames(Sc1e_fam), rownames(Sc1f_fam),rownames(Sc2_fam), rownames(Sc2a_fam), rownames(Sc2b_fam), rownames(Sc2d_fam),  rownames(Z12_fam),  rownames(Z15_fam), rownames(Z16_fam), rownames(Z17_fam), rownames(Z18_fam), rownames(Sc2e_fam),rownames(Z41_fam), rownames(Z43_fam)))

common_C3_sites <- Reduce(intersect, list(rownames(C3_fam), rownames(C3ssp3_fam),  rownames(Sc1a_fam), rownames(Sc1b_fam), rownames(Sc1c_fam), rownames(Sc1d_fam), rownames(Sc1e_fam), rownames(Sc1f_fam),rownames(Sc2_fam), rownames(Sc2a_fam), rownames(Sc2b_fam), rownames(Sc2d_fam), rownames(Sc2e_fam), rownames(Z12_fam),  rownames(Z15_fam), rownames(Z16_fam), rownames(Z17_fam), rownames(Z18_fam), rownames(Z19_fam),rownames(Z41_fam), rownames(Z43_fam)))

sites_per_RAD_sample <- "line_number_effect_common_sites.txt"
sites_per_RAD_sample <- read.table(file=sites_per_RAD_sample, sep = "", header = FALSE)
colnames(sites_per_RAD_sample) <- c("Parent", "Genetic_group", "Isolate_number","Nuclear_genotype", "Detected_sites")

monocommonisolates <- sites_per_RAD_sample[which(sites_per_RAD_sample$Nuclear_genotype == c("mono")),]
fit_loci_monoisolates <- lm(Detected_sites ~Isolate_number,data=monocommonisolates)
summary(fit_loci_monoisolates)

dicommonisolates <- sites_per_RAD_sample[which(sites_per_RAD_sample$Nuclear_genotype == c("di")),]
fit_loci_diisolates <- lm(Detected_sites ~Isolate_number,data=dicommonisolates)
summary(fit_loci_diisolates)


p <-ggplot(sites_per_RAD_sample, aes(x=Isolate_number, y=Detected_sites)) +
  geom_point( aes(colour=Parent)) +
  geom_smooth(aes(fill=Nuclear_genotype, color=Nuclear_genotype), method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


######################Loci locations################################
###A1
A1_gene_coordinates <- "A1_gene_coordinates.txt"
A1_gene_coords <- read.table(file=A1_gene_coordinates, sep = "", fill = TRUE,header = FALSE)
A1_gene_coords <- A1_gene_coords[,1:6]
colnames(A1_gene_coords) <- c("Sca", "Start", "End","Body", "Strand", "GeneID","Gene_length")
A1_gene_coords[,7] <- A1_gene_coords$End + 1 - A1_gene_coords$Start

A1_sites <- A1_common_sites[,1:2]
A1_site_list <- NULL
for(i in 1:nrow(A1_sites)){
  site <- A1_sites[i,]
  new_site <- A1_gene_coords[which(A1_gene_coords$Sca == site$Sca),]
  if(new_site$Strand=="+"){
    new_site[,8] <- (site[,2] + 1 - new_site$Start)/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }else{
    new_site[,8] <- (new_site$End + 1 - site[,2])/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }
  
  A1_site_list <- rbind(A1_site_list,new_site)
}

A1_gene_count <- A1_site_list[which(A1_site_list$Body==c("transcript")),]
A1_intron_count <- A1_site_list[which(A1_site_list$Body==c("intron")),]
A1_initial_count <- A1_site_list[which(A1_site_list$Body==c("initial")),]
A1_internal_count <- A1_site_list[which(A1_site_list$Body==c("internal")),]
A1_terminal_count <- A1_site_list[which(A1_site_list$Body==c("terminal")),]
A1_single_count <- A1_site_list[which(A1_site_list$Body==c("single")),]


###A5
A5_gene_coordinates <- "A5_gene_coordinates.txt"
A5_gene_coords <- read.table(file=A5_gene_coordinates, sep = "", fill = TRUE,header = FALSE)
A5_gene_coords <- A5_gene_coords[,1:6]
colnames(A5_gene_coords) <- c("Sca", "Start", "End","Body", "Strand", "GeneID","Gene_length")
A5_gene_coords[,7] <- A5_gene_coords$End + 1 - A5_gene_coords$Start

A5_sites <- A5_common_sites[,1:2]
A5_site_list <- NULL
for(i in 1:nrow(A5_sites)){
  site <- A5_sites[i,]
  new_site <- A5_gene_coords[which(A5_gene_coords$Sca == site$Sca),]
  if(new_site$Strand=="+"){
    new_site[,8] <- (site[,2] + 1 - new_site$Start)/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }else{
    new_site[,8] <- (new_site$End + 1 - site[,2])/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }
  
  A5_site_list <- rbind(A5_site_list,new_site)
}

A5_gene_count <- A5_site_list[which(A5_site_list$Body==c("transcript")),]
A5_intron_count <- A5_site_list[which(A5_site_list$Body==c("intron")),]
A5_initial_count <- A5_site_list[which(A5_site_list$Body==c("initial")),]
A5_internal_count <- A5_site_list[which(A5_site_list$Body==c("internal")),]
A5_terminal_count <- A5_site_list[which(A5_site_list$Body==c("terminal")),]
A5_single_count <- A5_site_list[which(A5_site_list$Body==c("single")),]


###B12
B12_gene_coordinates <- "B12_gene_coordinates.txt"
B12_gene_coords <- read.table(file=B12_gene_coordinates, sep = "", fill = TRUE,header = FALSE)
B12_gene_coords <- B12_gene_coords[,1:6]
colnames(B12_gene_coords) <- c("Sca", "Start", "End","Body", "Strand", "GeneID","Gene_length")
B12_gene_coords[,7] <- B12_gene_coords$End + 1 - B12_gene_coords$Start

B12_sites <- B12ssp1_common_sites[,1:2]
B12_site_list <- NULL
for(i in 1:nrow(B12_sites)){
  site <- B12_sites[i,]
  new_site <- B12_gene_coords[which(B12_gene_coords$Sca == site$Sca),]
  if(new_site$Strand=="+"){
    new_site[,8] <- (site[,2] + 1 - new_site$Start)/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }else{
    new_site[,8] <- (new_site$End + 1 - site[,2])/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }
  
  B12_site_list <- rbind(B12_site_list,new_site)
}

B12_gene_count <- B12_site_list[which(B12_site_list$Body==c("transcript")),]
B12_intron_count <- B12_site_list[which(B12_site_list$Body==c("intron")),]
B12_initial_count <- B12_site_list[which(B12_site_list$Body==c("initial")),]
B12_internal_count <- B12_site_list[which(B12_site_list$Body==c("internal")),]
B12_terminal_count <- B12_site_list[which(B12_site_list$Body==c("terminal")),]
B12_single_count <- B12_site_list[which(B12_site_list$Body==c("single")),]

###C2
C2_gene_coordinates <- "C2_gene_coordinates.txt"
C2_gene_coords <- read.table(file=C2_gene_coordinates, sep = "", fill = TRUE,header = FALSE)
C2_gene_coords <- C2_gene_coords[,1:6]
C2_gene_coords[,7] <- C2_gene_coords$End + 1 - C2_gene_coords$Start
colnames(C2_gene_coords) <- c("Sca", "Start", "End","Body", "Strand", "GeneID","Gene_length")

C2_sites <- C2_common_sites[,1:2]
C2_site_list <- NULL
for(i in 1:nrow(C2_sites)){
  site <- C2_sites[i,]
  new_site <- C2_gene_coords[which(C2_gene_coords$Sca == site$Sca),]
  if(new_site$Strand=="+"){
    new_site[,8] <- (site[,2] + 1 - new_site$Start)/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
      
  }else{
    new_site[,8] <- (new_site$End + 1 - site[,2])/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }
  
  C2_site_list <- rbind(C2_site_list,new_site)
}

C2_gene_countn <- C2_site_list[which(C2_site_list$Body==c("mRNA")),]
C2_intron_countn <- C2_site_list[which(C2_site_list$Body==c("intron")),]
C2_initial_countn <- C2_site_list[which(C2_site_list$Body==c("initial")),]
C2_internal_countn <- C2_site_list[which(C2_site_list$Body==c("internal")),]
C2_terminal_countn <- C2_site_list[which(C2_site_list$Body==c("terminal")),]
C2_single_countn <- C2_site_list[which(C2_site_list$Body==c("single")),]

###C3
C3_gene_coordinates <- "C3_gene_coordinates.txt"
C3_gene_coords <- read.table(file=C3_gene_coordinates, sep = "", fill = TRUE,header = FALSE)
C3_gene_coords <- C3_gene_coords[,1:6]
C3_gene_coords[,7] <- C3_gene_coords$End + 1 - C3_gene_coords$Start
colnames(C3_gene_coords) <- c("Sca", "Start", "End","Body", "Strand", "GeneID","Gene_length")

C3_sites <- C3_common_sites[,1:2]
C3_site_list <- NULL
for(i in 1:nrow(C3_sites)){
  site <- C3_sites[i,]
  new_site <- C3_gene_coords[which(C3_gene_coords$Sca == site$Sca),]
  if(new_site$Strand=="+"){
    new_site[,8] <- (site[,2] + 1 - new_site$Start)/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }else{
    new_site[,8] <- (new_site$End + 1 - site[,2])/new_site[,7]
    new_site[,9] <- site[,2]
    new_site <- new_site[which(new_site$Start <= site$Pos),]
    new_site <- new_site[which(new_site$End >= site$Pos),]
  }

  
  C3_site_list <- rbind(C3_site_list,new_site)
}

C3_gene_count <- C3_site_list[which(C3_site_list$Body==c("transcript")),]
C3_intron_count <- C3_site_list[which(C3_site_list$Body==c("intron")),]
C3_initial_count <- C3_site_list[which(C3_site_list$Body==c("initial")),]
C3_internal_count <- C3_site_list[which(C3_site_list$Body==c("internal")),]
C3_terminal_count <- C3_site_list[which(C3_site_list$Body==c("terminal")),]
C3_single_count <- C3_site_list[which(C3_site_list$Body==c("single")),]


####### plots for gene positions or shared sites for isolates 
plot(density(C3_gene_count$V8,bw=0.01), col="darkblue", xlab="Gene_proportion", xlim=c(0, 1))
lines(density(A5_gene_count$V8, bw=0.01), col="blue")

plot(density(C2_gene_countn$V8,bw=0.01), col="darkgreen", xlab="Gene_proportion", xlim=c(0, 1))
lines(density(A1_gene_count$V8, bw=0.01), col="green")
lines(density(B12_gene_count$V8, bw=0.01), col="lightgreen")


###C3
RAD_position_genes <- "RAD_loci_position_identifiers.txt"
RAD_pos_genes <- read.table(file=RAD_position_genes, sep = "",header = FALSE)
colnames(RAD_pos_genes) <- c("Isolate", "Nuclear_genotype", "Region","Body", "Feature", "Added_feature","Loci", "Percent")
RAD_pos_genes$Nuclear_genotype <- factor(RAD_pos_genes$Nuclear_genotype, levels = c("mono", "di"))
RAD_pos_genes$Region <- factor(RAD_pos_genes$Region, levels = c("noncoding", "coding"))
RAD_pos_genes$Isolate <- factor(RAD_pos_genes$Isolate, levels = c("A1", "B12","C2","A5","C3"))
RAD_pos_genes$Feature <- factor(RAD_pos_genes$Feature, levels = c("intron", "initial","internal","terminal","noncoding"))

RAD_pos <- RAD_pos_genes[which(RAD_pos_genes$Nuclear_genotype=="mono" & RAD_pos_genes$Region=="coding"),]

p <-ggplot(RAD_pos, aes(x=Isolate, y=Loci, fill = Nuclear_genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Region)+
  scale_y_continuous(limits = c(0,1750))+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


p <- ggplot(RAD_pos, aes(x = Isolate, y = Percent, fill = Feature)) +
  geom_col() +
  scale_x_discrete(limits = c(" "," ", "", "A1","B12","C2")) +
  scale_fill_viridis_d() +
  coord_polar("y")+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
       # legend.position='none',
        axis.line= element_line(colour="black"))
p

##### monokaryote and dikaryote box plot
###A1
#A1prog <- list(A1ssp1_common_sites, A1ssp2_common_sites, Z25_common_sites)
A1ssp1_sig_sites <- A1ssp1_common_sites[which(rownames(A1ssp1_common_sites) %in% rownames(A1_A1ssp1_chi2test)),]
A1ssp2_sig_sites <- A1ssp2_common_sites[which(rownames(A1ssp2_common_sites) %in% rownames(A1_A1ssp2_chi2test)),]
Z25_sig_sites <- Z25_common_sites[which(rownames(Z25_common_sites) %in% rownames(A1_Z25_chi2test)),]

A1_gene_coordinates <- "A1_gene_coordinates.txt"
A1_gene_coords <- read.table(file=A1_gene_coordinates, sep = "", fill = TRUE,header = FALSE)
A1_gene_coords <- A1_gene_coords[,1:6]
colnames(A1_gene_coords) <- c("Sca", "Start", "End","Body", "Strand", "GeneID")

Z25_sites <- Z25_sig_sites[,1:2]
Z25_site_list <- NULL
for(i in 1:nrow(Z25_sites)){
  site <- Z25_sites[i,]
  new_site <- A1_gene_coords[which(A1_gene_coords$Sca == site$Sca),]
  new_site <- new_site[which(new_site$Start <= site$Pos),]
  new_site <- new_site[which(new_site$End >= site$Pos),]
  Z25_site_list <- rbind(Z25_site_list,new_site)
}

A1ssp1_gene_count <- A1ssp1_site_list[which(A1ssp1_site_list$Body==c("transcript")),]
A1ssp1_intron_count <- A1ssp1_site_list[which(A1ssp1_site_list$Body==c("intron")),]
A1ssp1_initial_count <- A1ssp1_site_list[which(A1ssp1_site_list$Body==c("initial")),]
A1ssp1_internal_count <- A1ssp1_site_list[which(A1ssp1_site_list$Body==c("internal")),]
A1ssp1_terminal_count <- A1ssp1_site_list[which(A1ssp1_site_list$Body==c("terminal")),]
A1ssp1_single_count <- A1ssp1_site_list[which(A1ssp1_site_list$Body==c("single")),]

A1ssp2_gene_count <- A1ssp2_site_list[which(A1ssp2_site_list$Body==c("transcript")),]
A1ssp2_intron_count <- A1ssp2_site_list[which(A1ssp2_site_list$Body==c("intron")),]
A1ssp2_initial_count <- A1ssp2_site_list[which(A1ssp2_site_list$Body==c("initial")),]
A1ssp2_internal_count <- A1ssp2_site_list[which(A1ssp2_site_list$Body==c("internal")),]
A1ssp2_terminal_count <- A1ssp2_site_list[which(A1ssp2_site_list$Body==c("terminal")),]
A1ssp2_single_count <- A1ssp2_site_list[which(A1ssp2_site_list$Body==c("single")),]

Z25_gene_count <- Z25_site_list[which(Z25_site_list$Body==c("transcript")),]
Z25_intron_count <- Z25_site_list[which(Z25_site_list$Body==c("intron")),]
Z25_initial_count <- Z25_site_list[which(Z25_site_list$Body==c("initial")),]
Z25_internal_count <- Z25_site_list[which(Z25_site_list$Body==c("internal")),]
Z25_terminal_count <- Z25_site_list[which(Z25_site_list$Body==c("terminal")),]
Z25_single_count <- Z25_site_list[which(Z25_site_list$Body==c("single")),]

########## B12
#B12prog <- list( B12ssp3_common_sites, B12ssp4_common_sites, B12ssp5_common_sites, Z31_common_sites, Z34_common_sites, Z34ssp2_common_sites)

B12ssp3_sig_sites <- B12ssp3_common_sites[which(rownames(B12ssp3_common_sites) %in% rownames(B12ssp1_B12ssp3_chi2test)),]
B12ssp4_sig_sites <- B12ssp4_common_sites[which(rownames(B12ssp4_common_sites) %in% rownames(B12ssp1_B12ssp4_chi2test)),]
B12ssp5_sig_sites <- B12ssp5_common_sites[which(rownames(B12ssp5_common_sites) %in% rownames(B12ssp1_B12ssp5_chi2test)),]
Z31_sig_sites <- Z31_common_sites[which(rownames(Z31_common_sites) %in% rownames(B12ssp1_Z31_chi2test)),]
Z34_sig_sites <- Z34_common_sites[which(rownames(Z34_common_sites) %in% rownames(B12ssp1_Z34_chi2test)),]
Z34ssp2_sig_sites <- Z34ssp2_common_sites[which(rownames(Z34ssp2_common_sites) %in% rownames(B12ssp1_Z34ssp2_chi2test)),]

B12_gene_coordinates <- "B12_gene_coordinates.txt"
B12_gene_coords <- read.table(file=B12_gene_coordinates, sep = "", fill = TRUE,header = FALSE)
B12_gene_coords <- B12_gene_coords[,1:6]
colnames(B12_gene_coords) <- c("Sca", "Start", "End","Body", "Strand", "GeneID")

Z34ssp2_sites <- Z34ssp2_sig_sites[,1:2]
Z34ssp2_site_list <- NULL
for(i in 1:nrow(Z34ssp2_sites)){
  site <- Z34ssp2_sites[i,]
  new_site <- B12_gene_coords[which(B12_gene_coords$Sca == site$Sca),]
  new_site <- new_site[which(new_site$Start <= site$Pos),]
  new_site <- new_site[which(new_site$End >= site$Pos),]
  Z34ssp2_site_list <- rbind(Z34ssp2_site_list,new_site)
}

B12ssp3_gene_count <- B12ssp3_site_list[which(B12ssp3_site_list$Body==c("transcript")),]
B12ssp3_intron_count <- B12ssp3_site_list[which(B12ssp3_site_list$Body==c("intron")),]
B12ssp3_initial_count <- B12ssp3_site_list[which(B12ssp3_site_list$Body==c("initial")),]
B12ssp3_internal_count <- B12ssp3_site_list[which(B12ssp3_site_list$Body==c("internal")),]
B12ssp3_terminal_count <- B12ssp3_site_list[which(B12ssp3_site_list$Body==c("terminal")),]
B12ssp3_single_count <- B12ssp3_site_list[which(B12ssp3_site_list$Body==c("single")),]

B12ssp4_gene_count <- B12ssp4_site_list[which(B12ssp4_site_list$Body==c("transcript")),]
B12ssp4_intron_count <- B12ssp4_site_list[which(B12ssp4_site_list$Body==c("intron")),]
B12ssp4_initial_count <- B12ssp4_site_list[which(B12ssp4_site_list$Body==c("initial")),]
B12ssp4_internal_count <- B12ssp4_site_list[which(B12ssp4_site_list$Body==c("internal")),]
B12ssp4_terminal_count <- B12ssp4_site_list[which(B12ssp4_site_list$Body==c("terminal")),]
B12ssp4_single_count <- B12ssp4_site_list[which(B12ssp4_site_list$Body==c("single")),]

B12ssp5_gene_count <- B12ssp5_site_list[which(B12ssp5_site_list$Body==c("transcript")),]
B12ssp5_intron_count <- B12ssp5_site_list[which(B12ssp5_site_list$Body==c("intron")),]
B12ssp5_initial_count <- B12ssp5_site_list[which(B12ssp5_site_list$Body==c("initial")),]
B12ssp5_internal_count <- B12ssp5_site_list[which(B12ssp5_site_list$Body==c("internal")),]
B12ssp5_terminal_count <- B12ssp5_site_list[which(B12ssp5_site_list$Body==c("terminal")),]
B12ssp5_single_count <- B12ssp5_site_list[which(B12ssp5_site_list$Body==c("single")),]

Z31_gene_count <- Z31_site_list[which(Z31_site_list$Body==c("transcript")),]
Z31_intron_count <- Z31_site_list[which(Z31_site_list$Body==c("intron")),]
Z31_initial_count <- Z31_site_list[which(Z31_site_list$Body==c("initial")),]
Z31_internal_count <- Z31_site_list[which(Z31_site_list$Body==c("internal")),]
Z31_terminal_count <- Z31_site_list[which(Z31_site_list$Body==c("terminal")),]
Z31_single_count <- Z31_site_list[which(Z31_site_list$Body==c("single")),]

Z34_gene_count <- Z34_site_list[which(Z34_site_list$Body==c("transcript")),]
Z34_intron_count <- Z34_site_list[which(Z34_site_list$Body==c("intron")),]
Z34_initial_count <- Z34_site_list[which(Z34_site_list$Body==c("initial")),]
Z34_internal_count <- Z34_site_list[which(Z34_site_list$Body==c("internal")),]
Z34_terminal_count <- Z34_site_list[which(Z34_site_list$Body==c("terminal")),]
Z34_single_count <- Z34_site_list[which(Z34_site_list$Body==c("single")),]

Z34ssp2_gene_count <- Z34ssp2_site_list[which(Z34ssp2_site_list$Body==c("transcript")),]
Z34ssp2_intron_count <- Z34ssp2_site_list[which(Z34ssp2_site_list$Body==c("intron")),]
Z34ssp2_initial_count <- Z34ssp2_site_list[which(Z34ssp2_site_list$Body==c("initial")),]
Z34ssp2_internal_count <- Z34ssp2_site_list[which(Z34ssp2_site_list$Body==c("internal")),]
Z34ssp2_terminal_count <- Z34ssp2_site_list[which(Z34ssp2_site_list$Body==c("terminal")),]
Z34ssp2_single_count <- Z34ssp2_site_list[which(Z34ssp2_site_list$Body==c("single")),]

###########################   C2   #############################
C2ssp1_sig_sites <- C2ssp1_common_sites[which(rownames(C2ssp1_common_sites) %in% rownames(C2_C2ssp1_chi2test)),]
C2ssp2_sig_sites <- C2ssp2_common_sites[which(rownames(C2ssp2_common_sites) %in% rownames(C2_C2ssp2_chi2test)),] 
C2ssp3_sig_sites <- C2ssp3_common_sites[which(rownames(C2ssp3_common_sites) %in% rownames(C2_C2ssp3_chi2test)),] 
X35_sig_sites <- X35_common_sites[which(rownames(X35_common_sites) %in% rownames(C2_X35_chi2test)),]  
X36_sig_sites <- X36_common_sites[which(rownames(X36_common_sites) %in% rownames(C2_X36_chi2test)),]  
X38_sig_sites <- X38_common_sites[which(rownames(X38_common_sites) %in% rownames(C2_X38_chi2test)),]  
X39_sig_sites <- X39_common_sites[which(rownames(X39_common_sites) %in% rownames(C2_X39_chi2test)),] 

X39_sites <- X39_sig_sites[,1:2]
X39_site_list <- NULL
for(i in 1:nrow(X39_sites)){
  site <- X39_sites[i,]
  new_site <- C2_gene_coords[which(C2_gene_coords$Sca == site$Sca),]
  new_site <- new_site[which(new_site$Start <= site$Pos),]
  new_site <- new_site[which(new_site$End >= site$Pos),]
  X39_site_list <- rbind(X39_site_list,new_site)
}

C2ssp1_gene_count <- C2ssp1_site_list[which(C2ssp1_site_list$Body==c("mRNA")),]
C2ssp1_intron_count <- C2ssp1_site_list[which(C2ssp1_site_list$Body==c("intron")),]
C2ssp1_initial_count <- C2ssp1_site_list[which(C2ssp1_site_list$Body==c("initial")),]
C2ssp1_internal_count <- C2ssp1_site_list[which(C2ssp1_site_list$Body==c("internal")),]
C2ssp1_terminal_count <- C2ssp1_site_list[which(C2ssp1_site_list$Body==c("terminal")),]
C2ssp1_single_count <- C2ssp1_site_list[which(C2ssp1_site_list$Body==c("single")),]

C2ssp2_gene_count <- C2ssp2_site_list[which(C2ssp2_site_list$Body==c("mRNA")),]
C2ssp2_intron_count <- C2ssp2_site_list[which(C2ssp2_site_list$Body==c("intron")),]
C2ssp2_initial_count <- C2ssp2_site_list[which(C2ssp2_site_list$Body==c("initial")),]
C2ssp2_internal_count <- C2ssp2_site_list[which(C2ssp2_site_list$Body==c("internal")),]
C2ssp2_terminal_count <- C2ssp2_site_list[which(C2ssp2_site_list$Body==c("terminal")),]
C2ssp2_single_count <- C2ssp2_site_list[which(C2ssp2_site_list$Body==c("single")),]

C2ssp3_gene_count <- C2ssp3_site_list[which(C2ssp3_site_list$Body==c("mRNA")),]
C2ssp3_intron_count <- C2ssp3_site_list[which(C2ssp3_site_list$Body==c("intron")),]
C2ssp3_initial_count <- C2ssp3_site_list[which(C2ssp3_site_list$Body==c("initial")),]
C2ssp3_internal_count <- C2ssp3_site_list[which(C2ssp3_site_list$Body==c("internal")),]
C2ssp3_terminal_count <- C2ssp3_site_list[which(C2ssp3_site_list$Body==c("terminal")),]
C2ssp3_single_count <- C2ssp3_site_list[which(C2ssp3_site_list$Body==c("single")),]

X35_gene_count <- X35_site_list[which(X35_site_list$Body==c("mRNA")),]
X35_intron_count <- X35_site_list[which(X35_site_list$Body==c("intron")),]
X35_initial_count <- X35_site_list[which(X35_site_list$Body==c("initial")),]
X35_internal_count <- X35_site_list[which(X35_site_list$Body==c("internal")),]
X35_terminal_count <- X35_site_list[which(X35_site_list$Body==c("terminal")),]
X35_single_count <- X35_site_list[which(X35_site_list$Body==c("single")),]

X36_gene_count <- X36_site_list[which(X36_site_list$Body==c("mRNA")),]
X36_intron_count <- X36_site_list[which(X36_site_list$Body==c("intron")),]
X36_initial_count <- X36_site_list[which(X36_site_list$Body==c("initial")),]
X36_internal_count <- X36_site_list[which(X36_site_list$Body==c("internal")),]
X36_terminal_count <- X36_site_list[which(X36_site_list$Body==c("terminal")),]
X36_single_count <- X36_site_list[which(X36_site_list$Body==c("single")),]

X38_gene_count <- X38_site_list[which(X38_site_list$Body==c("mRNA")),]
X38_intron_count <- X38_site_list[which(X38_site_list$Body==c("intron")),]
X38_initial_count <- X38_site_list[which(X38_site_list$Body==c("initial")),]
X38_internal_count <- X38_site_list[which(X38_site_list$Body==c("internal")),]
X38_terminal_count <- X38_site_list[which(X38_site_list$Body==c("terminal")),]
X38_single_count <- X38_site_list[which(X38_site_list$Body==c("single")),]

X39_gene_count <- X39_site_list[which(X39_site_list$Body==c("mRNA")),]
X39_intron_count <- X39_site_list[which(X39_site_list$Body==c("intron")),]
X39_initial_count <- X39_site_list[which(X39_site_list$Body==c("initial")),]
X39_internal_count <- X39_site_list[which(X39_site_list$Body==c("internal")),]
X39_terminal_count <- X39_site_list[which(X39_site_list$Body==c("terminal")),]
X39_single_count <- X39_site_list[which(X39_site_list$Body==c("single")),]

###########################   A5   ############################
# A5_family -> A5ssp1, A5ssp2, A5ssp3, X21, X27, Z64, Z65, Z68
A5ssp1_sig_sites <- A5ssp1_common_sites[which(rownames(A5ssp1_common_sites) %in% rownames(A5_A5ssp1_chi2test)),]
A5ssp2_sig_sites <- A5ssp2_common_sites[which(rownames(A5ssp2_common_sites) %in% rownames(A5_A5ssp2_chi2test)),]
A5ssp3_sig_sites <- A5ssp3_common_sites[which(rownames(A5ssp3_common_sites) %in% rownames(A5_A5ssp3_chi2test)),]
X21_sig_sites <- X21_common_sites[which(rownames(X21_common_sites) %in% rownames(A5_X21_chi2test)),]
X27_sig_sites <- X27_common_sites[which(rownames(X27_common_sites) %in% rownames(A5_X27_chi2test)),]
Z64_sig_sites <- Z64_common_sites[which(rownames(Z64_common_sites) %in% rownames(A5_Z64_chi2test)),]
Z65_sig_sites <- Z65_common_sites[which(rownames(Z65_common_sites) %in% rownames(A5_Z65_chi2test)),]
Z68_sig_sites <- Z68_common_sites[which(rownames(Z68_common_sites) %in% rownames(A5_Z68_chi2test)),]

Z68_sites <- Z68_sig_sites[,1:2]
Z68_site_list <- NULL
for(i in 1:nrow(Z68_sites)){
  site <- Z68_sites[i,]
  new_site <- A5_gene_coords[which(A5_gene_coords$Sca == site$Sca),]
  new_site <- new_site[which(new_site$Start <= site$Pos),]
  new_site <- new_site[which(new_site$End >= site$Pos),]
  Z68_site_list <- rbind(Z68_site_list,new_site)
}

A5ssp1_gene_count <- A5ssp1_site_list[which(A5ssp1_site_list$Body==c("transcript")),]
A5ssp1_intron_count <- A5ssp1_site_list[which(A5ssp1_site_list$Body==c("intron")),]
A5ssp1_initial_count <- A5ssp1_site_list[which(A5ssp1_site_list$Body==c("initial")),]
A5ssp1_internal_count <- A5ssp1_site_list[which(A5ssp1_site_list$Body==c("internal")),]
A5ssp1_terminal_count <- A5ssp1_site_list[which(A5ssp1_site_list$Body==c("terminal")),]
A5ssp1_single_count <- A5ssp1_site_list[which(A5ssp1_site_list$Body==c("single")),]

A5ssp2_gene_count <- A5ssp2_site_list[which(A5ssp2_site_list$Body==c("transcript")),]
A5ssp2_intron_count <- A5ssp2_site_list[which(A5ssp2_site_list$Body==c("intron")),]
A5ssp2_initial_count <- A5ssp2_site_list[which(A5ssp2_site_list$Body==c("initial")),]
A5ssp2_internal_count <- A5ssp2_site_list[which(A5ssp2_site_list$Body==c("internal")),]
A5ssp2_terminal_count <- A5ssp2_site_list[which(A5ssp2_site_list$Body==c("terminal")),]
A5ssp2_single_count <- A5ssp2_site_list[which(A5ssp2_site_list$Body==c("single")),]

A5ssp3_gene_count <- A5ssp3_site_list[which(A5ssp3_site_list$Body==c("transcript")),]
A5ssp3_intron_count <- A5ssp3_site_list[which(A5ssp3_site_list$Body==c("intron")),]
A5ssp3_initial_count <- A5ssp3_site_list[which(A5ssp3_site_list$Body==c("initial")),]
A5ssp3_internal_count <- A5ssp3_site_list[which(A5ssp3_site_list$Body==c("internal")),]
A5ssp3_terminal_count <- A5ssp3_site_list[which(A5ssp3_site_list$Body==c("terminal")),]
A5ssp3_single_count <- A5ssp3_site_list[which(A5ssp3_site_list$Body==c("single")),]

X21_gene_count <- X21_site_list[which(X21_site_list$Body==c("transcript")),]
X21_intron_count <- X21_site_list[which(X21_site_list$Body==c("intron")),]
X21_initial_count <- X21_site_list[which(X21_site_list$Body==c("initial")),]
X21_internal_count <- X21_site_list[which(X21_site_list$Body==c("internal")),]
X21_terminal_count <- X21_site_list[which(X21_site_list$Body==c("terminal")),]
X21_single_count <- X21_site_list[which(X21_site_list$Body==c("single")),]

X27_gene_count <- X27_site_list[which(X27_site_list$Body==c("transcript")),]
X27_intron_count <- X27_site_list[which(X27_site_list$Body==c("intron")),]
X27_initial_count <- X27_site_list[which(X27_site_list$Body==c("initial")),]
X27_internal_count <- X27_site_list[which(X27_site_list$Body==c("internal")),]
X27_terminal_count <- X27_site_list[which(X27_site_list$Body==c("terminal")),]
X27_single_count <- X27_site_list[which(X27_site_list$Body==c("single")),]

Z64_gene_count <- Z64_site_list[which(Z64_site_list$Body==c("transcript")),]
Z64_intron_count <- Z64_site_list[which(Z64_site_list$Body==c("intron")),]
Z64_initial_count <- Z64_site_list[which(Z64_site_list$Body==c("initial")),]
Z64_internal_count <- Z64_site_list[which(Z64_site_list$Body==c("internal")),]
Z64_terminal_count <- Z64_site_list[which(Z64_site_list$Body==c("terminal")),]
Z64_single_count <- Z64_site_list[which(Z64_site_list$Body==c("single")),]

Z65_gene_count <- Z65_site_list[which(Z65_site_list$Body==c("transcript")),]
Z65_intron_count <- Z65_site_list[which(Z65_site_list$Body==c("intron")),]
Z65_initial_count <- Z65_site_list[which(Z65_site_list$Body==c("initial")),]
Z65_internal_count <- Z65_site_list[which(Z65_site_list$Body==c("internal")),]
Z65_terminal_count <- Z65_site_list[which(Z65_site_list$Body==c("terminal")),]
Z65_single_count <- Z65_site_list[which(Z65_site_list$Body==c("single")),]

Z68_gene_count <- Z68_site_list[which(Z68_site_list$Body==c("transcript")),]
Z68_intron_count <- Z68_site_list[which(Z68_site_list$Body==c("intron")),]
Z68_initial_count <- Z68_site_list[which(Z68_site_list$Body==c("initial")),]
Z68_internal_count <- Z68_site_list[which(Z68_site_list$Body==c("internal")),]
Z68_terminal_count <- Z68_site_list[which(Z68_site_list$Body==c("terminal")),]
Z68_single_count <- Z68_site_list[which(Z68_site_list$Body==c("single")),]

###########################   C3   ############################
# C3_family -> C3ssp3, Sc1a, Sc1b, Sc1c, Sc1d, Sc1e, Sc1f, Sc2, Sc2a, Sc2b, Sc2d, Sc2e, Z12, Z15, Z16, Z17, Z18, Z19, Z41, Z43 (X31)
C3ssp3_sig_sites <- C3ssp3_common_sites[which(rownames(C3ssp3_common_sites) %in% rownames(C3_C3ssp3_chi2test)),]
Sc1a_sig_sites <- Sc1a_common_sites[which(rownames(Sc1a_common_sites) %in% rownames(C3_Sc1a_chi2test)),]
Sc1b_sig_sites <- Sc1b_common_sites[which(rownames(Sc1b_common_sites) %in% rownames(C3_Sc1b_chi2test)),]
Sc1c_sig_sites <- Sc1c_common_sites[which(rownames(Sc1c_common_sites) %in% rownames(C3_Sc1c_chi2test)),]
Sc1d_sig_sites <- Sc1d_common_sites[which(rownames(Sc1d_common_sites) %in% rownames(C3_Sc1d_chi2test)),]
Sc1e_sig_sites <- Sc1e_common_sites[which(rownames(Sc1e_common_sites) %in% rownames(C3_Sc1e_chi2test)),]
Sc1f_sig_sites <- Sc1f_common_sites[which(rownames(Sc1f_common_sites) %in% rownames(C3_Sc1f_chi2test)),]
Sc2_sig_sites <- Sc2_common_sites[which(rownames(Sc2_common_sites) %in% rownames(C3_Sc2_chi2test)),]
Sc2a_sig_sites <- Sc2a_common_sites[which(rownames(Sc2a_common_sites) %in% rownames(C3_Sc2a_chi2test)),]
Sc2b_sig_sites <- Sc2b_common_sites[which(rownames(Sc2b_common_sites) %in% rownames(C3_Sc2b_chi2test)),]
Sc2d_sig_sites <- Sc2d_common_sites[which(rownames(Sc2d_common_sites) %in% rownames(C3_Sc2d_chi2test)),]
Sc2e_sig_sites <- Sc2e_common_sites[which(rownames(Sc2e_common_sites) %in% rownames(C3_Sc2e_chi2test)),]
Z12_sig_sites <- Z12_common_sites[which(rownames(Z12_common_sites) %in% rownames(C3_Z12_chi2test)),]
Z15_sig_sites <- Z15_common_sites[which(rownames(Z15_common_sites) %in% rownames(C3_Z15_chi2test)),]
Z16_sig_sites <- Z16_common_sites[which(rownames(Z16_common_sites) %in% rownames(C3_Z16_chi2test)),]
Z17_sig_sites <- Z17_common_sites[which(rownames(Z17_common_sites) %in% rownames(C3_Z17_chi2test)),]
Z18_sig_sites <- Z18_common_sites[which(rownames(Z18_common_sites) %in% rownames(C3_Z18_chi2test)),]
Z19_sig_sites <- Z19_common_sites[which(rownames(Z19_common_sites) %in% rownames(C3_Z19_chi2test)),]
Z41_sig_sites <- Z41_common_sites[which(rownames(Z41_common_sites) %in% rownames(C3_Z41_chi2test)),]
Z43_sig_sites <- Z43_common_sites[which(rownames(Z43_common_sites) %in% rownames(C3_Z43_chi2test)),]

Z43_sites <- Z43_sig_sites[,1:2]
Z43_site_list <- NULL
for(i in 1:nrow(Z43_sites)){
  site <- Z43_sites[i,]
  new_site <- C3_gene_coords[which(C3_gene_coords$Sca == site$Sca),]
  new_site <- new_site[which(new_site$Start <= site$Pos),]
  new_site <- new_site[which(new_site$End >= site$Pos),]
  Z43_site_list <- rbind(Z43_site_list,new_site)
}

C3ssp3_gene_count <- C3ssp3_site_list[which(C3ssp3_site_list$Body==c("transcript")),]
C3ssp3_intron_count <- C3ssp3_site_list[which(C3ssp3_site_list$Body==c("intron")),]
C3ssp3_initial_count <- C3ssp3_site_list[which(C3ssp3_site_list$Body==c("initial")),]
C3ssp3_internal_count <- C3ssp3_site_list[which(C3ssp3_site_list$Body==c("internal")),]
C3ssp3_terminal_count <- C3ssp3_site_list[which(C3ssp3_site_list$Body==c("terminal")),]
C3ssp3_single_count <- C3ssp3_site_list[which(C3ssp3_site_list$Body==c("single")),]

Sc1a_gene_count <- Sc1a_site_list[which(Sc1a_site_list$Body==c("transcript")),]
Sc1a_intron_count <- Sc1a_site_list[which(Sc1a_site_list$Body==c("intron")),]
Sc1a_initial_count <- Sc1a_site_list[which(Sc1a_site_list$Body==c("initial")),]
Sc1a_internal_count <- Sc1a_site_list[which(Sc1a_site_list$Body==c("internal")),]
Sc1a_terminal_count <- Sc1a_site_list[which(Sc1a_site_list$Body==c("terminal")),]
Sc1a_single_count <- Sc1a_site_list[which(Sc1a_site_list$Body==c("single")),]

Sc1b_gene_count <- Sc1b_site_list[which(Sc1b_site_list$Body==c("transcript")),]
Sc1b_intron_count <- Sc1b_site_list[which(Sc1b_site_list$Body==c("intron")),]
Sc1b_initial_count <- Sc1b_site_list[which(Sc1b_site_list$Body==c("initial")),]
Sc1b_internal_count <- Sc1b_site_list[which(Sc1b_site_list$Body==c("internal")),]
Sc1b_terminal_count <- Sc1b_site_list[which(Sc1b_site_list$Body==c("terminal")),]
Sc1b_single_count <- Sc1b_site_list[which(Sc1b_site_list$Body==c("single")),]

Sc1c_gene_count <- Sc1c_site_list[which(Sc1c_site_list$Body==c("transcript")),]
Sc1c_intron_count <- Sc1c_site_list[which(Sc1c_site_list$Body==c("intron")),]
Sc1c_initial_count <- Sc1c_site_list[which(Sc1c_site_list$Body==c("initial")),]
Sc1c_internal_count <- Sc1c_site_list[which(Sc1c_site_list$Body==c("internal")),]
Sc1c_terminal_count <- Sc1c_site_list[which(Sc1c_site_list$Body==c("terminal")),]
Sc1c_single_count <- Sc1c_site_list[which(Sc1c_site_list$Body==c("single")),]

Sc1d_gene_count <- Sc1d_site_list[which(Sc1d_site_list$Body==c("transcript")),]
Sc1d_intron_count <- Sc1d_site_list[which(Sc1d_site_list$Body==c("intron")),]
Sc1d_initial_count <- Sc1d_site_list[which(Sc1d_site_list$Body==c("initial")),]
Sc1d_internal_count <- Sc1d_site_list[which(Sc1d_site_list$Body==c("internal")),]
Sc1d_terminal_count <- Sc1d_site_list[which(Sc1d_site_list$Body==c("terminal")),]
Sc1d_single_count <- Sc1d_site_list[which(Sc1d_site_list$Body==c("single")),]

Sc1e_gene_count <- Sc1e_site_list[which(Sc1e_site_list$Body==c("transcript")),]
Sc1e_intron_count <- Sc1e_site_list[which(Sc1e_site_list$Body==c("intron")),]
Sc1e_initial_count <- Sc1e_site_list[which(Sc1e_site_list$Body==c("initial")),]
Sc1e_internal_count <- Sc1e_site_list[which(Sc1e_site_list$Body==c("internal")),]
Sc1e_terminal_count <- Sc1e_site_list[which(Sc1e_site_list$Body==c("terminal")),]
Sc1e_single_count <- Sc1e_site_list[which(Sc1e_site_list$Body==c("single")),]

Sc1f_gene_count <- Sc1f_site_list[which(Sc1f_site_list$Body==c("transcript")),]
Sc1f_intron_count <- Sc1f_site_list[which(Sc1f_site_list$Body==c("intron")),]
Sc1f_initial_count <- Sc1f_site_list[which(Sc1f_site_list$Body==c("initial")),]
Sc1f_internal_count <- Sc1f_site_list[which(Sc1f_site_list$Body==c("internal")),]
Sc1f_terminal_count <- Sc1f_site_list[which(Sc1f_site_list$Body==c("terminal")),]
Sc1f_single_count <- Sc1f_site_list[which(Sc1f_site_list$Body==c("single")),]

Sc2_gene_count <- Sc2_site_list[which(Sc2_site_list$Body==c("transcript")),]
Sc2_intron_count <- Sc2_site_list[which(Sc2_site_list$Body==c("intron")),]
Sc2_initial_count <- Sc2_site_list[which(Sc2_site_list$Body==c("initial")),]
Sc2_internal_count <- Sc2_site_list[which(Sc2_site_list$Body==c("internal")),]
Sc2_terminal_count <- Sc2_site_list[which(Sc2_site_list$Body==c("terminal")),]
Sc2_single_count <- Sc2_site_list[which(Sc2_site_list$Body==c("single")),]

Sc2a_gene_count <- Sc2a_site_list[which(Sc2a_site_list$Body==c("transcript")),]
Sc2a_intron_count <- Sc2a_site_list[which(Sc2a_site_list$Body==c("intron")),]
Sc2a_initial_count <- Sc2a_site_list[which(Sc2a_site_list$Body==c("initial")),]
Sc2a_internal_count <- Sc2a_site_list[which(Sc2a_site_list$Body==c("internal")),]
Sc2a_terminal_count <- Sc2a_site_list[which(Sc2a_site_list$Body==c("terminal")),]
Sc2a_single_count <- Sc2a_site_list[which(Sc2a_site_list$Body==c("single")),]

Sc2b_gene_count <- Sc2b_site_list[which(Sc2b_site_list$Body==c("transcript")),]
Sc2b_intron_count <- Sc2b_site_list[which(Sc2b_site_list$Body==c("intron")),]
Sc2b_initial_count <- Sc2b_site_list[which(Sc2b_site_list$Body==c("initial")),]
Sc2b_internal_count <- Sc2b_site_list[which(Sc2b_site_list$Body==c("internal")),]
Sc2b_terminal_count <- Sc2b_site_list[which(Sc2b_site_list$Body==c("terminal")),]
Sc2b_single_count <- Sc2b_site_list[which(Sc2b_site_list$Body==c("single")),]

Sc2d_gene_count <- Sc2d_site_list[which(Sc2d_site_list$Body==c("transcript")),]
Sc2d_intron_count <- Sc2d_site_list[which(Sc2d_site_list$Body==c("intron")),]
Sc2d_initial_count <- Sc2d_site_list[which(Sc2d_site_list$Body==c("initial")),]
Sc2d_internal_count <- Sc2d_site_list[which(Sc2d_site_list$Body==c("internal")),]
Sc2d_terminal_count <- Sc2d_site_list[which(Sc2d_site_list$Body==c("terminal")),]
Sc2d_single_count <- Sc2d_site_list[which(Sc2d_site_list$Body==c("single")),]

Sc2e_gene_count <- Sc2e_site_list[which(Sc2e_site_list$Body==c("transcript")),]
Sc2e_intron_count <- Sc2e_site_list[which(Sc2e_site_list$Body==c("intron")),]
Sc2e_initial_count <- Sc2e_site_list[which(Sc2e_site_list$Body==c("initial")),]
Sc2e_internal_count <- Sc2e_site_list[which(Sc2e_site_list$Body==c("internal")),]
Sc2e_terminal_count <- Sc2e_site_list[which(Sc2e_site_list$Body==c("terminal")),]
Sc2e_single_count <- Sc2e_site_list[which(Sc2e_site_list$Body==c("single")),]

Z12_gene_count <- Z12_site_list[which(Z12_site_list$Body==c("transcript")),]
Z12_intron_count <- Z12_site_list[which(Z12_site_list$Body==c("intron")),]
Z12_initial_count <- Z12_site_list[which(Z12_site_list$Body==c("initial")),]
Z12_internal_count <- Z12_site_list[which(Z12_site_list$Body==c("internal")),]
Z12_terminal_count <- Z12_site_list[which(Z12_site_list$Body==c("terminal")),]
Z12_single_count <- Z12_site_list[which(Z12_site_list$Body==c("single")),]

Z15_gene_count <- Z15_site_list[which(Z15_site_list$Body==c("transcript")),]
Z15_intron_count <- Z15_site_list[which(Z15_site_list$Body==c("intron")),]
Z15_initial_count <- Z15_site_list[which(Z15_site_list$Body==c("initial")),]
Z15_internal_count <- Z15_site_list[which(Z15_site_list$Body==c("internal")),]
Z15_terminal_count <- Z15_site_list[which(Z15_site_list$Body==c("terminal")),]
Z15_single_count <- Z15_site_list[which(Z15_site_list$Body==c("single")),]

Z16_gene_count <- Z16_site_list[which(Z16_site_list$Body==c("transcript")),]
Z16_intron_count <- Z16_site_list[which(Z16_site_list$Body==c("intron")),]
Z16_initial_count <- Z16_site_list[which(Z16_site_list$Body==c("initial")),]
Z16_internal_count <- Z16_site_list[which(Z16_site_list$Body==c("internal")),]
Z16_terminal_count <- Z16_site_list[which(Z16_site_list$Body==c("terminal")),]
Z16_single_count <- Z16_site_list[which(Z16_site_list$Body==c("single")),]

Z17_gene_count <- Z17_site_list[which(Z17_site_list$Body==c("transcript")),]
Z17_intron_count <- Z17_site_list[which(Z17_site_list$Body==c("intron")),]
Z17_initial_count <- Z17_site_list[which(Z17_site_list$Body==c("initial")),]
Z17_internal_count <- Z17_site_list[which(Z17_site_list$Body==c("internal")),]
Z17_terminal_count <- Z17_site_list[which(Z17_site_list$Body==c("terminal")),]
Z17_single_count <- Z17_site_list[which(Z17_site_list$Body==c("single")),]

Z18_gene_count <- Z18_site_list[which(Z18_site_list$Body==c("transcript")),]
Z18_intron_count <- Z18_site_list[which(Z18_site_list$Body==c("intron")),]
Z18_initial_count <- Z18_site_list[which(Z18_site_list$Body==c("initial")),]
Z18_internal_count <- Z18_site_list[which(Z18_site_list$Body==c("internal")),]
Z18_terminal_count <- Z18_site_list[which(Z18_site_list$Body==c("terminal")),]
Z18_single_count <- Z18_site_list[which(Z18_site_list$Body==c("single")),]

Z19_gene_count <- Z19_site_list[which(Z19_site_list$Body==c("transcript")),]
Z19_intron_count <- Z19_site_list[which(Z19_site_list$Body==c("intron")),]
Z19_initial_count <- Z19_site_list[which(Z19_site_list$Body==c("initial")),]
Z19_internal_count <- Z19_site_list[which(Z19_site_list$Body==c("internal")),]
Z19_terminal_count <- Z19_site_list[which(Z19_site_list$Body==c("terminal")),]
Z19_single_count <- Z19_site_list[which(Z19_site_list$Body==c("single")),]

Z41_gene_count <- Z41_site_list[which(Z41_site_list$Body==c("transcript")),]
Z41_intron_count <- Z41_site_list[which(Z41_site_list$Body==c("intron")),]
Z41_initial_count <- Z41_site_list[which(Z41_site_list$Body==c("initial")),]
Z41_internal_count <- Z41_site_list[which(Z41_site_list$Body==c("internal")),]
Z41_terminal_count <- Z41_site_list[which(Z41_site_list$Body==c("terminal")),]
Z41_single_count <- Z41_site_list[which(Z41_site_list$Body==c("single")),]

Z43_gene_count <- Z43_site_list[which(Z43_site_list$Body==c("transcript")),]
Z43_intron_count <- Z43_site_list[which(Z43_site_list$Body==c("intron")),]
Z43_initial_count <- Z43_site_list[which(Z43_site_list$Body==c("initial")),]
Z43_internal_count <- Z43_site_list[which(Z43_site_list$Body==c("internal")),]
Z43_terminal_count <- Z43_site_list[which(Z43_site_list$Body==c("terminal")),]
Z43_single_count <- Z43_site_list[which(Z43_site_list$Body==c("single")),]

RAD_sig_loci_feat <- "RAD_loci_position_identifiers_mono_di.txt"
Significant_loci_feat <- read.table(file=RAD_sig_loci_feat, sep = "", header = FALSE)
colnames(Significant_loci_feat) <- c("Isolate", "Parent", "Nuclear_genotype", "Region","Body", "Feature", "Added_feature","Loci")
Significant_loci_feat$Nuclear_genotype <- factor(Significant_loci_feat$Nuclear_genotype, levels = c("mono", "di"))
Significant_loci_feat$Region <- factor(Significant_loci_feat$Region, levels = c("noncoding", "coding"))
Significant_loci_feat$Added_feature <- factor(Significant_loci_feat$Added_feature, levels = c("noncoding","intron", "single","initial","internal", "terminal"))

loci_feat <- Significant_loci_feat[which(Significant_loci_feat$Region=="coding"),]

p <-ggplot(Significant_loci_feat, aes(x=Nuclear_genotype, y=log(Loci+1))) +
  geom_boxplot( aes(Nuclear_genotype)) +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2), alpha=0.75)+
  facet_wrap(~Added_feature)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


##### Frequency conguency
##### A1
A1_sig_sites <- A1_common_sites[which(rownames(A1_common_sites) %in% rownames(Z25_sig_sites)),]

fit_A1_depth_sd <- lm(AC_Mean ~AF_SD,data=A1_sig_sites)
summary(fit_A1_depth_sd)
fit_Z25_depth_sd <- lm(AC_Mean ~AF_SD,data=Z25_sig_sites)
summary(fit_Z25_depth_sd)

p <-ggplot() +
  geom_point(data=Z25_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black") +
  geom_point(data=A1_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "grey") +
  geom_smooth(data=Z25_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black", method = 'lm') +
  geom_smooth(data=A1_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "grey", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p

fit_A1_Z25_AF_Mean <- lm(A1_AF_Mean ~Z25_AF_Mean,data=A1_Z25_x2sites)
summary(fit_A1_Z25_AF_Mean)

p <-ggplot() +
  geom_point(data=A1_Z25_x2sites, aes(x=A1_AF_Mean, y=Z25_AF_Mean),color="black", alpha=0.25) +
  geom_smooth(data=A1_Z25_x2sites, aes(x=A1_AF_Mean, y=Z25_AF_Mean),color="black", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


plot(density(A1_Z25_x2sites$A1_AF_Mean,bw=0.1), col="grey", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
lines(density(A1_Z25_x2sites$Z25_AF_Mean,bw=0.1), col="black", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
abline(v=mean(A1_Z25_x2sites$A1_AF_Mean), col="grey",lty = 2)
abline(v=mean(A1_Z25_x2sites$Z25_AF_Mean), col="black",lty = 2)

#####B12
B12ssp1_sig_sites <- B12ssp1_common_sites[which(rownames(B12ssp1_common_sites) %in% rownames(B12ssp5_sig_sites)),]

fit_B12ssp1_depth_sd <- lm(AC_Mean ~AF_SD,data=B12ssp1_sig_sites)
summary(fit_B12ssp1_depth_sd)
cor(x=B12ssp1_sig_sites$AC_Mean, y= B12ssp1_sig_sites$AF_SD,method = "spearman")
fit_B12ssp4_depth_sd <- lm(AC_Mean ~AF_SD,data=B12ssp4_sig_sites)
summary(fit_B12ssp4_depth_sd)
cor(x=B12ssp4_sig_sites$AC_Mean, y= B12ssp4_sig_sites$AF_SD,method = "spearman")

fit_B12ssp1_depth_sd <- lm(AC_Mean ~AF_SD,data=B12ssp1_sig_sites)
summary(fit_B12ssp1_depth_sd)
cor(x=B12ssp1_sig_sites$AC_Mean, y= B12ssp1_sig_sites$AF_SD,method = "spearman")
fit_B12ssp5_depth_sd <- lm(AC_Mean ~AF_SD,data=B12ssp5_sig_sites)
summary(fit_B12ssp5_depth_sd)
cor(x=B12ssp5_sig_sites$AC_Mean, y= B12ssp5_sig_sites$AF_SD,method = "spearman")

p <-ggplot() +
  geom_point(data=B12ssp5_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black") +
  geom_point(data=B12ssp1_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "grey") +
  geom_smooth(data=B12ssp5_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black", method = 'lm') +
  geom_smooth(data=B12ssp1_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "white", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p

fit_B12ssp1_B12ssp5_AF_Mean <- lm(B12ssp1_AF_Mean ~B12ssp5_AF_Mean,data=B12ssp1_B12ssp5_x2sites)
summary(fit_B12ssp1_B12ssp5_AF_Mean)

p <-ggplot() +
  geom_point(data=B12ssp1_B12ssp5_x2sites, aes(x=B12ssp1_AF_Mean, y=B12ssp5_AF_Mean),color="black", alpha=0.25) +
  geom_smooth(data=B12ssp1_B12ssp5_x2sites, aes(x=B12ssp1_AF_Mean, y=B12ssp5_AF_Mean),color="black", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


plot(density(B12ssp1_B12ssp5_x2sites$B12ssp1_AF_Mean,bw=0.1), col="grey", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
lines(density(B12ssp1_B12ssp5_x2sites$B12ssp5_AF_Mean,bw=0.1), col="black", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
abline(v=mean(B12ssp1_B12ssp5_x2sites$B12ssp1_AF_Mean), col="grey",lty = 2)
abline(v=mean(B12ssp1_B12ssp5_x2sites$B12ssp5_AF_Mean), col="black",lty = 2)

#####C2
C2_sig_sites <- C2_common_sites[which(rownames(C2_common_sites) %in% rownames(C2ssp2_sig_sites)),]

fit_C2_depth_sd <- lm(AC_Mean ~AF_SD,data=C2_sig_sites)
summary(fit_C2_depth_sd)
cor(x=C2_sig_sites$AC_Mean, y= C2_sig_sites$AF_SD,method = "spearman")
fit_C2ssp2_depth_sd <- lm(AC_Mean ~AF_SD,data=C2ssp2_sig_sites)
summary(fit_C2ssp2_depth_sd)
cor(x=C2ssp2_sig_sites$AC_Mean, y= C2ssp2_sig_sites$AF_SD,method = "spearman")

p <-ggplot() +
  geom_point(data=C2ssp2_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black") +
  geom_point(data=C2_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "grey") +
  geom_smooth(data=C2ssp2_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black", method = 'lm') +
  geom_smooth(data=C2_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "white", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p

fit_C2_C2ssp2_AF_Mean <- lm(C2_AF_Mean ~C2ssp2_AF_Mean,data=C2_C2ssp2_x2sites)
summary(fit_C2_C2ssp2_AF_Mean)

p <-ggplot() +
  geom_point(data=C2_C2ssp2_x2sites, aes(x=C2_AF_Mean, y=C2ssp2_AF_Mean),color="black", alpha=0.25) +
  geom_smooth(data=C2_C2ssp2_x2sites, aes(x=C2_AF_Mean, y=C2ssp2_AF_Mean),color="black", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


plot(density(C2_C2ssp2_x2sites$C2_AF_Mean,bw=0.1), col="grey", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
lines(density(C2_C2ssp2_x2sites$C2ssp2_AF_Mean,bw=0.1), col="black", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
abline(v=mean(C2_C2ssp2_x2sites$C2_AF_Mean), col="grey",lty = 2)
abline(v=mean(C2_C2ssp2_x2sites$C2ssp2_AF_Mean), col="black",lty = 2)

#####A5
A5_sig_sites <- A5_common_sites[which(rownames(A5_common_sites) %in% rownames(Z68_sig_sites)),]

fit_A5_depth_sd <- lm(AC_Mean ~AF_SD,data=A5_sig_sites)
summary(fit_A5_depth_sd)
cor(x=A5_sig_sites$AC_Mean, y= A5_sig_sites$AF_SD,method = "spearman")
fit_Z65_depth_sd <- lm(AC_Mean ~AF_SD,data=Z65_sig_sites)
summary(fit_Z65_depth_sd)
cor(x=Z65_sig_sites$AC_Mean, y= Z65_sig_sites$AF_SD,method = "spearman")

fit_A5_depth_sd <- lm(AC_Mean ~AF_SD,data=A5_sig_sites)
summary(fit_A5_depth_sd)
cor(x=A5_sig_sites$AC_Mean, y= A5_sig_sites$AF_SD,method = "spearman")
fit_Z68_depth_sd <- lm(AC_Mean ~AF_SD,data=Z68_sig_sites)
summary(fit_Z68_depth_sd)
cor(x=Z68_sig_sites$AC_Mean, y= Z68_sig_sites$AF_SD,method = "spearman")


p <-ggplot() +
  geom_point(data=Z68_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black") +
  geom_point(data=A5_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "grey") +
  geom_smooth(data=Z68_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black", method = 'lm') +
  geom_smooth(data=A5_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "white", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p

fit_A5_Z68_AF_Mean <- lm(A5_AF_Mean ~Z68_AF_Mean,data=A5_Z68_x2sites)
summary(fit_A5_Z68_AF_Mean)

p <-ggplot() +
  geom_point(data=A5_Z68_x2sites, aes(x=A5_AF_Mean, y=Z68_AF_Mean),color="black", alpha=0.25) +
  geom_smooth(data=A5_Z68_x2sites, aes(x=A5_AF_Mean, y=Z68_AF_Mean),color="black", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


plot(density(A5_Z68_x2sites$A5_AF_Mean,bw=0.1), col="grey", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
lines(density(A5_Z68_x2sites$Z68_AF_Mean,bw=0.1), col="black", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
abline(v=mean(A5_Z68_x2sites$A5_AF_Mean), col="grey",lty = 2)
abline(v=mean(A5_Z68_x2sites$Z68_AF_Mean), col="black",lty = 2)

#####C3
C3_sig_sites <- C3_common_sites[which(rownames(C3_common_sites) %in% rownames(Z17_sig_sites)),]

fit_C3_depth_sd <- lm(AC_Mean ~AF_SD,data=C3_sig_sites)
summary(fit_C3_depth_sd)
cor(x=C3_sig_sites$AC_Mean, y= C3_sig_sites$AF_SD,method = "spearman")
fit_Z15_depth_sd <- lm(AC_Mean ~AF_SD,data=Z15_sig_sites)
summary(fit_Z15_depth_sd)
cor(x=Z15_sig_sites$AC_Mean, y= Z15_sig_sites$AF_SD,method = "spearman")

fit_C3_depth_sd <- lm(AC_Mean ~AF_SD,data=C3_sig_sites)
summary(fit_C3_depth_sd)
cor(x=C3_sig_sites$AC_Mean, y= C3_sig_sites$AF_SD,method = "spearman")
fit_Z17_depth_sd <- lm(AC_Mean ~AF_SD,data=Z17_sig_sites)
summary(fit_Z17_depth_sd)
cor(x=Z17_sig_sites$AC_Mean, y= Z17_sig_sites$AF_SD,method = "spearman")

p <-ggplot() +
  geom_point(data=Z17_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black") +
  geom_point(data=C3_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "grey") +
  geom_smooth(data=Z17_sig_sites, aes(x=AC_Mean, y=AF_SD),color="black", method = 'lm') +
  geom_smooth(data=C3_sig_sites, aes(x=AC_Mean, y=AF_SD),color= "white", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p

fit_C3_Z17_AF_Mean <- lm(C3_AF_Mean ~Z17_AF_Mean,data=C3_Z17_x2sites)
summary(fit_C3_Z17_AF_Mean)

p <-ggplot() +
  geom_point(data=C3_Z17_x2sites, aes(x=C3_AF_Mean, y=Z17_AF_Mean),color="black", alpha=0.25) +
  geom_smooth(data=C3_Z17_x2sites, aes(x=C3_AF_Mean, y=Z17_AF_Mean),color="black", method = 'lm') +
  #scale_y_continuous(limits=c(5,45))+
  #ylab("Drought susceptibility index (DSI)")+
  #geom_jitter(aes(color=factor(Parent)), shape=16, position=position_jitter(0.2))+
  #facet_wrap(~Treatment)+
  theme_bw()+
  theme(panel.border= element_blank(),
        strip.background=element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        legend.position='none',
        axis.line= element_line(colour="black"))
p


plot(density(C3_Z17_x2sites$C3_AF_Mean,bw=0.1), col="grey", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
lines(density(C3_Z17_x2sites$Z17_AF_Mean,bw=0.1), col="black", xlab="Frequency", ylab="Density", xlim=c(0, 1),ylim=c(0, 4))
abline(v=mean(C3_Z17_x2sites$C3_AF_Mean), col="grey",lty = 2)
abline(v=mean(C3_Z17_x2sites$Z17_AF_Mean), col="black",lty = 2)
























#####################################################     END     ############################################################


C2_C2ssp1_sig_RA <- site_test(C2_RA_COM,C2ssp1_RA_COM)
C2_C2ssp2_sig_RA <- site_test(C2_RA_COM,C2ssp2_RA_COM)
C2_C2ssp3_sig_RA <- site_test(C2_RA_COM,C2ssp3_RA_COM)
C2_X35_sig_RA <- site_test(C2_RA_COM,X35_RA_COM)
C2_X36_sig_RA <- site_test(C2_RA_COM,X36_RA_COM)
C2_X38_sig_RA <- site_test(C2_RA_COM,X38_RA_COM)
C2_X39_sig_RA <- site_test(C2_RA_COM,X39_RA_COM)
X35_X39_sig_RA <- site_test(X35_RA_COM,X39_RA_COM)



library("reshape2")
library("Rmisc")
A5_sig_RA_sites <- A5_common_sites[which(rownames(A5_common_sites) %in% rownames(common_A5_sites)),]
A5_A5p1_sig_RA_sites <- melt(A5_A5p1_sig_RA_sites[,c("Sca","Pos","A5_R1_RF","A5_R2_RF","A5_R3_RF","A5_R4_RF","A5_R5_RF")],id.vars = c(1,2))
af_sum <- summarySE(A5_A5p1_sig_RA_sites, measurevar = c("value"), groupvars = c("Sca","Pos"), na.rm = TRUE)
A5_A5p1_sig_RA_sites_plot <- merge(A5_A5p1_sig_RA_sites, af_sum, by = c("Sca","Pos"))
A5_A5p1_sig_RA_sites_plot$site <- paste0(A5_A5p1_sig_RA_sites_plot$Sca,"_",A5_A5p1_sig_RA_sites_plot$Pos)
A5_A5p1_sig_RA_sites_plot$variable <- c("A5")

A5p1_A5_sig_RA_sites <- A5ssp1_common_sites[which(rownames(A5ssp1_common_sites) %in% rownames(A5_A5ssp1_sig_RA)),]
A5p1_A5_sig_RA_sites <- melt(A5p1_A5_sig_RA_sites[,c("Sca","Pos","A5ssp1_R1_RF","A5ssp1_R2_RF","A5ssp1_R3_RF","A5ssp1_R4_RF","A5ssp1_R5_RF")],id.vars = c(1,2))
af_sum <- summarySE(A5p1_A5_sig_RA_sites, measurevar = c("value"), groupvars = c("Sca","Pos"), na.rm = TRUE)
A5p1_A5_sig_RA_sites_plot <- merge(A5p1_A5_sig_RA_sites, af_sum, by = c("Sca","Pos"))
A5p1_A5_sig_RA_sites_plot$site <- paste0(A5p1_A5_sig_RA_sites_plot$Sca,"_",A5p1_A5_sig_RA_sites_plot$Pos)
A5p1_A5_sig_RA_sites_plot$variable <- c("A5p1")

A5_A5p1_plot <- rbind(A5_A5p1_sig_RA_sites_plot, A5p1_A5_sig_RA_sites_plot)
#A5_X27_sig_RA_all <- rbind(A5_X27_sig_RA_sites_plot,X27_A5_sig_RA_sites_plot)
#A5_X27_sig_RA_all$site <- paste0(A5_X27_sig_RA_all$Sca, "_", A5_X27_sig_RA_all$Pos)
#A5_X27_sig_RA_all$site <- factor(A5_X27_sig_RA_all$site, levels = levels(A5_X27_sig_RA_all$site))

plot_sites <- sample(A5_A5p1_plot$site,50)

A5_A5p1_allele_frequency <- A5_A5p1_plot[which(A5_A5p1_plot$site %in% plot_sites),]

p<-ggplot(A5_A5p1_plot, aes(x=site, y=value.y, fill=variable)) +
  #geom_rect(aes(xmin=-Inf, xmax = Inf, ymin = 0.2, ymax = 0.8),fill = "lightgray", alpha = 0.05)+
  geom_bar(stat="identity",position = "dodge", width = 0.6)+
  geom_errorbar(aes( ymin =value.y-se, ymax =value.y+se),stat="identity", width = 0.6, position = "dodge")+
  geom_point( aes(x=site, y=value.x, fill=variable),position=position_jitterdodge(jitter.width=0.05, dodge.width = 0.6), size=0.25,color="white", alpha=0.5) +
  geom_point( aes(x=site, y=value.y, fill=variable),position=position_jitterdodge(jitter.width=0.0, dodge.width = 0.6),size=0.5,color="black", alpha=0.5) +
  #geom_jitter(position = position_jitter(0.2),size = 1, color="darkgray") +
  geom_hline(aes(yintercept = 0.75),linetype = "dashed", colour= 'black') +
  #geom_hline(aes(yintercept = 0.6),linetype = "dashed", colour= 'black') +
  geom_hline(aes(yintercept = 0.5),linetype = "dashed", colour= 'red') +
  #geom_hline(aes(yintercept = 0.4), linetype = "dashed", colour= 'black') +
  geom_hline(aes(yintercept = 0.25), linetype = "dashed", colour= 'gray') +
  scale_fill_manual(values = c("#00D4EA","#97F5FF"))+
  theme_bw()+
  theme(panel.border= element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
       # legend.position = "none",
        text = element_text(size=15),
        axis.line= element_line(colour="black"))
p <- p + ylab("Observed allele frequency matching reference" )
#p <- p + scale_color_manual(values=getPalette(colcount))

p


common_sig_RA <- Reduce(intersect, list(rownames(A5_A5ssp1_sig_RA), rownames(A5_A5ssp2_sig_RA)))
#rownames(common_RA) <- paste0(common_RA$Sca,"_" ,common_RA$Pos)
#names <- rownames(common_RA)




A5sig <- A5_common_sites[which(rownames(A5_common_sites) %in% rownames(A5_A5ssp1_sig_RA)),]
A5p1sig <- A5ssp1_common_sites[which(rownames(A5ssp1_common_sites) %in% rownames(A5_A5ssp1_sig_RA)),]

new <- NULL
X39_R1_RA <- X39_RA_COM[,11]
X39_R1_AA <- X39_RA_COM[,14]
X39allalleles <- c(X39_R1_RA,X39_R1_AA)

C2_R1_RA <- C2_RA_COM[,15]
C2_R1_AA <- C2_RA_COM[,18]
C2allalleles <- c(C2_R1_RA,C2_R1_AA)

A5p1_R1_RA <- A5ssp1_common_sites[which(rownames(A5ssp1_common_sites) %in% rownames(A5ssp1_common_sites)),(ncol(A5ssp1_common_sites)-5)]
#A5p1_R1_AA <- A5ssp1_common_sites[which(rownames(A5ssp1_common_sites) %in% rownames(A5ssp1_common_sites)),(ncol(A5ssp1_common_sites)-2)]
A5p1new <- c(A5p1_R1_RA)

A5_R1_RA <- A5_common_sites[which(rownames(A5_common_sites) %in% rownames(A5_A5ssp1_sig_RA)),(ncol(A5_common_sites)-5)]
#A5_R1_AA <- A5_common_sites[which(rownames(A5_common_sites) %in% rownames(A5_A5ssp1_sig_RA)),(ncol(A5_common_sites)-2)]
A5new <- c(A5_R1_RA)  

plot(density(A5_A5p1_x2sites$A5_RF_Mean,bw=0.1), col="darkblue",main="Significantly different alleles between A5 and A5ssp1", xlab="Frequency", ylab="Density", xlim=c(0, 1), ylim=c(0, 4))
lines(density(A5_A5p1_x2sites$A5ssp1_RF_Mean, bw=0.1), col="lightblue")
abline(v=0.5, col="red",lty=3)
abline(v=0.25, col="grey",lty=3)
abline(v=0.75, col="grey",lty=3)
