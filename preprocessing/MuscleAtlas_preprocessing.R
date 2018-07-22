setwd("C:/Dropbox/NICO/R/Shiny/Muscle_Atlas/preprocessing")
genename <- 'A4GNT'
muscle  <- read.table("C:/Dropbox/NICO/R/Muscle_Atlas/MuscleAtlas_GENENAME_norm.txt")
muscle <- muscle[order(colnames(muscle))]
colnames(muscle)
min(muscle, na.rm=T)
max(muscle, na.rm=T)


#list of gene names as RDS
list_genes <- rownames(muscle)
saveRDS(list_genes, '../data/MuscleAtlas_genelist.Rds')


#Make a list of data for ggplot
res <- muscle
x      <- c(rep('A1', length(grep('Blood',         colnames(muscle)))),  #list of sample types
            rep('A2', length(grep('Endothelial',   colnames(muscle)))),
            rep('A3', length(grep('Fibroblast',    colnames(muscle)))),
            rep('A4', length(grep('Macrophage',    colnames(muscle)))),
            rep('A5', length(grep('Muscle_Biopsy', colnames(muscle)))),
            rep('A6', length(grep('Muscle_Fiber',  colnames(muscle)))),
            rep('A7', length(grep('Myocyte',       colnames(muscle)))),
            rep('A8', length(grep('Smooth_muscle', colnames(muscle)))))
datalist <- vector("list", nrow(res))
names(datalist) <- rownames(res)
for (i in 1:nrow(res)){
  data   <- data.frame(x=factor(), y=numeric(), Gene=character(), stringsAsFactors=FALSE) #empty dataframe to collect data
  y     <- as.numeric(res[i,])            #collect data for gene name i
  data <- data.frame(x, y, rep(rownames(res[i,]))) #create table with x="sample type", y="data", "gene name"
  colnames(data) <- c("x","y","Gene")              #rename column names to make it possible to rbind later
  datalist[[i]] <- data
}
saveRDS(datalist, '../data/MuscleAtlas_data.Rds')


#stats for all genes
res <- muscle
library(matrixStats)
statslist <- vector("list", nrow(res))
names(statslist) <- rownames(res)

for (i in 1:nrow(res)){
mean <- cbind(
rowMeans(res[i, grepl('Blood',        colnames(res))], na.rm=T),
rowMeans(res[i, grepl('Endothelial',  colnames(res))], na.rm=T),
rowMeans(res[i, grepl('Fibroblast',   colnames(res))], na.rm=T),
rowMeans(res[i, grepl('Macrophage',   colnames(res))], na.rm=T),
rowMeans(res[i, grepl('Muscle_Biopsy',colnames(res))], na.rm=T),
rowMeans(res[i, grepl('Muscle_Fiber', colnames(res))], na.rm=T),
rowMeans(res[i, grepl('Myocyte',      colnames(res))], na.rm=T),
rowMeans(res[i, grepl('Smooth_muscle',colnames(res))], na.rm=T))
Sd <- cbind(
  rowSds(as.matrix(res[i, grepl('Blood',        colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('Endothelial',  colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('Fibroblast',   colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('Macrophage',   colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('Muscle_Biopsy',colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('Muscle_Fiber', colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('Myocyte',      colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('Smooth_muscle',colnames(res))]), na.rm=T))
nsize <- cbind(
  rowSums(!is.na(res[i, grepl('Blood',        colnames(res))])),
  rowSums(!is.na(res[i, grepl('Endothelial',  colnames(res))])),
  rowSums(!is.na(res[i, grepl('Fibroblast',   colnames(res))])),
  rowSums(!is.na(res[i, grepl('Macrophage',   colnames(res))])),
  rowSums(!is.na(res[i, grepl('Muscle_Biopsy',colnames(res))])),
  rowSums(!is.na(res[i, grepl('Muscle_Fiber', colnames(res))])),
  rowSums(!is.na(res[i, grepl('Myocyte',      colnames(res))])),
  rowSums(!is.na(res[i, grepl('Smooth_muscle',colnames(res))])))
stats <- data.frame(t(mean), t(Sd), t(nsize))
stats[,3] <- as.factor(stats[,3])
colnames(stats) <- c('Mean', 'Sd', 'n')
rownames(stats) <- c("Whole Blood", "Endothelial Cell", "Fibroblast", "Macrophage", 
                     "Muscle Biopsy", "Muscle Fiber", "Primary Myocyte", "Smooth Muscle Cell")
statslist[[i]] <- stats
}
saveRDS(statslist, '../data/MuscleAtlas_statslist.RDS')


