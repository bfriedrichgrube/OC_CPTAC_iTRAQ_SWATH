## After preprocessing, the output needs to be batch corrected and normalized ##

## Step 1: Prepare R session ##

library(preprocessCore)
library(pheatmap)
library(RColorBrewer)

setwd('')

sessionInfo()

## Step 2: Load in data ##

CPTAC_data <- read.delim2('mapDIA_input.tsv', dec='.', sep='\t', header=T)
dim(CPTAC_data)

#log2 transform data
CPTAC_data_log <- log2(as.matrix(CPTAC_data[,4:106]))
CPTAC_data_log[which(is.infinite(CPTAC_data_log))] <- NA # replace the infinite values with NA


batch_annotation <- read.delim2('batch_annotation.csv', sep=',', header=T)
dim(batch_annotation)

annotation.file <- read.delim2("annotation.txt", dec=".", sep="\t", header=T)
dim(annotation.file)

batch_info <- cbind(annotation.file, batch_annotation$batch)
head(batch_info)

batch_info <- batch_info[order(batch_info[,2]),]
head(batch_info)

## Step 3: Make heatmap with batch annotation above columns to make batch effect visible ##

anno_cols <- data.frame(
		batch = batch_info[,5]
		)
#dependent on the dataset that should be displayed, change the rownames accordingly SWATH sample names when swath data is displayed, 
rownames(anno_cols) <- colnames(CPTAC_data[,4:106])
#color_anno <- data.frame(
#		batch = c(A='blue', B='green', C='red', D='yellow', E='orange')
#		)		

pheatmap(CPTAC_data_log, show_rownames=F, show_colnames=F, filename='heatmap_data_before_batchcorr.pdf', border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_cols)


## Step 4: Normalization ##

CPTAC_data_norm <- normalize.quantiles(CPTAC_data_log)
colnames(CPTAC_data_norm) <- colnames(CPTAC_data_log)

## Step 5: Batch correction ##
CPTAC_data_norm <- cbind(CPTAC_data[,1:3], CPTAC_data_norm, CPTAC_data[,107])
colnames(CPTAC_data_norm)[107] <- 'RT'


batch_info[,c(2,3,5)]
batch_info[order(batch_info[,5]),4:5]


data_batchA <- CPTAC_data_norm[,c(5, 6, 8, 11, 20, 24, 36, 39, 40, 48, 49, 58, 63, 64, 70, 74, 76, 77, 81, 82, 84, 85, 102, 103)]
data_batchB <- CPTAC_data_norm[,c(16, 18, 19, 26, 27, 29, 31:33, 35, 46, 50:52, 54, 56, 61, 65, 68, 90, 92, 94, 97, 99)]
data_batchC <- CPTAC_data_norm[,c(4, 10, 12, 13, 15, 17, 22, 28, 30, 34, 38, 43:45, 53, 66, 73, 78, 86, 88, 93, 96, 98, 100, 106)]
data_batchD <- CPTAC_data_norm[,c(14, 21, 23, 25, 37, 60, 67, 71, 72, 75, 79, 80, 83, 89, 91, 105)]
data_batchE <- CPTAC_data_norm[,c(7, 9, 41, 42, 47, 55, 57, 59, 62, 69, 87, 95, 101, 104)]

data_batchA_corr <- as.matrix(data_batchA) - rowMeans(as.matrix(data_batchA), na.rm=T)
data_batchB_corr <- as.matrix(data_batchB) - rowMeans(as.matrix(data_batchB), na.rm=T)
data_batchC_corr <- as.matrix(data_batchC) - rowMeans(as.matrix(data_batchC), na.rm=T)
data_batchD_corr <- as.matrix(data_batchD) - rowMeans(as.matrix(data_batchD), na.rm=T)
data_batchE_corr <- as.matrix(data_batchE) - rowMeans(as.matrix(data_batchE), na.rm=T)
CPTAC_data_batchcorr <- cbind(data_batchA_corr, data_batchB_corr, data_batchC_corr, data_batchD_corr, data_batchE_corr)

#check heatmap again 
pheatmap(CPTAC_data_batchcorr, show_rownames=F, show_colnames=F, filename='heatmap_data_norm_after_batchcorr.pdf', border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_cols)

## Step 5: save for mapDIA 
#order for mapDIA input
colnames(CPTAC_data_batchcorr)
CPTAC_data_batchcorr_order <- CPTAC_data_batchcorr[,order(colnames(CPTAC_data_batchcorr))]

# adjustment for rowMean substraction
CPTAC_data_norm_mat <- as.matrix(CPTAC_data_norm[,4:106])
rm <- rowMeans(CPTAC_data_norm_mat,na.rm=T)
CPTAC_data_batchcorr_order_final <- CPTAC_data_batchcorr_order + rm

CPTAC_103HRD_mapDIA_input_norm_batchcorr <- cbind(CPTAC_data[,1:3], CPTAC_data_batchcorr_order_final, CPTAC_data[,107])
colnames(CPTAC_103HRD_mapDIA_input_norm_batchcorr)[107] <- 'RT'
dim(CPTAC_103HRD_mapDIA_input_norm_batchcorr)

write.table(CPTAC_103HRD_mapDIA_input_norm_batchcorr, file="mapDIA_input_norm_batchcorr.tsv", quote=F, row.names=F, sep="\t")

save(batch_info, CPTAC_103HRD_mapDIA_input_norm_batchcorr, CPTAC_data, CPTAC_data_log, CPTAC_data_norm, CPTAC_data_batchcorr, CPTAC_data_batchcorr_order, file='CPTAC_103HRD_quantnorm_batchcorr_20170208.Rdata')

length(unique(CPTAC_103HRD_mapDIA_input_norm_batchcorr$ProteinName))
[1] 4291

## feed into mapDIA without normalisation, without log2-transformation
## next step is documented in parameter file of mapDIA
