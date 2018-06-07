## Normalization and batch correction of mapDIA input ##

## Step 1: Prepare R session ##

library(preprocessCore)
library(pheatmap)
library(RColorBrewer)

setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_SWATH2stats_preprocessing/CPTAC_SWATH2stats_2stepfiltering_final/103files_newlib_jumbo_20170207_preprocess/')

sessionInfo()
R version 3.2.2 (2015-08-14)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 15.10

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] RColorBrewer_1.1-2    pheatmap_1.0.8        preprocessCore_1.32.0

loaded via a namespace (and not attached):
[1] colorspace_1.2-6 scales_0.4.0     plyr_1.8.3       gtable_0.2.0    
[5] Rcpp_0.12.7      grid_3.2.2       munsell_0.4.3   



## Step 2: Load in data ##

CPTAC_data <- read.delim2('CPTAC_103HRD_2stepfilter_20170126_mapDIA_input.tsv', dec='.', sep='\t', header=T)
dim(CPTAC_data)
[1] 150408    107


#log2 transform data
CPTAC_data_log <- log2(as.matrix(CPTAC_data[,4:106]))
CPTAC_data_log[which(is.infinite(CPTAC_data_log))] <- NA # replace the infinite values with NA


batch_annotation <- read.delim2('CPTAC_batch_annotation.csv', sep=',', header=T)
dim(batch_annotation)
[1] 103   3

annotation.file <- read.delim2("CPTAC_103HRD_annotation_20160929.txt", dec=".", sep="\t", header=T)
dim(annotation.file)
[1] 103   4

batch_info <- cbind(annotation.file, batch_annotation$batch)
head(batch_info)
                           Filename Condition BioReplicate Run
1 TCGA-09-1664-01A-01_S_W_JHU_141107       HRD            1   1
2 TCGA-13-1404-01A-01_S_W_JHU_141115   non_HRD            2   2
3 TCGA-13-1409-01A-01_S_W_JHU_141102   non_HRD            3   3
4 TCGA-13-1410-01A-01_S_W_JHU_141103   non_HRD            4   4
5 TCGA-13-1482-01A-01_S_W_JHU_141102       HRD            5   5
6 TCGA-13-1483-01A-01_S_W_JHU_141103       HRD            6   6
  batch_annotation$batch
1                      C
2                      E
3                      A
4                      A
5                      A
6                      A



batch_info <- batch_info[order(batch_info[,2]),]
head(batch_info)
                             Filename Condition BioReplicate Run
1  TCGA-09-1664-01A-01_S_W_JHU_141107       HRD            1   1
5  TCGA-13-1482-01A-01_S_W_JHU_141102       HRD            5   5
6  TCGA-13-1483-01A-01_S_W_JHU_141103       HRD            6   6
9  TCGA-13-1488-01A-01_S_W_JHU_141115       HRD            9   9
10 TCGA-13-1489-01A-01_S_W_JHU_141031       HRD           10  10
11 TCGA-13-1492-01A-01_S_W_JHU_141115       HRD           11  11
   batch_annotation$batch
1                       C
5                       A
6                       A
9                       E
10                      A
11                      E

## Step 3: Make heatmap with batch annotation above columns to make batch effect visible ##

anno_cols <- data.frame(
		batch = batch_info[,5]
		)
#dependent on the dataset that should be displayed, change the rownames accordingly SWATH sample names when swath data is displayed, 
rownames(anno_cols) <- colnames(CPTAC_data[,4:106])
#color_anno <- data.frame(
#		batch = c(A='blue', B='green', C='red', D='yellow', E='orange')
#		)		

#CPTAC_data_heat <- as.matrix(CPTAC_data[,4:106])

pheatmap(CPTAC_data_log, show_rownames=F, show_colnames=F, filename='heatmap_data_before_batchcorr.pdf', border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_cols)


## Step 4: Normalization ##

CPTAC_data_norm <- normalize.quantiles(CPTAC_data_log)
colnames(CPTAC_data_norm) <- colnames(CPTAC_data_log)

#check heatmap again
pheatmap(CPTAC_data_norm, show_rownames=F, show_colnames=F, filename='heatmap_data_norm_before_batchcorr.pdf', border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_cols)


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
[1] 150408    107


write.table(CPTAC_103HRD_mapDIA_input_norm_batchcorr, file="CPTAC_103HRD_jumbo_mapDIA_input_norm_batchcorr_20170208.tsv", quote=F, row.names=F, sep="\t")

save(batch_info, CPTAC_103HRD_mapDIA_input_norm_batchcorr, CPTAC_data, CPTAC_data_log, CPTAC_data_norm, CPTAC_data_batchcorr, CPTAC_data_batchcorr_order, file='CPTAC_103HRD_quantnorm_batchcorr_20170208.Rdata')

length(unique(CPTAC_103HRD_mapDIA_input_norm_batchcorr$ProteinName))
[1] 4291

## feed into mapDIA without normalisation, without log2-transformation

