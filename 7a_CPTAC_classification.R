## Mclust classification algorithm (For Figure 3, S2)

## Step 1: Prepare R session

#load packages
library(MSnbase)
library(mclust)

# set working directory
setwd('')

## Step 2: Load required data

#load protein data
load('CPTAC_SWATH_overlap.Rdata') 
load('CPTAC_DDA_overlap.Rdata')

#read in clinical and annotation table

clinical<- read.delim('CPTAC_clinical.csv', sep=',', dec='.', header=T, row.names=1)
dim(clinical)
[1] 103  23

## Step 3: Calculate z-scores for clustering

SWATH_z <- SWATH_overlap
SWATH_z <- sweep(SWATH_z, 1, apply(SWATH_z, 1, mean, na.rm=T), "-")
SWATH_z <- sweep(SWATH_z, 1, apply(SWATH_z, 1, sd, na.rm=T), "/")
dim(SWATH_z)
[1] 1599  103

DDA_z <- DDA_overlap
DDA_z <- sweep(DDA_z, 1, apply(DDA_z, 1, mean, na.rm=T), "-")
DDA_z <- sweep(DDA_z, 1, apply(DDA_z, 1, sd, na.rm=T), "/")
dim(DDA_z)
[1] 1599  103

DDA_z0 <- DDA_missing0
DDA_z0 <- sweep(DDA_z0, 1, apply(DDA_z0, 1, mean, na.rm=T), "-")
DDA_z0 <- sweep(DDA_z0, 1, apply(DDA_z0, 1, sd, na.rm=T), "/")
dim(DDA_z0)
[1] 4366  103


## Step 4: Classification

set.seed(0)

#SWATH common proteins
CPTAC_SWATH_cluster_final <- Mclust(t(SWATH_z))

summary(CPTAC_SWATH_cluster_final)
----------------------------------------------------
Gaussian finite mixture model fitted by EM algorithm 
----------------------------------------------------

Mclust VII (spherical, varying volume) model with 3 components:

 log.likelihood   n   df       BIC       ICL
      -212801.5 103 4802 -447858.9 -447858.9

Clustering table:
 1  2  3 
17 55 31 

SWATH_class <- CPTAC_SWATH_cluster_final$classification


#iTRAQ DDA common proteins
CPTAC_DDA_cluster_final <- Mclust(t(DDA_z))

summary(CPTAC_DDA_cluster_final)
----------------------------------------------------
Gaussian finite mixture model fitted by EM algorithm 
----------------------------------------------------

Mclust VII (spherical, varying volume) model with 3 components:

 log.likelihood   n   df       BIC       ICL
      -220097.7 103 4802 -462451.5 -462451.5

Clustering table:
 1  2  3 
52 20 31 

DDA_class <- CPTAC_DDA_cluster_final$classification


# iTRAQ DDA all without missing values
CPTAC_DDA_cluster_all <- Mclust(t(DDA_z0))

summary(CPTAC_DDA_cluster_all)
----------------------------------------------------
Gaussian finite mixture model fitted by EM algorithm 
----------------------------------------------------

Mclust VII (spherical, varying volume) model with 3 components:

 log.likelihood   n    df      BIC      ICL
      -604107.5 103 13103 -1268944 -1268944

Clustering table:
 1  2  3 
56 16 31 

DDA_class_all <- CPTAC_DDA_cluster_all$classification

#combine results into a table and calculate similarity
#update and upload table
tmp <- read.delim('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_Clustering/clustering_results_all_170214.csv', sep=',', header=T, as.is=T)
tmp <- tmp[,c(1:2,13:15)]
classification_result <- cbind(tmp, DDA_class, DDA_class_all)
tmp1 <- classification_result
rownames(tmp1)<- tmp1[,2]
tmp1 <- tmp1[sort(rownames(tmp1)),]
classification_result <- cbind(tmp1, SWATH_class)

rownames(classification_result) <- classification_result[,1]
classification_result <- classification_result[sort(rownames(classification_result)),]
head(classification_result)

save(DDA_z, CPTAC_DDA_cluster_final, file='DDA_classification_workflow.Rdata')
save(SWATH_z, CPTAC_SWATH_cluster_final,  file='SWATH_classification_workflow.Rdata')
save(classification_result, file='classification_result.Rdata')

## Step 5: Calculate Adjusted Rand Index (Figure 3B)

#ARI of results from common proteins
adjustedRandIndex(classification_result$DDA_class, classification_result$SWATH_class)
[1] 0.2146062

#ARI of results from SWATH common proteins vs. iTRAQ DDA original
adjustedRandIndex(classification_result$CPTAC_subgroup, classification_result$SWATH_class)
[1] 0.1400079

#ARI of results from iTRAQ DDA all vs. original
adjustedRandIndex(classification_result$CPTAC_subgroup, DDA_class_all)
[1] 0.3641888

#ARI of results from iTRAQ DDA common proteins vs. original
adjustedRandIndex(classification_result$CPTAC_subgroup, classification_result$DDA_class)
[1] 0.329218

#ARI of results from iTRAQ DDA all vs. common proteins
adjustedRandIndex(classification_result$DDA_class, DDA_class_all)
[1] 0.5831374

#ARI of results from iTRAQ DDA original vs. TCGA subgroup
adjustedRandIndex(classification_result$CPTAC_subgroup, classification_result$mRNA_subgroup)
[1] 0.201455

## Step 5: Prepare MSnSet Data for following step

#Prepare SWATH MSnSet
#add binary information of classification
bi_1 <- classification_result$SWATH_class==1
bi_2 <- classification_result$SWATH_class==2
bi_3 <- classification_result$SWATH_class==3
clinical_SWATH <- clinical[sort(rownames(clinical)),1:17]
clinical_SWATH <- cbind(clinical_SWATH, classification_result[,7], bi_1, bi_2, bi_3)
colnames(clinical_SWATH)[15]<-'HRD_status'
colnames(clinical_SWATH)[18:21]<- c('SWATH_class', '1.binary', '2.binary', '3.binary')

all(rownames(clinical_SWATH)==colnames(SWATH_overlap))
[1] TRUE

clinical_SWATH_anno <- new('AnnotatedDataFrame', data=clinical_SWATH)
CPTAC_SWATH_Set <- new('MSnSet', exprs=SWATH_overlap, phenoData=clinical_SWATH_anno)

#Prepare DDA MSnSet 
#add binary information of classification
bi_1 <- classification_result$DDA_class==1
bi_2 <- classification_result$DDA_class==2
bi_3 <- classification_result$DDA_class==3
clinical_DDA <- clinical[,1:17]
clinical_DDA <- cbind(clinical_DDA, classification_result[,6], bi_1, bi_2, bi_3)
colnames(clinical_DDA)[15]<-'HRD_status'
colnames(clinical_DDA)[18:21]<- c('DDA_class', '1.binary', '2.binary', '3.binary')

colnames(DDA_overlap)<- rownames(clinical_DDA)
all(rownames(clinical_DDA)==colnames(DDA_overlap))
[1] TRUE

clinical_DDA_anno <- new('AnnotatedDataFrame', data=clinical_DDA)
CPTAC_DDA_Set <- new('MSnSet', exprs=DDA_overlap, phenoData=clinical_DDA_anno)


save(CPTAC_SWATH_Set, file='CPTAC_SWATH_Set.Rdata')
save(CPTAC_DDA_Set, file='CPTAC_DDA_Set.Rdata')



