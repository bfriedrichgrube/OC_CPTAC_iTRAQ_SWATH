
setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Classification/')

library(MSnbase)
library(mclust)
library(RColorBrewer)
library(pheatmap)

## Load

load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_protein_level_filtered_20171101.Rdata') 
dim(SWATH_filtered)
[1] 1659  103

load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_DDA_protein_level_CDAP.Rdata')

## Impute

mean(SWATH_filtered, na.rm=T)
[1] -0.1921363
sd(SWATH_filtered, na.rm=T)
[1] 0.5570528
count_missing <- length(which(is.na(SWATH_filtered)))
count_missing
[1] 21914

new_mean <- mean(SWATH_filtered, na.rm=T)*9
new_sd <- sd(SWATH_filtered, na.rm=T)*0.48
set.seed(0)
background <- rnorm(count_missing, mean=new_mean, sd=new_sd)
data_hist <- hist(SWATH_filtered)
back_hist <- hist(background)
plot(data_hist, xlim=c(-8,8), col='lightblue', main='Distribution of data and distribution of background', xlab='log2 intensities')
plot(back_hist, xlim=c(-8,8), col='red', main='Distribution of data and distribution of background', xlab='log2 intensities', add=T)
SWATH_imputed <- SWATH_filtered
length(which(is.na(SWATH_imputed)))
[1] 21914
SWATH_imputed[which(is.na(SWATH_imputed))]<-background
length(which(is.na(SWATH_imputed)))
[1] 0





missing <- NULL
for(i in 1:nrow(DDA)){tmp <-length(which(is.na(DDA[i,])))
missing <- c(missing, tmp)
}

DDA_missing0 <- DDA[(missing==0),] 
dim(DDA_missing0)
[1] 4366  103

gene.ids <- intersect(rownames(SWATH_imputed), rownames(DDA_missing0))
length(gene.ids)
[1] 1599
SWATH_overlap <- SWATH_imputed[gene.ids,]
DDA_overlap <- DDA_missing0[gene.ids,]
DDA_overlap <- as.matrix(DDA_overlap)
dim(SWATH_overlap)
[1] 1599  103
dim(DDA_overlap)
[1] 1599  103

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

set.seed(0)
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

SWATH_class <- CPTAC_SWATH_cluster_final$classification
length(SWATH_class)
[1] 103
DDA_class <- CPTAC_DDA_cluster_final$classification
head(DDA_class)
X09.1664.01A.Unshared.Log.Ratio X13.1404.01A.Unshared.Log.Ratio 
                              1                               1 
X13.1409.01A.Unshared.Log.Ratio X13.1410.01A.Unshared.Log.Ratio 
                              2                               3 
X13.1482.01A.Unshared.Log.Ratio X13.1483.01A.Unshared.Log.Ratio 
                              3                               3 
head(SWATH_class)
  HRD_1  HRD_10 HRD_101 HRD_102 HRD_103  HRD_11 
      1       2       1       1       1       2 
tmp <- read.delim('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_Clustering/clustering_results_all_170214.csv', sep=',', header=T, as.is=T)
tmp <- tmp[,c(1:2,13:15)]
classification_result <- cbind(tmp, DDA_class)
tmp1 <- classification_result
rownames(tmp1)<- tmp1[,2]
tmp1 <- tmp1[sort(rownames(tmp1)),]
classification_result <- cbind(tmp1, SWATH_class)

head(classification_result)
                  TCGA_name sample_name1      mRNA_subgroup CPTAC_subgroup
HRD_1   TCGA-09-1664-01A-01        HRD_1 Data not available Differentiated
HRD_10  TCGA-13-1489-01A-01       HRD_10      Proliferative Differentiated
HRD_101 TCGA-61-2008-01A-02      HRD_101     Differentiated Differentiated
HRD_102 TCGA-61-2094-01A-01      HRD_102     Immunoreactive Immunoreactive
HRD_103 TCGA-61-2613-01A-01      HRD_103 Data not available        Stromal
HRD_11  TCGA-13-1492-01A-01       HRD_11      Proliferative  Proliferative
        HRD_status DDA_class SWATH_class
HRD_1          HRD         1           1
HRD_10         HRD         1           2
HRD_101        HRD         1           1
HRD_102        HRD         1           1
HRD_103        HRD         3           1
HRD_11         HRD         2           2
adjustedRandIndex(classification_result$DDA_class, classification_result$SWATH_class)
[1] 0.2015766
table(classification_result$DDA_class, classification_result$SWATH_class)  
   
     1  2  3
  1 10 34  8
  2  3 17  0
  3  4  4 23

clinical<- read.delim('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_Clustering/CPTAC_Clustering_SWATH_final/CPTAC_SWATH_clinical_cluster_final.csv', sep=',', dec='.', header=T, row.names=1)


dim(clinical)
[1] 103  23

bi_1 <- classification_result$SWATH_class==1
bi_2 <- classification_result$SWATH_class==2
bi_3 <- classification_result$SWATH_class==3
clinical_SWATH <- clinical[sort(rownames(clinical)),1:17]
clinical_SWATH <- cbind(clinical_SWATH, classification_result[,7], bi_1, bi_2, bi_3)

colnames(clinical_SWATH)[15]<-'HRD_status'
colnames(clinical_SWATH)[18:21]<- c('SWATH_class', '1.binary', '2.binary', '3.binary')

head(clinical_SWATH)
        bcr_patient_barcode age_at_diagnosis tumor_stage tumor_grade
HRD_1          TCGA-09-1664               37        IIIC          G1
HRD_10         TCGA-13-1489               70        IIIC          G2
HRD_101        TCGA-61-2008               40         IIC          G2
HRD_102        TCGA-61-2094               63        IIIC          G3
HRD_103        TCGA-61-2613               73        IIIC          G3
HRD_11         TCGA-13-1492               66        IIIC          G3
        tumor_residual_disease vital_status daystodeath_or_LFU Prog_or_Recur
HRD_1                  1-10 mm     DECEASED               2278       Unknown
HRD_10           Not available     DECEASED               2553           Yes
HRD_101 No Macroscopic disease       LIVING                931           Yes
HRD_102          Not available       LIVING               2182            No
HRD_103               11-20 mm     DECEASED                945           Yes
HRD_11           Not available     DECEASED               3805            No
        daystotumorprog.recur.LFU PlatinumStatus PlatinumFreeInterval..mos.
HRD_1               Not available  Not available                         NA
HRD_10                        806      Sensitive                      17.50
HRD_101                       817      Sensitive                      22.00
HRD_102                      2182      Sensitive                      66.20
HRD_103                       358  Not available                         NA
HRD_11                       3805      Sensitive                     117.67
                     ethnicity  race year_of_diagnosis HRD_status BioReplicate
HRD_1   NOT HISPANIC OR LATINO WHITE              2002        HRD            1
HRD_10  NOT HISPANIC OR LATINO WHITE              2002        HRD           10
HRD_101 NOT HISPANIC OR LATINO ASIAN              2007        HRD          101
HRD_102 NOT HISPANIC OR LATINO WHITE              2003        HRD          102
HRD_103        [Not Available] WHITE              1998        HRD          103
HRD_11  NOT HISPANIC OR LATINO WHITE              2002        HRD           11
              mRNA_subtype SWATH_class 1.binary 2.binary 3.binary
HRD_1   Data not available           1     TRUE    FALSE    FALSE
HRD_10       Proliferative           2    FALSE     TRUE    FALSE
HRD_101     Differentiated           1     TRUE    FALSE    FALSE
HRD_102     Immunoreactive           1     TRUE    FALSE    FALSE
HRD_103 Data not available           1     TRUE    FALSE    FALSE
HRD_11       Proliferative           2    FALSE     TRUE    FALSE


rownames(classification_result) <- classification_result[,1]
classification_result <- classification_result[sort(rownames(classification_result)),]
head(classification_result)
                              TCGA_name sample_name1      mRNA_subgroup
TCGA-09-1664-01A-01 TCGA-09-1664-01A-01        HRD_1 Data not available
TCGA-13-1404-01A-01 TCGA-13-1404-01A-01    non_HRD_2      Proliferative
TCGA-13-1409-01A-01 TCGA-13-1409-01A-01    non_HRD_3      Proliferative
TCGA-13-1410-01A-01 TCGA-13-1410-01A-01    non_HRD_4     Immunoreactive
TCGA-13-1482-01A-01 TCGA-13-1482-01A-01        HRD_5     Differentiated
TCGA-13-1483-01A-01 TCGA-13-1483-01A-01        HRD_6        Mesenchymal
                    CPTAC_subgroup HRD_status DDA_class SWATH_class
TCGA-09-1664-01A-01 Differentiated        HRD         1           1
TCGA-13-1404-01A-01 Immunoreactive    non_HRD         1           2
TCGA-13-1409-01A-01  Proliferative    non_HRD         2           1
TCGA-13-1410-01A-01    Mesenchymal    non_HRD         3           3
TCGA-13-1482-01A-01    Mesenchymal        HRD         3           3
TCGA-13-1483-01A-01    Mesenchymal        HRD         3           3


bi_1 <- classification_result$DDA_class==1
bi_2 <- classification_result$DDA_class==2
bi_3 <- classification_result$DDA_class==3
head(clinical)
          bcr_patient_barcode age_at_diagnosis tumor_stage tumor_grade
HRD_1            TCGA-09-1664               37        IIIC          G1
non_HRD_2        TCGA-13-1404               48        IIIC          G3
non_HRD_3        TCGA-13-1409               73        IIIC          G3
non_HRD_4        TCGA-13-1410               57          IV          G3
HRD_5            TCGA-13-1482               52          IV          G2
HRD_6            TCGA-13-1483               61        IIIC          G3
          tumor_residual_disease vital_status daystodeath_or_LFU Prog_or_Recur
HRD_1                    1-10 mm     DECEASED               2278       Unknown
non_HRD_2                1-10 mm       LIVING               2211           Yes
non_HRD_3                1-10 mm     DECEASED               1703           Yes
non_HRD_4          Not available       LIVING               2168           Yes
HRD_5                     >20 mm     DECEASED               1877           Yes
HRD_6                     >20 mm     DECEASED                894           Yes
          daystotumorprog.recur.LFU PlatinumStatus PlatinumFreeInterval..mos.
HRD_1                 Not available  Not available                         NA
non_HRD_2                      1581      Sensitive                      43.50
non_HRD_3                       737      Sensitive                      17.53
non_HRD_4                       319      Sensitive                       6.37
HRD_5                           602      Sensitive                      14.80
HRD_6                           296      Resistant                       1.90
                       ethnicity  race year_of_diagnosis Condition BioReplicate
HRD_1     NOT HISPANIC OR LATINO WHITE              2002       HRD            1
non_HRD_2 NOT HISPANIC OR LATINO WHITE              2008   non_HRD            2
non_HRD_3 NOT HISPANIC OR LATINO WHITE              2008   non_HRD            3
non_HRD_4 NOT HISPANIC OR LATINO WHITE              2008   non_HRD            4
HRD_5     NOT HISPANIC OR LATINO WHITE              2000       HRD            5
HRD_6     NOT HISPANIC OR LATINO WHITE              2000       HRD            6
                mRNA_subtype  CPTAC_subtype SWATH_class  bi_1  bi_2  bi_3  bi_4
HRD_1     Data not available Differentiated           1  TRUE FALSE FALSE FALSE
non_HRD_2      Proliferative Immunoreactive           2 FALSE  TRUE FALSE FALSE
non_HRD_3      Proliferative  Proliferative           1  TRUE FALSE FALSE FALSE
non_HRD_4     Immunoreactive    Mesenchymal           4 FALSE FALSE FALSE  TRUE
HRD_5         Differentiated    Mesenchymal           4 FALSE FALSE FALSE  TRUE
HRD_6            Mesenchymal    Mesenchymal           2 FALSE  TRUE FALSE FALSE

clinical_DDA <- clinical[,1:17]
clinical_DDA <- cbind(clinical_DDA, classification_result[,6], bi_1, bi_2, bi_3)
 
colnames(clinical_DDA)[15]<-'HRD_status'
colnames(clinical_DDA)[18:21]<- c('DDA_class', '1.binary', '2.binary', '3.binary')
save(classification_result, file='classification_result_20171101.Rdata')
# Achtung previous results were overwitten save(DDA_overlap, DDA_missing0, DDA_z, CPTAC_DDA_cluster_final, file='DDA_classification_workflow_20171027.Rdata')
save(DDA_overlap, DDA_missing0, DDA_z, CPTAC_DDA_cluster_final, file='DDA_classification_workflow_20171101.Rdata')
save(SWATH_overlap, SWATH_z, SWATH_imputed, SWATH_filtered, CPTAC_SWATH_cluster_final,  file='SWATH_classification_workflow_20171101.Rdata')
save(clinical_DDA, clinical_SWATH, file='clinical_annotation_tables_20161101.Rdata')
q()

