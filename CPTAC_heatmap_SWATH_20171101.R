
setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Classification/')
load('CPTAC_SWATH_final_Set_20171101.Rdata')
library(RColorBrewer)
library(pheatmap)
library(MSnbase)
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:parallel’:

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from ‘package:stats’:

    IQR, mad, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, append, as.data.frame, as.vector, cbind, colnames,
    do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl,
    intersect, is.unsorted, lapply, lengths, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unlist, unsplit

Loading required package: Biobase
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.

Loading required package: mzR
Loading required package: Rcpp
Loading required package: BiocParallel
Loading required package: ProtGenerics

This is MSnbase version 1.18.1 
  Read '?MSnbase' and references therein for information
  about the package and how to get started.


Attaching package: ‘MSnbase’

The following object is masked from ‘package:stats’:

    smooth
load('CPTAC_DDA_final_Set_20171101.Rdata')
load('CPTAC_SWATH_WGCNA_20171101.Rdata')
> ls()
[1] "clustEG"         "CPTAC_DDA_Set"   "CPTAC_SWATH_Set"
> summary(clustEG)
          Length Class  Mode     
black       86   -none- character
green       42   -none- character
grey       360   -none- character
red         41   -none- character
turquoise 1063   -none- character

heatmap_black <- exprs(CPTAC_SWATH_Set)[names(clustEG$black),]
heatmap_green <- exprs(CPTAC_SWATH_Set)[names(clustEG$green),]
heatmap_grey <- exprs(CPTAC_SWATH_Set)[names(clustEG$grey),]
heatmap_grey <- exprs(CPTAC_SWATH_Set)[names(clustEG$grey),]
heatmap_red <- exprs(CPTAC_SWATH_Set)[names(clustEG$red),]
heatmap_turquoise <- exprs(CPTAC_SWATH_Set)[names(clustEG$turquoise),]
heatmap_input <- rbind(heatmap_black, heatmap_green, heatmap_grey, heatmap_red, heatmap_turquoise)
SWATH_class <- pData(CPTAC_SWATH_Set)[,18]
names(SWATH_class)<- rownames(pData(CPTAC_SWATH_Set))
SWATH_class <- sort(SWATH_class)
heatmap_input <- as.matrix(heatmap_input[,names(SWATH_class)])
DDA_class <- pData(CPTAC_DDA_Set)[,18]

names(DDA_class) <- rownames(pData(CPTAC_DDA_Set))
head(SWATH_class)
  HRD_1 HRD_101 HRD_102 HRD_103  HRD_29  HRD_34 
      1       1       1       1       1       1 

head(DDA_class)
    HRD_1 non_HRD_2 non_HRD_3 non_HRD_4     HRD_5     HRD_6 
        1         1         2         3         3         3 
head(pData(CPTAC_DDA_Set))
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
                       ethnicity  race year_of_diagnosis HRD_status
HRD_1     NOT HISPANIC OR LATINO WHITE              2002        HRD
non_HRD_2 NOT HISPANIC OR LATINO WHITE              2008    non_HRD
non_HRD_3 NOT HISPANIC OR LATINO WHITE              2008    non_HRD
non_HRD_4 NOT HISPANIC OR LATINO WHITE              2008    non_HRD
HRD_5     NOT HISPANIC OR LATINO WHITE              2000        HRD
HRD_6     NOT HISPANIC OR LATINO WHITE              2000        HRD
          BioReplicate       mRNA_subtype DDA_class 1.binary 2.binary 3.binary
HRD_1                1 Data not available         1     TRUE    FALSE    FALSE
non_HRD_2            2      Proliferative         1     TRUE    FALSE    FALSE
non_HRD_3            3      Proliferative         2    FALSE     TRUE    FALSE
non_HRD_4            4     Immunoreactive         3    FALSE    FALSE     TRUE
HRD_5                5     Differentiated         3    FALSE    FALSE     TRUE
HRD_6                6        Mesenchymal         3    FALSE    FALSE     TRUE
> head(pData(CPTAC_SWATH_Set))
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
tmp <- pData(CPTAC_DDA_Set)
dim(tmp)
[1] 103  21

tmp <- tmp[sort(rownames(tmp)),]
head(tmp)
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
              mRNA_subtype DDA_class 1.binary 2.binary 3.binary
HRD_1   Data not available         1     TRUE    FALSE    FALSE
HRD_10       Proliferative         1     TRUE    FALSE    FALSE
HRD_101     Differentiated         1     TRUE    FALSE    FALSE
HRD_102     Immunoreactive         1     TRUE    FALSE    FALSE
HRD_103 Data not available         3    FALSE    FALSE     TRUE
HRD_11       Proliferative         2    FALSE     TRUE    FALSE


anno_cols_all <- data.frame(
	SWATH_class  = pData(CPTAC_SWATH_Set)[,18],
	DDA_class = tmp[,18],
	mRNA_class = tmp[,17])
rownames(anno_cols_all) <- rownames(tmp)
gene_modules <- c(rep.int('ECM Organization', times=86), rep.int('Immune response', times=42), rep.int('Metabolism',times=360), rep.int('Complement cascade', times=41), rep.int('Translation',times=1063))
length(gene_modules)
[1] 1592
dim(heatmap_input)
[1] 1592  103
anno_row <- data.frame(Gene_module = gene_modules)
rownames(anno_row) <- rownames(heatmap_input)
color_anno_all <- list(
	SWATH_class = c('1'='palegreen', '2'='lightpink', '3'='lightblue1'),
	DDA_class = c('1' = 'lightgoldenrod', '2'='lightpink', '3'='lightblue1'),
	mRNA_class = c('Differentiated'='red', 'Immunoreactive'='chartreuse', 'Mesenchymal'='deepskyblue3', 'Proliferative'='darkorchid4', 'Data not available'='white'),
	Original_Class = c('Differentiated'='indianred1', 'Immunoreactive'='darkolivegreen1', 'Mesenchymal'='skyblue2', 'Proliferative'='mediumorchid1', 'Stromal'='gold'),
	Gene_module=c('Complement cascade'='red', 'Immune response'='green', 'ECM Organization'='black', 'Translation'='lightblue', 'Metabolism'='grey'))

min(heatmap_input)
[1] -8.311899
max(heatmap_input)
[1] 3.702412

colors <- colorRampPalette(c('blue', 'white', 'red'))(n=10)
breaks <- c(seq(-8, -0.5, length.out=5), 0, seq(0.5, 8, length.out=5))
pheatmap(heatmap_input, color=colors, breaks=breaks, show_rownames=F, show_colnames=F,
filename='heatmap_SWATH_overlap_WGCNA_20171101.pdf', border_color=NA, cluster_cols=F, cluster_rows=F,
annotation_col=anno_cols_all, annotation_colors=color_anno_all, annotation_row=anno_row)


