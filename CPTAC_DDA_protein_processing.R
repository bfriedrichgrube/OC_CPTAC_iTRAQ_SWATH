## DDA protein data preparation ##

## Step 1: Prepare R session ##

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

setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/')

## Step 2: Load in data ##

DDA <- read.delim2('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_raw_data/TCGA_Ovarian_JHU_Proteome.itraq.tsv', sep='\t', dec='.', header=T, as.is=T)
dim(DDA)
[1] 8600  271

## Step 3: Filter just for ratios calculated from unshared peptides

a <- seq(1,271, 2)
DDA <- DDA[,a]
dim(DDA)
[1] 8600  136

## Step 4: Filter out columns that were not in common with SWATH

rownames(DDA)<- DDA[,1]
DDA[,1]<- NULL
DDA <- DDA[,1:132]
DDA <-DDA[,order(colnames(DDA))] 
DDA <- DDA[,11:132]
dim(DDA)
[1] 8600  122
DDA <- DDA[,c(1,3:12, 14:16, 18:33, 37:42, 44:51, 53:58, 60:69, 71:79, 81, 83:85, 87:96, 98:103, 106:113, 115:118, 121:122)]
dim(DDA)
[1] 8600  103
