## DDA iTRAQ protein data preparation ##

## Step 1: Prepare R session ##

sessionInfo()

setwd('')

## Step 2: Load in data ##

# load in iTRAQ data downloaded from CPTAC data portal data acquired from John Hopkins University
DDA <- read.delim2('TCGA_Ovarian_JHU_Proteome.itraq.tsv', sep='\t', dec='.', header=T, as.is=T)
dim(DDA)
[1] 8600  271

## Step 3: Filter just for ratios calculated from unshared peptides

#extract every second column which contains measurements from proteotypic peptides only
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

## Step 5: Filter for proteins without missing values

missing <- NULL
for(i in 1:nrow(DDA)){tmp <-length(which(is.na(DDA[i,])))
missing <- c(missing, tmp)
}

DDA_missing0 <- DDA[(missing==0),] 
dim(DDA_missing0)

## Step 6: Save output

save(DDA, DDA_missing0, file='CPTAC_DDA_protein_level_CDAP.Rdata')
