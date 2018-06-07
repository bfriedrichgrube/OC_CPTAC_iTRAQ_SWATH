## SWATH protein data preparation ##


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

data <- read.delim('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_mapDIA_analysis/CPTAC_mapDIA302/CPTAC_103_wocomp_newlib_jumbo_min2max3_quantnorm_batchcorr/fragments_for_protein_quantification.txt', sep='\t', dec='.', as.is=T, header=T)
dim(data)
[1] 150408    106


## Step 3: Filter data and select top3 fragments


data_filtered <- data[1:90421,]
dim(data_filtered)
[1] 90421   106

tail(data_filtered[,1])

data_filtered$fragment_sum <- apply(as.matrix(data_filtered[,4:106]), 1, sum, na.rm=T)

l <- data_filtered[order(data_filtered$fragment_sum, decreasing=T),]
d <- by(l, l$Peptide, head, n=3)
data_top3 <- Reduce(rbind, d)
dim(data_top3)
[1] 47052   107

## Step 4: Aggregate fragments to peptides by summing the top3 fragments

length(unique(data_top3$Peptide))
[1] 15684

data_peptide <- aggregate(data_top3[,4:106], by=list(data_top3$Protein, data_top3$Peptide), FUN=sum)

dim(data_peptide)

length(unique(data_peptide$Group.1))


## Step 5: Calculate peptide ratio to internal reference (peptide mean)

peptide_mean <- apply(as.matrix(data_peptide[,3:105]), 1, mean, na.rm=T)
length(peptide_mean)

data_pepratio <- as.matrix(data_peptide[,3:105])/peptide_mean

data_logratio <- log2(data_pepratio)


## Step 6: Calculate protein ratios by averaging peptide ratios


data_protein <- aggregate(data_logratio, by=list(data_peptide$Group.1), FUN=mean, na.rm=T)

dim(data_protein)

rownames(data_protein) <- data_protein[,1]
data_protein[,1]<- NULL
dim(data_protein)
[1] 2914  103


ls()


save(data_peptide, data_pepratio, data_logratio, file='CPTAC_SWATH_peptide_level.Rdata')
save(data_protein, file='CPTAC_SWATH_protein_level.Rdata')
save(data_top3, file='CPTAC_SWATH_transition_level_top3.Rdata')
q()


