## SWATH protein data preparation ##
## after mapDIA we use the Matrix fragments_for_protein_quantification for further steps

## Step 1: Prepare R session ##

sessionInfo()

setwd('')

## Step 2: Load in data ##

data <- read.delim('fragments_for_protein_quantification.txt', sep='\t', dec='.', as.is=T, header=T)
dim(data)

## Step 3: Filter data and select top3 fragments (we select overall top3 to have the same fragments in each sample)
# mapDIA filtered out fragments and set all their values to 0 and put them on the end

data_filtered <- data[1:90421,]
dim(data_filtered)

tail(data_filtered[,1])

data_filtered$fragment_sum <- apply(as.matrix(data_filtered[,4:106]), 1, sum, na.rm=T)

l <- data_filtered[order(data_filtered$fragment_sum, decreasing=T),]
d <- by(l, l$Peptide, head, n=3)
data_top3 <- Reduce(rbind, d)
dim(data_top3)

## Step 4: Aggregate fragments to peptides by summing the top3 fragments

length(unique(data_top3$Peptide))

data_peptide <- aggregate(data_top3[,4:106], by=list(data_top3$Protein, data_top3$Peptide), FUN=sum)

dim(data_peptide)

length(unique(data_peptide$Group.1))

## Step 5: Calculate peptide ratio to internal reference (peptide mean) this is comparable to the peptide ratios in iTRAQ

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

save(data_peptide, data_pepratio, data_logratio, file='CPTAC_SWATH_peptide_level.Rdata')
save(data_protein, file='CPTAC_SWATH_protein_level.Rdata')
save(data_top3, file='CPTAC_SWATH_transition_level_top3.Rdata')

## Step 7: Filtering with adaption to abundance level, high abundant proteins are allowed more missing values, while low abundant protein might have missing values due to technical reason and are allowed less missing values

#mean protein abundance will be calculated from peptide level, and protein matrix will be sorted accordingly
data_protein_abund <- aggregate(data_peptide[,3:105], by=list(data_peptide$Group.1), FUN=mean, na.rm=T)

mean_abund <- apply(data_protein_abund[,2:104], 1, mean, na.rm=T)
names(mean_abund)<- data_protein_abund[,1]
head(mean_abund)

mean_abund <- sort(mean_abund)
head(mean_abund)

data_protein <- as.matrix(data_protein)
data_protein_sorted <- data_protein[names(mean_abund),]

#count missing values per protein
SWATH_missing <- NULL
for (i in 1:nrow(data_protein_sorted)) {
tmp <- length(which(is.nan(as.matrix(data_protein_sorted[i,]))))
SWATH_missing <- c(SWATH_missing, tmp)
}
length(SWATH_missing)

names(SWATH_missing) <- rownames(data_protein_sorted)


#plot missing values vs. abundance (Figure S1C)
library(ggplot2)

pdf('missing_vs_abund_20180227.pdf', width=10, height=10)
ggplot(a)+
geom_smooth(aes(index, SWATH_missing), se=FALSE, size=2)+
geom_smooth(aes(index, 5*(log2(mean_abund))), colour='red', size=2)+
coord_cartesian(ylim=c(1,100))+
scale_y_continuous(sec.axis= sec_axis(~./5))+theme_bw()+
theme(axis.text.x = element_text(colour='black', size=35, angle=45, hjust=1), axis.text.y=element_text(colour='black', size=35), axis.title.x=element_text(size=0), axis.title.y=element_text(size=0))
dev.off()

#filter now based on percentiles
#lowest 10% abundant protein are allowed no missing values, the next 10% are allowed less than 10% missing values, the highest abundant 10% proteins are allowed 90% missing values
quantile(mean_abund, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
length(which(mean_abund<=3861.542)) 
#[1] 291
...
a <- names(which(SWATH_missing[1:291]==0))
...
SWATH_filtered <- data_protein_sorted[c(a,b,c,d,e,f,g,h,i,j),]
dim(SWATH_filtered)

#one protein could not be mapped to gene name and therefore removed
SWATH_filtered <- SWATH_filtered[grep('P62158', rownames(SWATH_filtered), invert=T),]
mapping_table <- read.delim('M2017110183C3DD8CE55183C76102DC5D3A26728BF22CD96.tab', sep='\t', dec='.', header=T, as.is=T)
mapping_table <- mapping_table[!duplicated(mapping_table$From),]
dim(mapping_table)
rownames(SWATH_filtered) <- mapping_table$To

#save output
save(data_protein_abund, mean_abund, data_protein_sorted, SWATH_missing, a, b, c, d, e, f, g, h, i, j, file='CPTAC_protein_filter_continuous.Rdata')
save(SWATH_filtered, file='CPTAC_protein_level_filtered.Rdata')

