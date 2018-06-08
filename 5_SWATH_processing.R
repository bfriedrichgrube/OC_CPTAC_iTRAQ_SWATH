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


