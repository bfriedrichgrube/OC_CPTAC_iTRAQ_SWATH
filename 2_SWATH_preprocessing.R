## OpenSWATH, pyprophet and TRIC as described in the publication ##

## Step 1 - Prepare R session ##

#load required packages
library(SWATH2stats)
library(data.table)

# set working directory
setwd('')

#document information about versions for future reproducibility
sessionInfo()

## Step 2 - Read in data ##

#read in data
data <- fread('feature_alignment.tsv', sep='\t', header=T)

data <- as.data.frame(data)
dim(data)

#read in annotation file (Condition column needs to be adjusted regarding potential group comparisons)
annotation.file <- read.delim2('annotation.txt', dec='.', sep='\t', header=T)
dim(annotation.file)
[1] 103   4
head(annotation.file)
                            Filename Condition BioReplicate Run
1 TCGA-09-1664-01A-01_S_W_JHU_141107       HRD            1   1
2 TCGA-13-1404-01A-01_S_W_JHU_141115   non_HRD            2   2
3 TCGA-13-1409-01A-01_S_W_JHU_141102   non_HRD            3   3
4 TCGA-13-1410-01A-01_S_W_JHU_141103   non_HRD            4   4
5 TCGA-13-1482-01A-01_S_W_JHU_141102       HRD            5   5
6 TCGA-13-1483-01A-01_S_W_JHU_141103       HRD            6   6
> 
## Step 3 - Reduce data ##

#remove unnecessary columns
data.reduced <- reduce_OpenSWATH_output(data)
dim(data.reduced)

#remove unneeded proteins
data.reduced <- data.reduced[grep("iRT", data.reduced$ProteinName, invert=TRUE),]
dim(data.reduced)

data.reduced <- data.reduced[grep("DECOY", data.reduced$ProteinName, invert=TRUE),]
dim(data.reduced)

data.reduced$ProteinName <-gsub("^1\\/", "", data.reduced$ProteinName)
dim(data.reduced)

## Step 4 - Annotate data ##
#annotate the data
data.annotated <- sample_annotation(data.reduced, annotation.file)
dim(data.annotated)

## Step 5 - Transfer to data output ##

#splitting transitions needed for MSstats, mapDIA, aLFQ
data.transition <- disaggregate(data.annotated)
dim(data.transition)

# convert in mapDIA output
mapDIA.input<-convert4mapDIA(data.transition, RT=T)
dim(mapDIA.input)

write.table(mapDIA.input, file="mapDIA_input.tsv", quote=F, row.names=F, sep="\t")


