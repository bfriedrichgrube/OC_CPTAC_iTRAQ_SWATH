#### Ovarian Cancer - John Hopkins #####
# New OpenSWATH run using new Library
# Library without window optimization, ciRTs
# OpenSWATH with iRTs, MST80, Jumboprophet
# no filtering needed anymore


## Step 1 - Prepare R session ##

#load required packages
library(SWATH2stats)
library(data.table)

# set working directory
setwd('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_SWATH2stats_preprocessing/CPTAC_SWATH2stats_2stepfiltering_final/103files_newlib_jumbo_20170207_preprocess')

#document information about versions for future reproducibility
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
[1] data.table_1.9.6  SWATH2stats_1.0.3

loaded via a namespace (and not attached):
[1] magrittr_1.5   plyr_1.8.3     tools_3.2.2    reshape2_1.4.1 Rcpp_0.12.7   
[6] stringi_1.0-1  grid_3.2.2     stringr_1.0.0  chron_2.3-47 

## Step 2 - Read in data ##

#read in data
data <- fread('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_SWATHsearch/103files_newlib_jumbo_20170206/feature_alignment.tsv', sep='\t', header=T)
Read 1077768 rows and 69 (of 69) columns from 1.233 GB file in 00:00:20

data <- as.data.frame(data)
dim(data)
[1] 1077768      69

#read in annotation file
annotation.file <- read.delim2('CPTAC_103HRD_annotation_20160929.txt', dec='.', sep='\t', header=T)
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
[1] 1077768      12

#remove unneeded proteins
data.reduced <- data.reduced[grep("iRT", data.reduced$ProteinName, invert=TRUE),]
dim(data.reduced)
[1] 1076428      12

data.reduced <- data.reduced[grep("DECOY", data.reduced$ProteinName, invert=TRUE),]
dim(data.reduced)
[1] 1067498      12

data.reduced$ProteinName <-gsub("^1\\/", "", data.reduced$ProteinName)
dim(data.reduced)
[1] 1067498      12

## Step 4 - Annotate data ##
#annotate the data
data.annotated <- sample_annotation(data.reduced, annotation.file)
dim(data.annotated)
[1] 1067498      15

## Step 5 - Transfer to data output ##

#splitting transitions needed for MSstats, mapDIA, aLFQ
data.transition <- disaggregate(data.annotated)
dim(data.transition)
[1] 6404988      10

# convert in mapDIA output
mapDIA.input<-convert4mapDIA(data.transition, RT=T)

dim(mapDIA.input)
[1] 150408    107
write.table(mapDIA.input, file="CPTAC_103HRD_jumbo_20170208_mapDIA_input.tsv", quote=F, row.names=F, sep="\t")


 Condition <- c('low', 'high', 'high', 'high', 'intermediate', 'low', 'low', 'intermediate', 'low', 'high',      'high', 'intermediate', 'intermediate', 'intermediate', 'low', 'low', 'intermediate', 'intermediate', 'intermediate', 'intermediate',        'intermediate', 'high', 'intermediate', 'high', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate',    'intermediate', 'low', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'high', 'intermediate', 'intermediate', 'intermediate',   'intermediate', 'intermediate', 'intermediate', 'intermediate', 'high', 'intermediate', 'intermediate', 'intermediate', 'high', 'intermediate',      'low', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'low', 'intermediate','high', 'intermediate', 'high', 'intermediate', 'low', 'high', 'high', 'intermediate', 'intermediate', 'intermediate','low', 'intermediate', 'intermediate', 'low', 'intermediate', 'low', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'high', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'low', 'intermediate', 'intermediate', 'intermediate', 'intermediate', 'low', 'low')


