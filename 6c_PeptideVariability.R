## Peptide vaiablility comparison

## Step 1: Peptide values need to be calculated manually from PSM files downloaded from CPTAC data portal

#first PSM files need to be mapped to biomart genes
library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
biomart_data <- getBM(attributes=c('refseq_peptide','uniprotswissprot','hgnc_symbol'), mart=mart)
#mapping was performed parallel in batches of 4 to 5 PSM files
result <- list()
files <- list.files()[74:78]

for (m in 1:5){
	print(Sys.time())
	data <- read.delim2(files[m], sep='\t', dec='.', header=T, as.is=T)
	print(dim(data))
	data <- data[which(data$AmbiguousMatch==0),]
	tmp <- strsplit(data[1,1], '_')

	data <- data[,c(12,14, 23:26)]

	colnames(data)[3:6] <- tmp[[1]][2:5]

	tmp2 <- strsplit(data$Protein, split="\\|")
	uniprot<- vector()
	tmp3 <- list()
	tmp4<- list()
	tmp5<- list()
	tmp6 <- NULL
	print(dim(data))
	print(length(tmp2))
	for (i in 1:length(tmp2)){
		tmp3[[i]]<- grep('NP_', tmp2[[i]], value=T)
		
		if (length(tmp3[[i]])==0){
			uniprot[i] <- NA} else {
			
			tmp4[[i]]<- unlist(strsplit(tmp3[[i]], '\\.'))
			tmp5[[i]]<- grep('NP_', tmp4[[i]], value=T)
			tmp_annotation <- NULL
			for (j in 1:length(tmp5[[i]])){
				tmp6 <- biomart_data[grep(tmp5[[i]][j], biomart_data$refseq_peptide),]
				tmp_annotation <- c(tmp_annotation, tmp6$uniprotswissprot)	
			}
			if (length(unique(tmp_annotation))==1){
				uniprot[i] <- unique(tmp_annotation)
				} else { 
				uniprot[i] <- NA 
			}
		}
	}
	data_x <- cbind(data, uniprot)
	print(dim(data_x))
	result[[m]]<- data_x
}
save(data_x, ...)

#data_x is replaced by data_1, data_2, ...
setwd('')
files <- list.files()
for (i in 1: length(files)){
load(files[i])
}

## make list of all data.frames (44)
list_data <- list(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17, data_18, data_19, data_20, data_21, data_22, data_23, data_24, data_25, data_26, data_27, data_28, data_29, data_30, data_31, data_32, data_33, data_34, data_35, data_36, data_37, data_38, data_39, data_40, data_41, data_42, data_43, data_44)

## remove all tag flags from peptide sequence and take out just the peak area from columns
## then create new data.frame consisting of uniprot, peptide and 4 iTRAQ measurement columns
## then aggregate peptides by summing
## then calculate peptide ratios for each iTRAQ batch
## store again in List

data_ratios <- list()
for (i in 1:length(list_data)){
	tmp_data <- list_data[[i]]
	tmp_filtered <- tmp_data[which(!is.na(as.character(tmp_data$uniprot))),c(1,3:7)]
	tmp_names <- c('Protein', 'Peptide', colnames(tmp_filtered)[3:5])
	tmp_pep <- gsub('[[:digit:]]', '', tmp_filtered[,1])
	tmp_pep <- gsub('[[:punct:]]', '', tmp_pep)
	tmp1 <- as.numeric(gsub('(.*)\\/(.*)','\\1', tmp_filtered[,2]))
	tmp2 <- as.numeric(gsub('(.*)\\/(.*)','\\1', tmp_filtered[,3]))
	tmp3 <- as.numeric(gsub('(.*)\\/(.*)','\\1', tmp_filtered[,4]))
	tmp4 <- as.numeric(gsub('(.*)\\/(.*)','\\1', tmp_filtered[,5]))
	tmp_new <- data.frame(as.character(tmp_filtered[,6]),as.character(tmp_pep),as.numeric(tmp1),as.numeric(tmp2),as.numeric(tmp3),as.numeric(tmp4), stringsAsFactors=F)
	tmp_data_pep <- aggregate(tmp_new[,3:6], by=list(tmp_new[,1], tmp_new[,2]), sum, na.rm=T)
	tmp_ratios <- data.frame(tmp_data_pep[,1:2], tmp_data_pep[,4]/tmp_data_pep[,3], tmp_data_pep[,5]/tmp_data_pep[,3], tmp_data_pep[,6]/tmp_data_pep[,3])
	colnames(tmp_ratios) <- tmp_names	
	data_ratios[[i]] <- tmp_ratios
}

length(data_ratios)
## merge all data.frame together
data_merged <- data_ratios[[1]]

for (j in 2:length(data_ratios)){
	data_merged <- merge(data_merged, data_ratios[[j]], all=T)
}
dim(data_merged)
## remove all rows where no Uniprot ID was mapped
data_merged_filtered <- data_merged[which(data_merged[,1]!=''),]
dim(data_merged_filtered)
## aggregate again, to combine separate rows from merging
data_merged_filter_aggr <- aggregate(data_merged_filtered[,3:125], by=list(data_merged_filtered[,1], data_merged_filtered[,2]), mean, na.rm=T)
dim(data_merged_filter_aggr)
data_peptide <- data_merged_filter_aggr[,3:125]
data_peptide <- data_peptide[,order(colnames(data_peptide))]
data_peptide <- cbind(data_merged_filter_aggr[,1:2], data_peptide[,c(1, 3:12, 14:16, 18:33, 37:42, 44:51, 53:58, 60:69, 71:79, 81, 83:85, 87:96, 98:103,  106:113, 115:118, 121:122)])
## log2 ratio
data_log_pepratio <- log2(as.matrix(data_peptide[,3:105]))
length(which(is.infinite(data_log_pepratio)))
length(which(is.na(data_log_pepratio)))
dim(data_log_pepratio)

colnames(data_log_pepratio)[1:2]<- c('Protein', 'Peptide')



## Step 2: Load SWATH data additionally and select common peptides

load('CPTAC_SWATH_peptide_level.Rdata')
SWATH_log_pepratio <- cbind(data_peptide[,1:2], data_logratio)
gene.ids <- intersect(SWATH_log_pepratio[,1], data_log_pepratio[,1])

SWATH_overlap <- SWATH_log_pepratio[SWATH_log_pepratio[,1] %in% gene.ids,]
dim(SWATH_overlap)
head(gene.ids)
length(gene.ids)
length(unique(SWATH_overlap[,1]))

DDA_overlap <- data_log_pepratio[data_log_pepratio[,1] %in% gene.ids, ]

## Step 3: SWATH and iTRAQ DDA data are centered to fit the same range between 0 and 1 to be compared directly protein-wise comparsion (not published)

SWATH_range <- as.matrix(SWATH_overlap[,3:105])
SWATH_range <- (SWATH_range-min(SWATH_range, na.rm=T))/(max(SWATH_range, na.rm=T)-min(SWATH_range, na.rm=T))
SWATH_range <- SWATH_range[,c(1,66,74,80,23,27,92,97,51,2,6:7,60,8,61:65,67:69,9:10,70:73,
11,75:76,12:14,77,15,78:79,16:17,81:82,18:20,83,21,22,84:90,24,91,25:26,28:41,93,42,94,43,95:96,98:99,44:50,100:101,52,102,53:58,103,59,3:5)]

DDA_range <- as.matrix(DDA_overlap[,3:105])
DDA_range[which(is.infinite(DDA_range))]<- NA
DDA_range <- (DDA_range-min(DDA_range, na.rm=T))/(max(DDA_range, na.rm=T)-min(DDA_range, na.rm=T))
colnames(SWATH_range) <- colnames(DDA_range)

DDA_range <- cbind(DDA_overlap[,1:2], DDA_range)
SWATH_range <- cbind(SWATH_overlap[,1:2], SWATH_range)
SWATH_order <- SWATH_range[order(as.character(SWATH_range$Group.1)),]
DDA_order <- DDA_range[order(as.character(DDA_range$Protein)),]
SWATH_split <- split(SWATH_order, SWATH_order[,1])
DDA_split <- split(DDA_order, DDA_order[,1])


library(reshape2)
plot_final <- list()
for (i in 1:2694){
SWATH_input <- melt(SWATH_split[[i]][,2:105], id=c('Group.2'))
SWATH_input$value <- SWATH_input$value/(-1)
DDA_input <- melt(DDA_split[[i]][,2:105], id=c('Peptide'))
colnames(SWATH_input)[1]<- 'Peptide'
input <- rbind(SWATH_input, DDA_input)
tmp <- ggplot(input)+geom_point(aes(x=variable, y=value, colour=Peptide))+theme(legend.position='none', axis.text.x=element_blank(), plot.title=element_text(size=24), axis.title=element_text(size=20))+coord_cartesian(ylim=c(-1,1))+xlab('Samples')+ylab('SWATH                                             iTRAQ-DDA')+labs(title=paste(unique(SWATH_order[,1])[i]))
plot_final[[i]]<- tmp
}

## Step 4: SWATH and iTRAQ DDA data is unlogged to compared peptide variability in an uncompressed format
# a CV-like score is calculated by dividing the standard deviation of each protein by its mean

SWATH_unlog <- 2**(as.matrix(SWATH_overlap[,3:105]))
SWATH_mean <- aggregate(SWATH_unlog, by=list(SWATH_overlap$Group.1), FUN=mean,  na.rm=T)
SWATH_sd <- aggregate(SWATH_unlog, by=list(SWATH_overlap$Group.1), FUN=sd,  na.rm=T)
dim(SWATH_mean)
SWATH_CV <- as.matrix(SWATH_sd[,2:104])/as.matrix(SWATH_mean[,2:104])
min(SWATH_CV, na.rm=T)
max(SWATH_CV, na.rm=T)
vioplot(as.numeric(na.omit(SWATH_CV))

DDA_unlog <- 2**(as.matrix(DDA_overlap[,3:105]))
dim(DDA_unlog)
DDA_mean <- aggregate(DDA_unlog, by=list(DDA_overlap$Protein), FUN=mean,  na.rm=T)
dim(DDA_mean)
DDA_sd <- aggregate(DDA_unlog, by=list(DDA_overlap$Protein), FUN=sd,  na.rm=T)
DDA_CV <- as.matrix(DDA_sd[,2:104])/as.matrix(DDA_mean[,2:104])

#plot Figure 2C
pdf('peptide_per_protein_CVs_20180314.pdf', height=10, width=10)
par(cex.axis=2.5, las=1)
vioplot(as.numeric(na.omit(DDA_CV)), as.numeric(na.omit(SWATH_CV)), col='grey', names=c('', ''))
dev.off()
