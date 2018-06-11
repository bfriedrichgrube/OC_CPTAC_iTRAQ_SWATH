
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_raw_data/CPTAC_peptide_ratios_103.Rdata')
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_SWATH_peptide_level.Rdata')

SWATH_log_pepratio <- cbind(data_peptide[,1:2], data_logratio)

gene.ids <- intersect(SWATH_log_pepratio[,1], data_log_pepratio[,1])




SWATH_overlap <- SWATH_log_pepratio[SWATH_log_pepratio[,1] %in% gene.ids,]
dim(SWATH_overlap)
[1] 14542   105
head(gene.ids)
[1] "Q53GG5" "O94826" "Q8N4T8" "P35221" "P13807" "Q8NDA8"
length(gene.ids)
[1] 2694
length(unique(SWATH_overlap[,1]))
[1] 2694

DDA_overlap <- data_log_pepratio[data_log_pepratio[,1] %in% gene.ids, ]


SWATH_range <- as.matrix(SWATH_overlap[,3:105])
SWATH_range <- (SWATH_range-min(SWATH_range, na.rm=T))/(max(SWATH_range, na.rm=T)-min(SWATH_range, na.rm=T))
SWATH_range <- SWATH_range[,c(1,66,74,80,23,27,92,97,51,2,6:7,60,8,61:65,67:69,9:10,70:73,
11,75:76,12:14,77,15,78:79,16:17,81:82,18:20,83,21,22,84:90,24,91,25:26,28:41,93,42,94,43,95:96,98:99,44:50,100:101,52,102,53:58,103,59,3:5)]

DDA_range <- as.matrix(DDA_overlap[,3:105])
DDA_range[which(is.infinite(DDA_range))]<-NA
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



