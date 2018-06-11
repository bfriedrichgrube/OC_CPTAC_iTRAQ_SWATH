### Calculate missing values of unfiltered data of SWATH and iTRAQ DDA (Figure S1A-B)

## Step 1: Prepare R session

setwd('')
load('CPTAC_DDA_protein_level_CDAP.Rdata')
load('CPTAC_SWATH_protein_level.Rdata')

library(ggplot2)

sessionInfo()

# Step 2: Calculate number of missing values per percentile of samples and plot missing values numbers and percentage

SWATH_missing <- NULL
for (i in 1:nrow(data_protein)) {
tmp <- length(which(is.nan(as.matrix(data_protein[i,]))))
SWATH_missing <- c(SWATH_missing, tmp)}
length(SWATH_missing)
length(which(SWATH_missing<93))
...

SWATH_samples <- c("0%", "1-10%", "11-20%", "21-30%", "31-40%", "41-50%", "51-60%", "61-70%", "71-80%", "81-90%", "91-100%")
SWATH_ratio <- c(11.9, 34.2, 43.6, 52.3, 61.3, 68.3, 74.4, 82.4, 89.8, 96.3, 100)
SWATH_absolut <- c(347, 651, 274, 251, 263, 204, 179, 232, 215, 191, 107)
SWATH_complete<- data.frame(SWATH_samples, SWATH_ratio, SWATH_absolut)

#plot
pdf('CPTAC_SWATH_missingvalue.pdf', width=10, height=10)
ggplot(SWATH_complete)+
geom_bar(aes(x=SWATH_samples, y=SWATH_absolut), stat='identity', fill='lightskyblue', colour='lightskyblue')+
coord_cartesian(ylim=c(0, 1000))+
geom_point(aes(x=SWATH_samples, y=(SWATH_ratio*10), group=1),shape=15, size=3,stat='identity', colour='lightsalmon')+
geom_line(aes(x=SWATH_samples, y=(SWATH_ratio*10), group=1),stat='identity', colour='lightsalmon', size=2)+
scale_y_continuous(SWATH_complete, sec.axis = sec_axis(~./10))+scale_x_discrete(name='% missing values')+theme_bw()+
geom_text(aes(label=SWATH_absolut, x=SWATH_samples, y=SWATH_absolut-30), colour='dodgerblue4', size=9)+
geom_text(aes(label=paste(SWATH_ratio), x=SWATH_samples, y=(SWATH_ratio*10+30)), colour='salmon1', size=9)+
theme(axis.text.x = element_text(angle=45, hjust=1, colour='black', size=32), axis.text.y=element_text(colour='black', size=32), axis.title.x=element_text(size=35), axis.title.y=element_text(size=0))
dev.off()


# substract 3 because first 3 lines are mean, std, median
DDA <- DDA[4:8600, ]

missing <- NULL
for(i in 1:nrow(DDA)){tmp <-length(which(is.na(DDA[i,])))
missing <- c(missing, tmp)
}


length(missing)
length(which(missing<93))
...

DDA_samples <- c("0%", "1-10%", "11-20%", "21-30%", "31-40%", "41-50%", "51-60%", "61-70%", "71-80%", "81-90%", "91-100%")
DDA_ratio <- c(50.7, 67.1, 73.1, 77.8, 81.8, 85.6, 89, 92.7, 96.3, 99.3, 100)
DDA_absolute <- c(4363, 1407, 513, 410, 342, 322, 291, 326, 306, 258, 59)
DDA_complete<- data.frame(DDA_samples, DDA_ratio, DDA_absolute)


pdf('CPTAC_DDA_missingvalue.pdf', width=10, height=10)
ggplot(DDA_complete)+
geom_bar(aes(x=DDA_samples, y=DDA_absolute), stat='identity', fill='lightskyblue', colour='lightskyblue')+
coord_cartesian(ylim=c(0, 5000))+
geom_point(aes(x=DDA_samples, y=(DDA_ratio*50), group=1),shape=15, size=3,stat='identity', colour='lightsalmon')+
geom_line(aes(x=DDA_samples, y=(DDA_ratio*50), group=1),stat='identity', colour='lightsalmon', size=2)+
scale_y_continuous(DDA_complete, sec.axis = sec_axis(~./50))+scale_x_discrete(name='% missing values')+theme_bw()+
geom_text(aes(label=DDA_absolute, x=DDA_samples, y=DDA_absolute-120), size=9, colour='dodgerblue4')+
geom_text(aes(label=paste(DDA_ratio, sep=''), x=DDA_samples, y=(DDA_ratio*50-180)), colour='salmon1', size=9)+
theme(axis.text.x = element_text(angle=45, hjust=1, colour='black', size=32), axis.text.y=element_text(colour='black', size=32), axis.title.x=element_text(size=35), axis.title.y=element_text(size=0))
dev.off()
