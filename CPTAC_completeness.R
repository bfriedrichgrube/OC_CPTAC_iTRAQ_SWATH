setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_completeness/')
load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_DDA_protein_level_CDAP.Rdata')

load('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_Data/CPTAC_processed_data/CPTAC_SWATH_protein_level.Rdata')


SWATH_missing <- NULL
for (i in 1:nrow(data_protein)) {
tmp <- length(which(is.nan(as.matrix(data_protein[i,]))))
SWATH_missing <- c(SWATH_missing, tmp)}
length(SWATH_missing)
[1] 2914
length(which(SWATH_missing<93))
[1] 2807
length(which(SWATH_missing<83))
[1] 2616
length(which(SWATH_missing<73))
[1] 2401
length(which(SWATH_missing<62))
[1] 2169
length(which(SWATH_missing<52))
[1] 1990
length(which(SWATH_missing<42))
[1] 1786
length(which(SWATH_missing<31))
[1] 1523
length(which(SWATH_missing<21))
[1] 1272
length(which(SWATH_missing<11))
[1] 998
length(which(SWATH_missing==0))
[1] 347

SWATH_samples <- c("0%", "1-10%", "11-20%", "21-30%", "31-40%", "41-50%", "51-60%", "61-70%", "71-80%", "81-90%", "91-100%")
SWATH_ratio <- c(11.9, 34.2, 43.6, 52.3, 61.3, 68.3, 74.4, 82.4, 89.8, 96.3, 100)
SWATH_absolut <- c(347, 651, 274, 251, 263, 204, 179, 232, 215, 191, 107)
SWATH_complete<- data.frame(SWATH_samples, SWATH_ratio, SWATH_absolut)

library(ggplot2)

pdf('CPTAC_SWATH_missingvalue_20180205.pdf', width=10, height=10)
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


missing <- NULL
for(i in 1:nrow(DDA)){tmp <-length(which(is.na(DDA[i,])))
missing <- c(missing, tmp)
}

# substract 3 because first 3 lines are mean, std, median

dim(DDA)
[1] 8600  103
length(missing)
[1] 8600
8597
length(which(missing<93))
[1] 8541
8538
length(which(missing<83))
[1] 8283
8280
length(which(missing<73))
[1] 7977
7974
length(which(missing<62))
[1] 7651
7648
length(which(missing<52))
[1] 7360
7357
length(which(missing<42))
[1] 7038
7035
length(which(missing<31))
[1] 6696
6693
length(which(missing<21))
[1] 6286
6283
length(which(missing<11))
[1] 5773
5770
length(which(missing==0))
[1] 4366
4363

DDA_samples <- c("0%", "1-10%", "11-20%", "21-30%", "31-40%", "41-50%", "51-60%", "61-70%", "71-80%", "81-90%", "91-100%")
DDA_ratio <- c(50.7, 67.1, 73.1, 77.8, 81.8, 85.6, 89, 92.7, 96.3, 99.3, 100)
DDA_absolute <- c(4363, 1407, 513, 410, 342, 322, 291, 326, 306, 258, 59)
DDA_complete<- data.frame(DDA_samples, DDA_ratio, DDA_absolute)


pdf('CPTAC_DDA_missingvalue_20180205.pdf', width=10, height=10)
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
