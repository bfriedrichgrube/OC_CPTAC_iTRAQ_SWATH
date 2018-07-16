## Differential expression analysis of 103 HGS-OC samples with diffferent HRD status ##

## After mapDIA ##

## Step 1: Prepare R session ## 

library(calibrate)


setwd('Documents/Biomarker_Cancer/Ovarian_JohnHopkins/CPTAC_mapDIA_analysis/CPTAC_mapDIA302/CPTAC_103HRD_newlib_jumbo_min2max3_quantnorm_batchcorr/')

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
[1] calibrate_1.7.2 MASS_7.3-43  

## Step 2: Load in data ##

output <- read.table('analysis_output.txt', sep='\t', dec='.', header=T, as.is=T)

## Step 3: Plot differentially expressed proteins in Volcano plot ##

pdf('Volcano.pdf', width=10, height=10)
with(output, plot(log2FC, -log10(FDR), pch=20, main='Volcano Plot HRD vs. non-HRD', xlim=c(-1.2, 1.2), ylim=c(0,20)))
with(subset(output, FDR<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(FDR), pch=20, col='red'))
with(subset(output, FDR<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(FDR), pch=20, col='blue'))
with(subset(output, FDR<0.05 & abs(log2FC)>log2(1.3)),textxy(log2FC, -log10(FDR), labs=Protein, cex=0.9, offset=0.5))
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()

pdf('Volcano_HRD_20180315.pdf', width=10, height=10)
with(output, plot(log2FC, -log10(FDR), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-1, 1)))
with(subset(output, FDR<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(FDR), pch=20, cex=2.2, col='red'))
with(subset(output, FDR<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(FDR), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()





output[which(output$FDR<0.05 & output$log2FC>log2(1.3)),]
     Protein nPeptide nFragment Label      Label2   log2FC    score SignedScore
876   P31323        2        11   0/1 HRD/non_HRD 0.515806 0.999655    0.999655
1053  P48449        2        12   0/1 HRD/non_HRD 0.516672 0.897385    0.897385
1113  P50583        2        12   0/1 HRD/non_HRD 0.626117 0.950158    0.950158
1879  Q8IU81        3        17   0/1 HRD/non_HRD 0.414399 0.741266    0.741266
2384  Q9HCG8        2        11   0/1 HRD/non_HRD 0.487901 0.999640    0.999640
2702  Q9Y653        3        16   0/1 HRD/non_HRD 0.501883 0.914642    0.914642
             FDR log_oddsDE
876  1.45267e-05    7.97089
1053 1.70606e-02    2.16850
1113 8.88357e-03    2.94777
1879 4.94975e-02    1.05256
2384 2.83070e-05    7.92919
2702 1.49751e-02    2.37167

> output[which(output$FDR<0.05 & output$log2FC<(-log2(1.3))),]
     Protein nPeptide nFragment Label      Label2    log2FC    score
7     A2VCL2        2        12   0/1 HRD/non_HRD -0.515245 0.999880
367   P02511        2        12   0/1 HRD/non_HRD -0.508338 0.999982
378   P02743        2        12   0/1 HRD/non_HRD -0.390664 0.997624
397   P03915        2        11   0/1 HRD/non_HRD -0.579959 0.976663
399   P04004        4        24   0/1 HRD/non_HRD -0.515707 1.000000
630   P14174        2        12   0/1 HRD/non_HRD -0.568619 0.999889
683   P17661        6        34   0/1 HRD/non_HRD -0.831290 1.000000
756   P22676        2        12   0/1 HRD/non_HRD -0.442724 0.956766
760   P23141        4        23   0/1 HRD/non_HRD -0.398062 0.891371
1042  P47712        2        11   0/1 HRD/non_HRD -0.664734 0.965656
1141  P51854        3        18   0/1 HRD/non_HRD -0.787940 1.000000
1622  Q15124        2        11   0/1 HRD/non_HRD -0.538398 0.995556
2005  Q8WXE9        2        10   0/1 HRD/non_HRD -0.426719 0.999988
2297  Q9GZM7        3        18   0/1 HRD/non_HRD -0.578226 0.999187
     SignedScore         FDR log_oddsDE
7      -0.999880 9.72280e-06    9.02613
367    -0.999982 2.64006e-06   10.92240
378    -0.997624 3.56036e-04    6.03977
397    -0.976663 2.72929e-03    3.73409
399    -1.000000 7.40149e-17   35.08110
630    -0.999889 4.90624e-06    9.10704
683    -1.000000 0.00000e+00   37.87170
756    -0.956766 7.58659e-03    3.09693
760    -0.891371 2.07880e-02    2.10483
1042   -0.965656 5.42341e-03    3.33637
1141   -1.000000 3.04921e-09   17.36950
1622   -0.995556 5.76182e-04    5.41171
2005   -0.999988 2.10885e-06   11.31390
2297   -0.999187 8.06516e-05    7.11358



P31323 low association with cancer
P48449 low association with cancer
P50583 no association with cancer
Q8IU81 no association with cancer
Q9HCG8 low association with cancer
Q9Y653 low association with cancer, know to activate VEGFA prodoction, angiogenesis


A2VCL2 no association with cancer
P02511 no association with cancer
P02743 low association with cancer, mutated in breast cancer
P03915 low association with cancer
P04004 discussed association, medium
P14174 medium association with cancer
P17661 low association with cancer
P22676 low association with cancer
P23141 low association with cancer
P47712 low association with cancer
P51854 low association with cancer
Q15124 low association with cancer
Q8WXE9 low association with cancer
Q9GZM7 low association with cancer


