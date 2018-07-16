## Differential expression analysis of 103 HGS-OC samples with diffferent HRD status (Figure 6A)

## After mapDIA ##

## Step 1: Prepare R session ## 

library(calibrate)

setwd('')

## Step 2: Load in data ##

output <- read.table('analysis_output.txt', sep='\t', dec='.', header=T, as.is=T)

## Step 3: Plot differentially expressed proteins in Volcano plot ##

pdf('Volcano_HRD.pdf', width=10, height=10)
with(output, plot(log2FC, -log10(FDR), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-1, 1)))
with(subset(output, FDR<0.05 & log2FC>log2(1.3)), points(log2FC, -log10(FDR), pch=20, cex=2.2, col='red'))
with(subset(output, FDR<0.05 & log2FC<(-log2(1.3))), points(log2FC, -log10(FDR), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()

