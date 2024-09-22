## ----setup, echo=FALSE------------------------------------
library(knitr)
opts_chunk$set(cache=FALSE, message=FALSE, echo=TRUE, results="markup")
options(width=60)

## ----data-------------------------------------------------
library(EDASeq)
library(yeastRNASeq)
library(leeBamViews)

input <- function(inputfile) {
	mypref <<- inputfile
}

run <- function() {}

output <- function(outputfile) {
## ----import-raw-------------------------------------------
#files <- list.files(file.path(system.file(package = "yeastRNASeq"),
#                              "reads"), pattern = "fastq", full.names = TRUE)
files <- list.files(file.path(mypref), pattern="fastq", full.names=TRUE)
names(files) <- gsub("\\.fastq.*", "", basename(files))
met <- DataFrame(conditions=c(rep("mut",2), rep("wt",2)),
                 row.names=names(files))
fastq <- FastqFileList(files)
elementMetadata(fastq) <- met
fastq

## ----import-aligned---------------------------------------
#files <- list.files(file.path(system.file(package = "leeBamViews"), "bam"),
#                    pattern = "bam$", full.names = TRUE)
files <- list.files(file.path(mypref), pattern="bam$", full.names=TRUE)

names(files) <- gsub("\\.bam", "", basename(files))

gt <- gsub(".*/", "", files)
gt <- gsub("_.*", "", gt)
lane <- gsub(".*(.)$", "\\1", gt)
geno <- gsub(".$", "", gt)

pd <- DataFrame(geno=geno, lane=lane,
                row.names=paste(geno,lane,sep="."))

bfs <- BamFileList(files)
elementMetadata(bfs) <- pd
bfs

## ----plot-total-------------------------------------------
colors <- c(rep(rgb(1,0,0,alpha=0.7),2),
            rep(rgb(0,0,1,alpha=0.7),2),
            rep(rgb(0,1,0,alpha=0.7),2),
            rep(rgb(0,1,1,alpha=0.7),2))
barplot(bfs,las=2,col=colors)

## ----plot-total-2-----------------------------------------
plotQuality(bfs,col=colors,lty=1)
legend("topright",unique(elementMetadata(bfs)[,1]), fill=unique(colors))

## ----plot-qual--------------------------------------------
plotQuality(bfs[[1]],cex.axis=.8)
barplot(bfs[[1]],las=2)

## ----plot-nt----------------------------------------------
plotNtFrequency(bfs[[1]])

## ----load-gene-level--------------------------------------
data(geneLevelData)
head(geneLevelData)

## ----load-lgc---------------------------------------------
data(yeastGC)
head(yeastGC)
data(yeastLength)
head(yeastLength)

## ----filter-----------------------------------------------
filter <- apply(geneLevelData,1,function(x) mean(x)>10)
table(filter)
common <- intersect(names(yeastGC),
                    rownames(geneLevelData[filter,]))
length(common)

## ----create-object----------------------------------------
feature <- data.frame(gc=yeastGC,length=yeastLength)
data <- newSeqExpressionSet(counts=as.matrix(geneLevelData[common,]),
                            featureData=feature[common,],
                            phenoData=data.frame(
                              conditions=factor(c(rep("mut",2),rep("wt",2))),
                              row.names=colnames(geneLevelData)))
data

## ----show-data--------------------------------------------
head(counts(data))
pData(data)
head(fData(data))

## ----show-offset------------------------------------------
head(offst(data))

## ----boxplot-genelevel------------------------------------
boxplot(data,col=colors[1:4])

## ----md-plot----------------------------------------------
MDPlot(data,c(1,3))

## ----plot-mean-var----------------------------------------
meanVarPlot(data[,1:2], log=TRUE, ylim=c(0,16))
meanVarPlot(data, log=TRUE, ylim=c(0,16))

## ----plot-gc----------------------------------------------
biasPlot(data, "gc", log=TRUE, ylim=c(1,5))

## ----boxplot-gc-------------------------------------------
lfc <- log(counts(data)[,3]+0.1) - log(counts(data)[,1]+0.1)
biasBoxplot(lfc, fData(data)$gc)

## ----normalization----------------------------------------
dataWithin <- withinLaneNormalization(data,"gc", which="full")
dataNorm <- betweenLaneNormalization(dataWithin, which="full")

## ----plot-gc-norm-----------------------------------------
biasPlot(dataNorm, "gc", log=TRUE, ylim=c(1,5))
boxplot(dataNorm, col=colors)

## ----norm-offset------------------------------------------
dataOffset <- withinLaneNormalization(data,"gc",
                                      which="full",offset=TRUE)
dataOffset <- betweenLaneNormalization(dataOffset,
                                       which="full",offset=TRUE)

## ----edger------------------------------------------------
library(edgeR)
design <- model.matrix(~conditions, data=pData(dataOffset))

y <- DGEList(counts=counts(dataOffset),
             group=pData(dataOffset)$conditions)
y$offset <- -offst(dataOffset)
y <- estimateDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

## ----deseq------------------------------------------------
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(dataOffset),
                              colData = pData(dataOffset),
                              design = ~ conditions)

normFactors <- exp(-1 * offst(dataOffset))
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors

dds <- DESeq(dds)
res <- results(dds)
res

## ----unrounded--------------------------------------------
dataNorm <- betweenLaneNormalization(data, round=FALSE, offset=TRUE)

norm1 <- normCounts(dataNorm)
norm2 <- exp(log(counts(dataNorm) + 0.1 ) + offst(dataNorm)) - 0.1

#head(norm1 - norm2)
## ----rounded----------------------------------------------
write.csv(round(normCounts(dataNorm)) - round(counts(dataNorm) * exp(offst(dataNorm))), paste(outputfile, "rounded", "csv", sep="."))

write.csv(norm1-norm2, paste(outputfile, "diff", "csv", sep="."))

}


## ----getLengthAndGC, eval=FALSE---------------------------
#  getGeneLengthAndGCContent(id=c("ENSG00000012048", "ENSG00000139618"), org="hsa")

## ----getLengthAndGC-full, eval=FALSE----------------------
#  fData(data) <- getGeneLengthAndGCContent(featureNames(data),
#                                                org="sacCer3", mode="org.db")

## ----sessionInfo------------------------------------------
#sessionInfo()


