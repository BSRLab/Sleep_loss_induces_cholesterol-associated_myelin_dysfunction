#Affymetrix MG 430 2.0 packages
http://rss.acs.unt.edu/Rdoc/library/mouse4302/html/00Index.html
#A basic example is at
 http://bioinf.wehi.edu.au/marray/ibc2004/lab2/lab2.html#Estrogen


1) Install the required R packages using

source("http://bioconductor.org/biocLite.R")
biocLite("mouse4302cdf")
biocLite("affy")
biocLite("limma")

#there might be many dependencies among packages. In that case, R will compain that it is lacking some other packages. You can then install those.


library(mouse4302cdf)
library(affy)
library(limma)


2) Prepare a text file describing the files. This will look like as follows:

Name	FileName	Target
1BS		BS_08172011_(Mouse430_2).CEL	BS
1UBS	1UBS_08172011_(Mouse430_2).CEL	UBS
2BS		2BS_08172011_(Mouse430_2).CEL	BS
2UBS	2UBS_08172011_(Mouse430_2).CEL	UBS
3BW		3BW_08172011_(Mouse430_2).CEL	BW
3UBW	3UBW_08172011_(Mouse430_2).CEL	UBW
4BW		4BW_08172011_(Mouse430_2).CEL	BW
4UBW	4UBW_08172011_(Mouse430_2).CEL	UBW
5BSD	5BSD_08172011_(Mouse430_2).CEL	BSD
5UBSD	5UBSD_08172011_(Mouse430_2).CEL	UBSD
6BSD	6BSD_08172011_(Mouse430_2).CEL	BSD
6UBSD	6UBSD_08172011_(Mouse430_2).CEL	UBSD



I assumed BS, BW, BSD are the three experimental conditions and U* are their total versions. If this is not the case, 

3) Use ReadAffy function of the affy library to read data into R. This requires the cel files and the corresponding cdf file. cdf should be 
available after step 1.  cel files should be saved in a directory, with  the above phenotype file describing the samples.

setwd("")
pd <- read.AnnotatedDataFrame("file_desc.TXT", header = TRUE, row.names = 1, as.is = TRUE)
rawAffyData <- ReadAffy(filenames=pData(pd)$FileName, phenoData=pd)


4)  Quantile normalize within each group.

bs = rma(rawAffyData[, c(1, 3)])
ubs = rma(rawAffyData[, c(2, 4)])

bw = rma(rawAffyData[, c(5, 7)])
ubw = rma(rawAffyData[, c(6, 8)])

bsd = rma(rawAffyData[, c(9, 11)])
ubsd = rma(rawAffyData[, c(10, 12)])



5) Some diagnostic plots


# Before RMA normalization:
pdf("boxplot_before.pdf")
boxplot(rawAffyData, col="red", las = 2, names = pd@data[, 2])
dev.off()


pdf("hist_before.pdf")
hist(rawAffyData, type = "l", col = c(1, 2, 1, 2, 3, 4, 3, 4, 5, 6, 5, 6), lwd = 3)
dev.off()

# After RMA normalization:
pdf("boxplot_after_within_group.pdf")
boxplot(data.frame(cbind(exprs(bs), exprs(ubs), exprs(bw), exprs(ubw), exprs(bsd), exprs(ubsd))), col="blue", las = 2, names = 
c("BS", "BS", "UBS", "UBS", "BW", "BW", "UBW", "UBW", "BSD", "BSD", "UBSD", "UBSD"))
dev.off()


#Can also look at pairwise correlations of the replicates etc.


6) Combine all the samples
#pd1 = new("AnnotatedDataFrame", data= pd@data[ c(1, 3, 2, 4, 5, 7, 6, 8, 9, 11, 10, 12), ], varMetadata=pd@varMetadata, dimLabels=c("rowNames", "columnNames"))

expr <- new("ExpressionSet", exprs = cbind(exprs(bs), exprs(ubs), exprs(bw), exprs(ubw), exprs(bsd), exprs(ubsd)))


7) Median center all the samples
medcenter = function(unNorm){
# simple median centering
	medians = apply(unNorm, 2, median)
	Norm.tmp = sweep(unNorm, 2 ,medians)
	return(Norm.tmp)
}


exprs.temp = medcenter(exprs(expr))
exprs(expr) = exprs.temp


8) Simple clustering. Interestingly, good clustering of the treatment replicates, totals tend to cluster together.
c1 = hclust(dist(t(exprs(expr)),  method = "maximum"))
pdf("cluster_samples.pdf")
plot(c1)
dev.off()


9-a) We could filter low expressed genes...

9) Identify genes specific to each condition by comparing IP versus total within each treatment group.


#To identify genes expressed in BS compared to UBS
eset = exprs(expr)[, c(1:4)]
designBS = model.matrix(~ factor(c("BS", "BS", "UBS", "UBS"), levels = c("BS", "UBS")))
fit = lmFit(eset, designBS)
fit1 = eBayes(fit)
topTable(fit1, coef = 1, adjust = "BH")

## gives a list of probe sets that are DE between BS vs UBS.

10) Perform step 9 for BW vs. UBW and BSD vs. UBSD and take the union of probe sets that are DE in either of the three comparisons. 
Lets call the row
indices of the union as "index".

11) Perform differential expression analysis on the union probe set from 9.


designM = model.matrix(~ factor(c("BS", "BS", "BW", "BW", "BSD", "BSD"), levels = c("BS", "BW", "BSD")))
colnames(designM) = c("BS", "BW", "BSD")

eset = exprs(expr)[index, c(1:2, 5:6, 9:10)] #expression of all the probes sets that are specific to one or more of te treatment groups.

#perform all pairwise comparisons
fit <- lmFit(eset, designM)
contrast.matrix <- makeContrasts(BSD-BW, BSD-BS, BW-BS, levels = designM)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# alist of top genes for BSD versus BW comparison
topTable(fit2, coef = 1, adjust = "fdr")

results <- decideTests(fit2)
#A Venn diagram showing numbers of genes signi?cant in each comparison
vennDiagram(results)



