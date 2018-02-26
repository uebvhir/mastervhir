#---------------------------------------------------------------------------------------------
###THIS IS AN EXAMPLE CODE FOR THE ANALYSIS OF AFFYMETRIX GENE MICROARRAYS
#---------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
###FOLDER DESTINATION DEFINITIONS
#---------------------------------------------------------------------------------------------
# In Rstudio do: 
# Start doing "Session --> Set working directory --> To source file location" 
# 
workingDir <-getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir, "results")
celfilesDir <- file.path(workingDir, "celfiles")

#---------------------------------------------------------------------------------------------
###INSTALLATION OF PACKAGES NEEDED
#---------------------------------------------------------------------------------------------

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(pkg)
}else{
  require(pkg, character.only=T)
  }
}
installifnot("org.Mm.eg.db")
installifnot("pd.mogene.1.1.st.v1")
installifnot("mogene11sttranscriptcluster.db")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")

#-----------------------------------------------------------------
# RANDOM SELECTION OF SAMPLES
#-----------------------------------------------------------------
targetsAll <-read.csv(file=file.path(dataDir,"fileList.csv"), header = TRUE, sep=",", row.names=1) 
targetsAll

ID <- 1234567
muscle <- ifelse (ID %% 2 == 0, "quadriceps", "soleus") 
set.seed <- ID
mostres <-  c(sample (1:6,4), sample (7:12,4))
targets <- targetsAll[mostres,]

#------------------------------------------------------------------
### LOAD DATA: READ CEL FILES. 
#------------------------------------------------------------------

#CELFILES
require(oligo)
CELfiles <- file.path(celfilesDir, rownames(targets))
CELfiles
rawData <- read.celfiles(CELfiles)

#DEFINE SOME VARIABLES FOR PLOTS
sampleNames <- as.character(targets$samplename)
sampleColor <- as.character(as.integer(targets$group)+1)

#---------------------------------------------------------------------------------------------
###QUALITY CONTROL OF ARRAYS: RAW DATA
#---------------------------------------------------------------------------------------------

#BOXPLOT
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)

#HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)

#PRINCIPAL COMPONENT ANALYSIS
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

#SAVE TO A FILE
pdf(file.path(resultsDir, "QCPlots_Raw.pdf"))
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples of RawData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()


#---------------------------------------------------------------------------------------------
###DATA NORMALIZATION
#---------------------------------------------------------------------------------------------
eset<-rma(rawData)

write.exprs(eset, file.path(resultsDir, "NormData.txt"))


#---------------------------------------------------------------------------------------------
###QUALITY CONTROL OF ARRAYS: NORMALIZED DATA
#---------------------------------------------------------------------------------------------

#BOXPLOT
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)

#HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)

#PRINCIPAL COMPONENT ANALYSIS

plotPCA(exprs(eset), labels=sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

#SAVE TO A FILE
pdf(file.path(resultsDir, "QCPlots_Norm.pdf"))
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(eset), labels=sampleNames, dataDesc="selected samples", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()

#ARRAY QUALITY METRICS
# arrayQualityMetrics(eset,  reporttitle="QualityControl", force=TRUE)


#---------------------------------------------------------------------------------------------
###FILTER OUT THE DATA 
#---------------------------------------------------------------------------------------------

annotation(eset) <- "org.Mm.eg.db"
require(org.Mm.eg.db)
require(genefilter)
eset_filtered <- nsFilter(eset, var.func=IQR,
         var.cutoff=0.75, var.filter=TRUE,
         filterByQuantile=TRUE)
#NUMBER OF GENES OUT
print(eset_filtered$filter.log$numLowVar)

#NUMBER OF GENES IN
print(eset_filtered$eset)


#---------------------------------------------------------------------------------------------
###DIFERENTIAL EXPRESSED GENES SELECTION. LINEAR MODELS. COMPARITIONS
#---------------------------------------------------------------------------------------------

#CONTRAST MATRIX.lINEAR MODEL
require(limma)
treat <- targets$group
lev <- factor(treat, levels = unique(treat))
design <-model.matrix(~0+lev)
colnames(design) <- levels(lev)
rownames(design) <- sampleNames
print(design)

#COMPARISON
cont.matrix1 <- makeContrasts( 
        KO.vs.CTL = cKO-CTL,
        levels = design)
comparison1 <- "Effect of Knocking Out gene"

#MODEL FIT
fit1 <- lmFit(eset_filtered$eset, design)
fit.main1 <- contrasts.fit(fit1, cont.matrix1)
fit.main1 <- eBayes(fit.main1)

#---------------------------------------------------------------------------------------------
### DIFERENTIAL EXPRESSED GENES LISTS.TOPTABLES
#---------------------------------------------------------------------------------------------

#FILTER BY FALSE DISCOVERY RATE AND FOLD CHANGE =2.
topTab <-  topTable (fit.main1, number=nrow(fit.main1), coef=1, adjust="fdr", lfc=abs(2))
dim(topTab)

#EXPORTED TO CSV AND HTML FILE
write.csv2(topTab, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison1, ".csv", sep = "")))

print(xtable(topTab,align="lllllll"),type="html",html.table.attributes="",
      file=paste("Selected.Genes.in.comparison.",comparison1,".html", sep=""))

#---------------------------------------------------------------------------------------------
###VOLCANO PLOTS
#---------------------------------------------------------------------------------------------

volcanoplot(fit.main1, highlight=10, names=fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep="\n"))
abline(v = c(-2, 2))


pdf(file.path(resultsDir,"Volcanos.pdf"))
volcanoplot(fit.main1, highlight = 10, names = fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep = "\n"))
abline(v = c(-2, 2))
dev.off()


#---------------------------------------------------------------------------------------------
###HEATMAP PLOTS
#---------------------------------------------------------------------------------------------

#PREPARE THE DATA
my_frame <- data.frame(exprs(eset))
head(my_frame)
HMdata <- merge(my_frame, topTab, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
head(HMdata)
HMdata2 <- data.matrix(HMdata, rownames.force=TRUE)
head(HMdata2)
write.csv2(HMdata2, file = file.path(resultsDir,"Data2HM.csv"))

#HEATMAP PLOT
my_palette <- colorRampPalette(c("blue", "red"))(n = 299)
require(gplots)
heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap Induced.vs.WT FC>=3",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=c(rep("red",4),rep("blue",4)),
          tracecol=NULL,
          srtCol=30)

#EXPORT TO PDF FILE
pdf(file.path(resultsDir,"HeatMap InducedvsWT.pdf"))
heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap Induced.vs.WT FC>=3",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=c(rep("red",4),rep("blue",4)),
          tracecol=NULL,
          srtCol=30)
dev.off()

#---------------------------------------------------------------------------------------------
### DATA ANNOTATION
#---------------------------------------------------------------------------------------------

myProbes <- rownames(exprs(eset))
head(myProbes)

## ----mappings------------------------------------------------------------
require(mogene11sttranscriptcluster.db)
keytypes(mogene11sttranscriptcluster.db)

## ------------------------------------------------------------------------
geneAnots <- select(mogene11sttranscriptcluster.db, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
geneAnots <- geneAnots[!is.na(geneAnots$SYMBOL),]
head(geneAnots)

##-------------------------------------------------------------------------
selected<- topTab[,"adj.P.Val"]<0.05 & abs(topTab[,"logFC"]) > 1
sum(selected)
selectedTopTab <- topTab[selected,]
head(selectedTopTab)

##------------------------------------------------------------------------
selectedProbes <- rownames(selectedTopTab)
selectedAnots <-  select(mogene11sttranscriptcluster.db, selectedProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
write.csv2(selectedAnots, file="selectedTopTab_Annotations.csv2")

##------------------------------------------------------------------------
# GOAnalysis

selectedEntrezs <-  select(mogene11sttranscriptcluster.db, selectedProbes, "ENTREZID")[,"ENTREZID"]
selectedEntrezs <- unique(selectedEntrezs[!is.na(selectedEntrezs)])

entrezUniverse = select(mogene11sttranscriptcluster.db, myProbes, "ENTREZID",  keytype="PROBEID")[,"ENTREZID"]
entrezUniverse <- unique(entrezUniverse[!is.na(entrezUniverse)])

require(GOstats)

## Creamos los "hiperparametros" en que se basa el analisis
GOparams = new("GOHyperGParams",
               geneIds=selectedEntrezs, 
               universeGeneIds=entrezUniverse,
               annotation="mogene11sttranscriptcluster.db", 
               ontology="BP",
               pvalueCutoff=0.001, conditional=FALSE,
               testDirection="over")

## Ejecutamos los analisis

GOhyper = hyperGTest(GOparams)
print(head(summary(GOhyper)))

## ------------------------------------------------------------------------
GOfilename =file.path(paste("GOResults.",".html", sep=""))
htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))

#---------------------------------------------------------------------------------------------
#END OF SCRIPT
#---------------------------------------------------------------------------------------------






