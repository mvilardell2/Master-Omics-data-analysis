
# ------------ Libraries ------------
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(DMRcate)
library(FactoMineR)

#library(doParallel)
#registerDoParallel(cores = 3)

# ----- working directory and path to folder ----
setwd("C:/Users/MARINA/Documents/MASTER/EPIGENOMICS/Ilumina_data_analysis/Practical")
idat.folder <- "/Users/MARINA/Documents/MASTER/EPIGENOMICS/Ilumina_data_analysis/Practical"

targets <- read.metharray.sheet(base=idat.folder)


# ---- loading data ----
rgset <- read.metharray.exp(targets = targets,verbose = T)
rgset #RGChannelSet object where we have the RAW DATA
# EPIC array
##The class of RGSet is a RGChannelSet object. 
#This is the initial object of a minfi analysis that contains the raw intensities in the green and red channels. 
#Note that this object contains the intensities of the internal control probes as well. 

# --- Data exploration ----
# Extract phenodata
phenoData <- pData(rgset)

#The RGChannelSet stores also a manifest object that contains the probe design information of the array:

manifest <- getManifest(rgset)
#EPIC array and type of probes


####Probes information
getProbeInfo(manifest, type = "I")
getProbeInfo(manifest, type = "II")

##A MethylSet objects contains only the methylated and unmethylated signals
MSet <- preprocessRaw(rgset) #13 samples 866836 probes

Meth<-getMeth(MSet)
dim(Meth)
UnMeth<-getUnmeth(MSet)
dim(UnMeth)

# --- Get BETA values ---

#A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 
#An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. 
#Mapping a MethylSet to a RatioSet may be irreversible, i.e. one cannot be guranteed to retrieve the methylated and unmethylated signals from a RatioSet.

#A RatioSet can be created with the function ratioConvert:
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
#The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.
beta <- getBeta(RSet) 
head(beta)
# Beta values are not normalized yet, this is only the storage of beta values. 


# --GenomicRatioSet--Map to the genome

#The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation information.
#The output object is a GenomicRatioSet (class holding M or/and Beta values together with associated genomic coordinates). 
#It is possible to merge the manifest object with the genomic locations by setting the option mergeManifest to TRUE.

GRset <- mapToGenome(RSet)
#Get the probe location to the genome

beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)

sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)

gr <- granges(GRset) #assign the position of the probe to the genome
head(gr, n= 3)


# ---- Full annotation for each probe ---
annotation <- getAnnotation(GRset)
head(annotation)
#If we have type I--> we will have 2 sequences (unmethylated and methylated)
#If it is type II--> only 1 sequence

levels(as.factor(unlist(strsplit(annotation$UCSC_RefGene_Group,";"))))



####################################

# ------ Quality control -----------

#minfi provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels.
#When plotting these two medians against each other, it has been observed that good samples cluster together, 
#while failed samples tend to separate and have lower median intensities. 
#In order to obtain the methylated and unmethylated signals, we need to convert the RGChannelSet to an object containing the methylated and unmethylated signals using the function preprocessRaw.
#It takes as input a RGChannelSet and converts the red and green intensities to methylated and unmethylated signals according to the special 450K probe design, and returns the converted signals in a new object of class MethylSet. It does not perform any normalization.

#getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:
head(getMeth(MSet))
head(getUnmeth(MSet))

#The functions getQC and plotQC are designed to extract and plot the quality control information from the MethylSet:
qc <- getQC(MSet)
qc #median of methylated and median of unmethylated

plotQC(qc)
plotQC(qc, badSampleCutoff = 10) #to change the cut off
# There are samples that are below the threshold, which means that some problem may have occured with the hibridation
#and this samples have not quite good quality. 


#samples to remove: 13,5,12
remove<-c(13,5,12)
rgset<-rgset[,-remove]
targets = targets[-remove,]
#keep 10 samples
dim(rgset)

#To further explore the quality of the samples, it is useful to look at the Beta value densities of the samples, with the option to color the densities by group:
densityPlot(MSet, sampGroups = targets$patient.diagnosis.ch1)
#each line of the plot represents a sample (6 for Chron disease and 4 for healthy), and they are colored depending on their phenotype.
#The samples have a good shape/behaviour. They follow a B-modal shape.


#detection pvals
detP<-detectionP(rgset)
dim(detP) #all the probes in the EPIC
head(detP) #colums are samples and probes are rows
#for each prove(rows) and sample there is a detection p-value to know if the probe is well hybridized or not


barplot(colMeans(detP),col=factor(targets$patient.diagnosis.ch1),las=2,cex.names=0.8,
        main="Mean detection p-values")


## No samples exceed the cutoff of 0.05, so now no samples are going to be removed.
keep<-colMeans(detP) < 0.01
table(keep)

### to remove probe
#if the row mean of the detection p-value for a probe is higher than the threshold, it means that the probe has 
# high p-values in all the samples, hence the probe is not wee hibridized and has to be removed.
keep<-rowMeans(detP) < 0.01
table(keep)
#keep with 861342 probes

#keep with the well hibrizided probes in the row of the rgset
rgset<-rgset[keep,]



# ---- Sex prediction of the samples -----

##A MethylSet objects contains only the methylated and unmethylated signals
MSet <- preprocessRaw(rgset) #9 
#A RatioSet can be created with the function ratioConvert:
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)

predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
predictedSex <- getSex(GRset, cutoff = -2)


head(predictedSex)
#To choose the cutoff to separate the two gender clusters, 
#one can plot med(Y)med(Y) against med(Y)med(Y) with the function plotSex:

GRset<-addSex(GRset)####add sex to GRatio object

pdf("sex.pdf")
plotSex(GRset)
dev.off()






#####analysis

###Normalization

gRatioSet.quantile <- preprocessQuantile(rgset,verbose = T)##SQN


densityPlot(getBeta(gRatioSet.quantile), sampGroups = targets$patient.diagnosis.ch1)




#####Annotation

ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k<-as.data.frame(ann850k)


## remove probes with SNPs at CpG or SBE site
#####get snp info
snp.info<-getSnpInfo(gRatioSet.quantile)
####mapping

gRatioSet.quantile<- dropLociWithSnps(gRatioSet.quantile)


####get values to proceed with analysis
betas<-getBeta(gRatioSet.quantile)



#####PCA
spcol <- c(rep("blue1",4),rep("red",4))

pca <- prcomp(t(na.omit(betas)))
plot(pca$x[,1],pca$x[,2],col=spcol,cex=.8)
text(pca$x[,1],pca$x[,2],labels=rownames(pca$x),col=spcol,cex=1)

pca <- FactoMineR::PCA(t(betas), scale.unit = T, graph = F, ncp = 40)
factoextra::fviz_pca_ind(pca, axes = c(1,2), habillage=as.factor(targets$Status), repel = F)

Mvalues<-logit2(betas)

##DMP finding

dmp <- dmpFinder(Mvalues, pheno = targets$patient.diagnosis.ch1, type = "categorical")

dmp<-subset(dmp,dmp$qval<.05)

library(limma)

group<-as.factor(pData(gRatioSet.quantile)$patient.diagnosis.ch1)###targets$Status
group<-relevel(group,ref="WT")

design <- model.matrix(~group)##KO up
colnames(design)<-c("Intercept","Comp")
####generate adjusted beta matrix
fit <- lmFit(Mvalues, design)

fit2 <- eBayes(fit)

results <- topTable(fit2,coef="Comp",num=dim(fit2)[1],sort.by="P",adjust.method = "BH") ####Pval list 1

res<-subset(results,results$adj.P.Val<.05)
##subset by log fold chahe
res<-res[res$logFC>2 | res$logFC<  -2, ]


dmp <- merge(res,ann850k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
dmp<-as.data.frame(dmp)
rownames(dmp)<-dmp$Row.names
dmp<-dmp[,-1]
## Outputs
proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
proms.island<-proms[proms$Relation_to_Island=="Island",]
#proms.island<-merge(proms.island,betas,by="row.names")
#rownames(proms.island)<-proms.island$Row.names
#proms.island<-proms.island[,-1]
write.csv(dmp,file="KO.vs.WT.csv",row.names=T)

colnames(betas)<-targets$Status

spcol <- c(rep("blue1",4),rep("red",4))

library(gplots)


heatmap.2(as.matrix(betas[rownames(betas) %in% rownames(proms.island),]),
          main="Diffmeth CpG's ",
          labRow=NA,
          trace="none",
          na.rm=T,
          col=greenred,
          ColSideColors=spcol,
          distfun=function(x) dist(x,method="euclidean"),
          dendrogram = "column")



#####
# install.packages("devtools")
#devtools::install_github("TomKellyGenetics/heatmap.2x", ref="test")

library(heatmap.2x)## approach
library(RColorBrewer)

samples<-ifelse(targets$Status=="KO","red","blue1")

#####M3C
my.matrix<-as.matrix(betas[rownames(betas) %in% rownames(proms.island),])
my.matrix<-as.matrix(betas[rownames(betas) %in% rownames(proms),])

data <- t(scale(t(my.matrix))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol <- c(rep("blue1",4),rep("red",4))

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="both", 
           cexRow=1, cexCol=1,
           main="KO vs WT",
           #labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)

legend("bottomleft",c("2","1"),pch=20:20,col=c("red","blue1"))







##GO enrichment
#####GO enrichment analysis

####we want to split our results
##in KO enriched and WT enriched


ko<-proms.island[proms.island$logFC>0,]
wt<-proms.island[proms.island$logFC<0,]

genesid<-ko$UCSC_RefGene_Name

genesid<- strsplit(as.character(genesid),';')
genesid<-unique(unlist(genesid))

genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
genesid <- genesid[!is.na(genesid)]
#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)


eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego2 <- enrichGO(gene         = eg$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.1,
                 qvalueCutoff  = 0.05,
                 readable=T
)

dotplot(ego2, title=" DMPs KO vs WT ",showCategory=25)

ego2$Description[grep("metabolism",ego2$Description)]
########missMethyl
library(missMethyl)

enrichment_GO <- gometh(rownames(dmp),all.cpg = rownames(betas),collection = "GO", 
                        array.type = "EPIC",plot.bias = T,prior.prob = T,equiv.cpg = T,anno = ann850k) 

enrichment_GO<-enrichment_GO[enrichment_GO$ONTOLOGY=="BP",]

enrichment_GO<-enrichment_GO[enrichment_GO$FDR<.05,]
enrichment_GO<-enrichment_GO[order(enrichment_GO$FDR,decreasing = F),]
enrichment_KEGG <- gometh(rownames(proms.island),all.cpg = rownames(betas),collection = "KEGG", 
                          array.type = "450K",plot.bias = T,prior.prob = T,equiv.cpg = T,anno = ann450k) 

enrichment_KEGG<-enrichment_KEGG[enrichment_KEGG$FDR<.05,]
enrichment_KEGG[enrichment_KEGG$P.DE<.05,]

#####DOSE

data(DO2EG)

eg.ids <- toTable(org.Hs.egALIAS2EG)

gene.set<-eg.ids[eg.ids$alias_symbol %in% genesid,"gene_id"]

univ <- Lkeys(org.Hs.egGO)
x <- enrichDO(gene          = gene.set,
              ont           = "DO", 
              pvalueCutoff  = 1,
              pAdjustMethod = "fdr",
              universe      = univ,
              qvalueCutoff  = .1,
              readable      = FALSE)
results.dose<-as.data.frame(x)


dotplot(x)
##myMs<-getM(gRatioSet.quantile)
#colnames(myMs)<-targets$Name
myMs <- logit2(betas)###convert back the Betas to M
colnames(myMs)<-targets$SampleName
###establish the groups


pheno <- pData(gRatioSet.quantile)$Status
designMatrix <- model.matrix(~ pheno)

myannotation <- cpg.annotate("array", myMs, what="M", arraytype = "EPIC",
                             analysis.type="differential", design=designMatrix, coef=2)
###DMR finding

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2,betacutoff = .1)

results.ranges <- extractRanges(dmrcoutput, genome = "hg19")


dmrs.res<-as.data.frame(results.ranges)
###filter those with less than 5 CpGs inside DMR
my.dmrs<-subset(dmrs.res,dmrs.res$no.cpgs >= 10)


write.csv(my.dmrs,"dmrs.proms.csv")

###select specifical

dmrs.res[grep("GFRA1",dmrs.res$overlapping.promoters),]

#GFRA1 and GSTM2,SEPT9
#GSTM2

groups <- c(WT="magenta", KO="forestgreen")
type<-pheno
cols <- groups[as.character(type)]

DMR.plot(ranges=results.ranges, dmr=9, CpGs=betas, phen.col=cols, genome="hg19")
dev.off()

DMR.plot(ranges=results.ranges, dmr=9, CpGs=myMs,what="M",
         arraytype = "EPIC", phen.col=cols, genome="hg19")

DMR.plot(ranges=results.ranges, dmr=598, CpGs=myBetas, what="Beta", arraytype = "450K",
         phen.col=cols, genome="hg19")
