library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(DMRcate)


#library(doParallel)
#registerDoParallel(cores = 3)

# ----- working directory and path to folder ----
setwd("C:/Users/MARINA/Documents/MASTER/EPIGENOMICS/Ilumina_data_analysis/idats")

idat.folder <- "/Users/MARINA/Documents/MASTER/EPIGENOMICS/Ilumina_data_analysis/idats" 

targets <- read.metharray.sheet(base=idat.folder)
targets$Basename<-gsub("Class/","",targets$Basename)

# ---- loading data ----
rgset <- read.metharray.exp(targets = targets,verbose = T)
rgset #information about the array
#RGChannelSet object where we have the RAW DATA, NOT NORMALIZED

#The class of RGSet is a RGChannelSet object. 
#This is the initial object of a minfi analysis that contains the raw intensities in the green and red channels. 
#Note that this object contains the intensities of the internal control probes as well. 


# --- Data exploration ----
# Extract phenodata
phenoData <- pData(rgset)

#The RGChannelSet stores also a manifest object that contains the probe design information of the array:
manifest <- getManifest(rgset)
manifest


#Probes information
probesI<-getProbeInfo(manifest, type = "I")
probesII<<-getProbeInfo(manifest, type = "II")

##A MethylSet objects contains only the methylated and unmethylated signals
MSet <- preprocessRaw(rgset) 

Meth<-getMeth(MSet)
UnMeth<-getUnmeth(MSet)


# --- Get BETA values ---

#A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 
#An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. 
#Mapping a MethylSet to a RatioSet may be irreversible, i.e. one cannot be guranteed to retrieve the methylated and unmethylated signals from a RatioSet.

#A RatioSet can be created with the function ratioConvert:
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
#The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.
beta <- getBeta(RSet) 

# Beta values are not normalized yet, this is only the storage of beta values. 


# --GenomicRatioSet--Map to the genome

#The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation information.
#The output object is a GenomicRatioSet (class holding M or/and Beta values together with associated genomic coordinates). 
#It is possible to merge the manifest object with the genomic locations by setting the option mergeManifest to TRUE.

GRset <- mapToGenome(RSet)
#Get the probe location to the genome

beta <- getBeta(GRset)
M <- getM(GRset) 
CN <- getCN(GRset) #get the copy number

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


#Type I have 4 position (more than type 2). 
#we have the 2 adress for type I, and one for typeII
names(annotation)

annotation[2,]

levels(as.factor(unlist(strsplit(annotation$UCSC_RefGene_Group,";"))))

####################################

# ------ Quality control -----------

#minfi provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels.
#When plotting these two medians against each other, it has been observed that good samples cluster together, 
#while failed samples tend to separate and have lower median intensities. 
#In order to obtain the methylated and unmethylated signals, we need to convert the RGChannelSet to an object containing the methylated and unmethylated signals using the function preprocessRaw.
#It takes as input a RGChannelSet and converts the red and green intensities to methylated and unmethylated signals according to the special 450K probe design, and returns the converted signals in a new object of class MethylSet. It does not perform any normalization.

#getMeth() and getUnmeth() are used to obtain the intensities signal matrices:
head(getMeth(MSet))
head(getUnmeth(MSet))

#The functions getQC and plotQC are designed to extract and plot the quality control information from the MethylSet:
qc <- getQC(MSet)
qc #median of methylated and median of unmethylated

png("QC.plot2.png")
plotQC(qc)
# if the dots are below the dash, it is bad quality, and are bad hibridized.
#In this case, as the dots are up the line, it is good quality

plotQC(qc, badSampleCutoff = 10) #to change the cut off

dev.off()

#To further explore the quality of the samples, it is useful to look at the Beta value densities of the samples, with the option to color the densities by group:
densityPlot(MSet, sampGroups = pData(MSet)$Status)
phenoData$Status
#each line is a sample (6 for cancer, and 6 for normal), and we don't have any problem with the probes. 
# Si hi hagues un problema amb una probe, tindriem un peak al mig del grafic, per exemple, i llavros s'ha d'eliminar la probe. 

#The 450k array contains several internal control probes that can be used to assess the quality control of different sample 
#preparation steps (bisulfite conversion, hybridization, etc.). 
#The values of these control probes are stored in the initial RGChannelSet 
#and can be plotted by using the function controlStripPlot and by specifying the control probe type:
png("Bis_con.png")
controlStripPlot(rgset, controls="BISULFITE CONVERSION II")
dev.off()
#This plot shows the hibridization in the red and green channels. 


# ---- Sex prediction of the samples -----
#In the sample sheet we don't have the sex information, and may be interested to predict it. 

#By looking at the median total intensity of the X chromosome-mapped probes, denoted med(X)med(X), 
#and the median total intensity of the Y-chromosome-mapped probes, denoted med(Y)med(Y),
#one can observe two different clusters of points corresponding to which gender the samples belong to.
#To predict the gender, minfi separates the points by using a cutoff on log2med(Y)

#The default cutoff is ???2???2. Since the algorithm needs to map probes to the X-chr and to the Y-chr, the input of the function 
#getSex needs to be a GenomicMethylSet or a GenomicRatioSet.

predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
predictedSex <- getSex(GRset, cutoff = -2)

head(predictedSex)
#To choose the cutoff to separate the two gender clusters, 
#one can plot med(Y)med(Y) against med(Y)med(Y) with the function plotSex:

GRset<-addSex(GRset)####add sex to GRatio object

pdf("sex.pdf")
plotSex(GRset)
dev.off()

# --- continue with the quality control ----
# detection pvals

detP<-detectionP(rgset)
dim(detP)
head(detP) #colums are samples and probes are rows
#each prove(rows) has a detection p-value to know if the probe is well hybridized or not

#Plotting the mean detection p-value for each sample will allow us to gauge whether any samples have many failed probes 
#- this will be indicated by a large mean detection p-value. Samples with mean detection p-values exceeding a cutoff 
#such as 0.05 can be excluded from further analysis.
pdf("pvals.detect_1.pdf")

barplot(colMeans(detP),col=factor(targets$Status),las=2,cex.names=0.8,
        main="Mean detection p-values")
dev.off()
#we are plotting the mean of the detection p-value of each sample
#the threshold is 0.05, so if the detection p-value is over 0.05, the sample has to be removed, because it will contain probes that are not weel hybridized. 
# In this case, the mean of the detection p-values for each sample is under the threshold. It doesn't mean that all the probes are 
#well hibridized, we can have may be a probe that its detection p-value is higher than the threshold, but the average of all the probes for the sample
# are well hibridized, so we do not remove this sample. But mainly all the probes of the sample are well hibridized


###in case you want to remove samples

keep<-colMeans(detP) < 0.01 #keep the columns whose colmeans is less than the threshold
rgset<-rgset[,keep] #remove the ones we do not keep
targets = targets[keep,]
detP <- detP[,keep]

### to remove probe

keep<-rowMeans(detP) < 0.01
rgset<-rgset[keep,]

######################################

# ------------ Normalization ------------------

gRatioSet.quantile <- preprocessQuantile(rgset,verbose = T)##SQN #normalization using quantiles #1.mapping to the genome 2.fixing outliesrs. 3. quantile normalization
gRatioSet.Illumina <- preprocessIllumina(rgset)##Illumina
gRatioSet.SWAN <- preprocessSWAN(rgset,verbose = T)##SWAN
gRatio.Fun<-preprocessFunnorm(rgset,verbose = T)

#Which is the best method?
#Depending on what we want to study. Now we will use ILLUMINA

densityPlot(getBeta(gRatioSet.quantile), sampGroups = targets$Status)


####correlations
sample.quantile<-colMeans(getBeta(gRatioSet.quantile),na.rm = T)
sample.illumina<-colMeans(getBeta(gRatioSet.Illumina),na.rm = T)
sample.swan<-colMeans(getBeta(gRatioSet.SWAN),na.rm = T)


mod1 <-lm(as.numeric(sample.quantile)~as.numeric(sample.illumina))## 
modsum<-summary(mod1)

r2 <- cor.test(as.numeric(sample.quantile),as.numeric(sample.illumina),method="spearman")$estimate
my.p<-cor.test(as.numeric(sample.quantile),as.numeric(sample.illumina),method="spearman")$p.value
my.p<-signif(my.p, digits=3)

#png("Quantiles_Illumina.png")
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
plot(as.numeric(sample.quantile)~as.numeric(sample.illumina), main=paste("Correlation", "p-value", my.p, sep=" "), pch=20,col="grey40")

#    ,xlim=c(0,10000),ylim=c(0,10000))
#  plot(as.numeric(cnts[,1])~as.numeric(cnts[,2]), main=paste("Correlation","Pulse Replicates Global", "p-value", my.p, sep=" "), xlab="Ctrl1", ylab="Ctrl2", pch=20)
abline(mod1, col="red")
# abline(h=0,col="blue1")
#abline(v=0,col="blue1")
legend('topleft', legend = mylabel, bty = 'n')


dev.off()


#####

ann450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
# for each probe (rows) we get the genomic position, the probe sequence, the name..., the type..
#one probe can match to several positions (ex/two genes)

ann450k<-as.data.frame(ann450k)

#remove CpG that are in the sexual chromosomes, because methylations in these chr are very variable and come from impinting, which is not important for the analysis. 
keep <- !(featureNames(gRatioSet.quantile) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep) #nearly 11000 probes are located in the chr X and Y

gRatioSet.quantile <- gRatioSet.quantile[keep,]

## remove probes with SNPs at CpG or SBE site
#####get snp info
#snp.info<-getSnpInfo(gRatioSet.quantile)
gRatioSet.quantile<- dropLociWithSnps(gRatioSet.quantile)


####get values to proceed with analysis
betas<-getBeta(gRatioSet.quantile) #betas are normalized
M <- getM(gRatioSet.quantile)


###put same order

betas<-betas[,match(rownames(pData(gRatioSet.quantile)),colnames(betas))]
#####

colnames(betas)<-pData(gRatioSet.quantile)$Name
#betas<-DMRcate::rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)

####remove probes

#probes2remove<-read.csv("48639-non-specific-probes-Illumina450k.txt")

#betas<-betas[!(rownames(betas) %in% probes2remove$TargetID), ]


########################################################

# ------------ Identify Differential methylation positions -----------
##DMP finding
# 3 things: betas, phenotype (the variable we are considering (cancer vs normal)) and its type
# be sure that the order of betas is the same as the pheno
colnames(betas)
targets$Name #pheno
dmp <- dmpFinder(betas, pheno = targets$Status, type = "categorical")
dmp<-subset(dmp,dmp$qval<.05)
head(dmp)

dim(dmp)

####LIMMA
library(limma) #used to build a linear model

group<-as.factor(pData(gRatioSet.quantile)$Status)
group<-relevel(group,ref="Normal")

design <- model.matrix(~group+pData(gRatioSet.quantile)$predictedSex)##Cancer up
design <- model.matrix(~pData(gRatioSet.quantile)$Status)##
colnames(design)<-c("Intercept","Comp")
####generate adjusted beta matrix
fit <- lmFit(logit2(betas), design)

fit2 <- eBayes(fit)

results <- topTable(fit2,coef="Comp",num=dim(fit2)[1],sort.by="P",adjust.method = "BH") ####Pval list 1

res<-subset(results,results$adj.P.Val<.05)
################
####

dmp <- merge(res,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
rownames(dmp)<-dmp$Row.names
dmp<-dmp[,-1]
head(dmp)

## Outputs
#subset DM positions 
proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
dim(proms) #from the 82000 CpG we have nearly 27000 CpG in promoters

#now, we want to know which of these are in Island
#proms.island<-proms[grep("Island",proms$Relation_to_Island),]
proms.island<-proms[proms$Relation_to_Island=="Island",]
dim(proms.island)#only 9000
head(proms.island)#now we can merge with the betas

#proms.island<-merge(proms.island,betas,by="row.names")
#rownames(proms.island)<-proms.island$Row.names
#proms.island<-proms.island[,-1]

write.csv(dmp,file="PT.vs.T.csv",row.names=T)

#DM positions can split the samples in normal and in cancer?
head(targets)#cancer has a name with T at the end, and normalt have a P
Normal <- grep("P",colnames(betas),perl=T)
Cancer <- grep("T",colnames(betas),perl=T)

spcol <- c(rep("blue1",length(Normal)),rep("red",length(Cancer)))

library(gplots)

#We want a matrix of those CpG in betas and those that are in proms.islands results, and first all the norlmals and then all the cancers
l <- as.matrix(betas[rownames(betas) %in% rownames(proms.island),c(Normal,Cancer)])
head(l)

heatmap.2(l,
          main="Diffmeth CpG's",
          labRow=NA,
          trace="none",
          na.rm=T,
          col=greenred,
          ColSideColors=spcol,
          distfun=function(x) dist(x,method="euclidean"),
          dendrogram = "column")



pca <- prcomp(t(na.omit(l)))
plot(pca$x[,1],pca$x[,2],col=spcol,cex=.2)
text(pca$x[,1],pca$x[,2],labels=rownames(pca$x),col=spcol,cex=1)

dev.off()

# -------- testing candidates -------

#select GFRA1 and GSTM2,SEPT9
tested<-proms.island[grep("GFRA1|GSTM2|SEPT9",proms.island$UCSC_RefGene_Name),]
head(tested)
dim(tested)# 35 CpG with these genes selected

Normal <- grep("P",colnames(betas),perl=T)
Cancer <- grep("T",colnames(betas),perl=T)

spcol <- c(rep("blue1",length(Normal)),rep("red",length(Cancer)))

l <- as.matrix(betas[rownames(betas) %in% rownames(tested),c(Normal,Cancer)])

heatmap.2(l,
          main="Diffmeth CpG's biomarkers",
          labRow=NA,
          trace="none",
          na.rm=T,
          col=greenred,
          ColSideColors=spcol,
          distfun=function(x) dist(x,method="manhattan"),
          dendrogram = "column")

#With only these 3 genes, the differenciated methylated positions are more dramatic.
#The normal samples are low methylated while the cancer samples are highly methlyated. So these genes could be biomarkers. 



###############################################
# ------- GO enrichment ---------

#####GO enrichment analysis
genesid<-proms.island$UCSC_RefGene_Name

###split
#This is because one CpG can map to different isoforms of same genes or overlapping genes (two positions of different genes)--- we can not deal with this type of name vector
#So, remove the ; and get the uniques genes_id
genesid<- strsplit(as.character(genesid),';')
genesid<-unique(unlist(genesid))

##remove empty elements from a vector
genesid<-genesid[ genesid != "" ]
genesid <- genesid[!is.na(genesid)]

#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)


eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego2 <- enrichGO(gene         = eg$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable=T
)



dotplot(ego2, title="Control vs Cancer",showCategory=35)
dev.off()


ego2 <-enrichKEGG(
  eg$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  #universe,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)


my.go.res<-as.data.frame(ego2)


#####DOSE (Disease ontology database)

library(DOSE)
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
#the bigger the cercle, the more genes. The more red color, more significally statistical. 


########missMethyl
library(missMethyl)

enrichment_GO <- gometh(rownames(proms.island),all.cpg = rownames(betas),collection = "GO", 
                        array.type = "450K",plot.bias = T,prior.prob = T,equiv.cpg = T,anno = ann450k) 

enrichment_GO<-enrichment_GO[enrichment_GO$ONTOLOGY=="BP",]

enrichment_GO<-enrichment_GO[enrichment_GO$FDR<.05,]

enrichment_KEGG <- gometh(rownames(proms.island),all.cpg = rownames(betas),collection = "KEGG", 
                          array.type = "450K",plot.bias = T,prior.prob = T,equiv.cpg = T,anno = ann450k) 

enrichment_KEGG<-enrichment_KEGG[enrichment_KEGG$FDR<.05,]
enrichment_KEGG[enrichment_KEGG$P.DE<.05,]



###########################################

# ---------- DM regions ---------------


#GFRA1 and GSTM2,SEPT9
###DMR

#####bumphunter----consume lot of memory---not use in a laptop!!!
#library(bumphunter)
#pheno <- pData(gRatioSet.quantile)$Status
#designMatrix <- model.matrix(~ pheno)
#Run the algorithm with a large number of permutations, say B=1000:

# dmrs <- bumphunter(gRatioSet.quantile, design = designMatrix, cutoff = 0.01, B=1000, type="Beta")

# library(doParallel)
#registerDoParallel(cores = 3)
#The results of bumphunter are stored in a data frame with the rows being the different differentially methylated regions (DMRs):
#  
#names(dmrs)
#head(dmrs$table, n=3)
#-------


####DMRCate procedure
#dmp[grep("GFRA1",dmp$UCSC_RefGene_Name),]

library(DMRcate)

data(dmrcatedata)

myMs<-getM(gRatioSet.quantile)
#myBetas<-getBeta(gRatioSet.quantile)
colnames(myMs)<-targets$Name
#myMs <- logit2(myBetas)###convert back the Betas to M

myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05,rmXY = T,rmcrosshyb = TRUE) ### get out SNPs of my samples
#myBetas <- rmSNPandCH(myBetas, dist=2, mafcut=0.05) ### get out SNPs of my samples
###establish the groups

Normal <- grep("P",colnames(myMs.noSNPs),perl=T)
Cancer <- grep("T",colnames(myMs.noSNPs),perl=T)

pheno <- pData(gRatioSet.quantile)$Status
designMatrix <- model.matrix(~ pheno)

myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "450K",
                             analysis.type="differential", design=designMatrix, coef=2)

head(myannotation@ranges@elementMetadata@listData$stat)

str(myannotation)


###DMR finding
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

str(dmrcoutput)
my.dmrs<-dmrcoutput$results

##filter by Stouffervpval

#my.dmrs<-my.dmrs[my.dmrs$Stouffer<.05,]
#my.dmrs<-my.dmrs[my.dmrs$no.cpgs >= 5,]

results.ranges <- extractRanges(dmrcoutput, genome = "hg19")


dmrs.res<-as.data.frame(results.ranges)
###filter those with less than 5 CpGs inside DMR
my.dmrs<-subset(dmrs.res,dmrs.res$no.cpgs >= 5)


write.csv(my.dmrs,"dmrs.proms.csv")

###select specifical

my.dmrs[grep("GFRA1",my.dmrs$overlapping.promoters),]

#GFRA1 and GSTM2,SEPT9
#GSTM2

groups <- c(Cancer="red", Normal="blue1")

cols <- groups[as.character(pheno)]

DMR.plot(ranges=results.ranges, dmr=1, CpGs=myMs.noSNPs,what="M",
         arraytype = "450K", phen.col=cols, genome="hg19")

DMR.plot(ranges=results.ranges, dmr=598, CpGs=myBetas, phen.col=cols, genome="hg19")
dev.off()

DMR.plot(ranges=results.ranges, dmr=598, CpGs=myBetas, what="Beta", arraytype = "450K",
         phen.col=cols, genome="hg19")


################PLOTS

####DMRff

library(dmrff)
library(limma)
group<-as.factor(pData(gRatioSet.quantile)$Status)
group<-relevel(group,ref="Normal")

design <- model.matrix(~group)##Cancer up
colnames(design)<-c("Intercept","Comp")
####generate adjusted beta matrix
fit <- lmFit(betas, design)

fit2 <- eBayes(fit)


stats <- data.frame(estimate=fit2$coefficients[,"Comp"],
                    se=sqrt(fit2$s2.post) * fit2$stdev.unscaled[,"Comp"],
                    p.value=fit2$p.value[,"Comp"])

####add annotation to stats
stats<-merge(stats,ann450k[,c("chr","pos")],by="row.names")

dmrs <- dmrff(estimate=stats$estimate,
              se=stats$se,
              p.value=stats$p.value,
              methylation=betas,
              chr=stats$chr,
              pos=stats$pos,
              maxgap=500,
              verbose=T)


dmrs <- dmrs[which(dmrs$p.value < 0.05 & dmrs$n > 1),]
annot<-dmrff.sites(dmrs, stats$chr, stats$pos)
library(GenomicRanges)

dmrs<-makeGRangesFromDataFrame(dmrs)

#####
###expression
####
expres<-read.csv("../Expression_data.csv")


genesid<-proms.island$UCSC_RefGene_Name

###split

genesid<- strsplit(as.character(genesid),';')
genesid<-unique(unlist(genesid))

genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
genesid <- genesid[!is.na(genesid)]


###select the genes of the proms+island
expres<-expres[expres$Gene %in% genesid,]

##mean
expres<-aggregate(.~Gene, data=expres[,c(2:13,15)], mean)
rownames(expres)<-expres$Gene
expres<-expres[,-1]


Normal <- grep("\\d+P",colnames(expres),perl=T)
Cancer <- grep("\\d+T",colnames(expres),perl=T)

pdf("ANOVA_PLOT_expression_signature_proms.cgi.pdf")
#GFRA1 and GSTM2,SEPT9

i<-"SEPT9"

for (i in 1:(dim(expres)[1])){
  
  Norm<-as.numeric(expres[i,Normal])
  Canc<-as.numeric(expres[i,Cancer])
  
  class<-c(rep("Norm",length(Norm)), rep("Canc",length(Canc)))
  y<-c(Norm,Canc)
  
  result<-aov(y~as.factor(class))
  
  res<-summary(result)
  
  pval<-res[[1]][["Pr(>F)"]][1]
  pval<-signif(pval, digits=3)
  
  legend=c("Normal","Cancer")
  
  boxplot(Norm,Canc,names=legend, col="lightgrey",main=paste(rownames(expres)[i],"p-value",pval, sep=" "), ylab="Expression")
}

dev.off()



###methylation

##add beta values to proms+islands results

proms.island<-merge(proms.island,betas,by="row.names")
rownames(proms.island)<-proms.island$Row.names
proms.island<-proms.island[,-1]



Normal <- grep("\\d+P",colnames(proms.island),perl=T)
Cancer <- grep("\\d+T",colnames(proms.island),perl=T)

pdf("ANOVA_PLOT_methylation_signature_proms.cgi.pdf")

for (i in 1:(dim(proms.island)[1])){
  
  Norm<-as.numeric(proms.island[i,Normal])
  Canc<-as.numeric(proms.island[i,Cancer])
  
  class<-c(rep("Norm",length(Norm)), rep("Canc",length(Canc)))
  y<-c(Norm,Canc)
  
  result<-aov(y~as.factor(class))
  
  res<-summary(result)
  
  pval<-res[[1]][["Pr(>F)"]][1]
  pval<-signif(pval, digits=3)
  
  legend=c("Normal","Cancer")
  
  boxplot(Norm,Canc,names=legend, col="lightgrey",main=paste(rownames(proms.island)[i],proms.island$UCSC_RefGene_Name[i],"p-value",pval, sep=" "), ylab="Methylation")
}

dev.off()