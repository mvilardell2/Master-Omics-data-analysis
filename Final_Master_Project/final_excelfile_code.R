###################################################################
###########      GENERATION OF THE GTF FILE       #################
##################################################################

setwd("C:/Users/MARINA/Documents/MASTER/TFM/Data_analysis")

# This gtf file will be the annotation file for the generation of the excel

library(stringr)

geneList<-rownames(brca)

chrlist<-list()
strand2list<-list()
startlist<-list()
endlist<-list()
lenlist<-list()


for (i in 1:nrow(brca)){ #nrow(brca)
  genes<-geneList[i]
  test<-annotated_probes_TF[str_detect(annotated_probes_TF$`Gene Symbol`,genes),]
  if(nrow(test)>1){
    if(test$Alignments[1]=='---'){
      gene_line<-test[2,]
    }else{
      gene_line<-test[1,]
    }
  }else{
    gene_line<-test
  }
  full<-str_split(gene_line$Alignments,'//')[[1]][1]
  
  chr<-str_split(full,':')[[1]][1]
  chrlist<-append(chrlist,chr)
  
  resta<-str_split(full,':')[[1]][2]
  coordinates<-str_split(resta,' ')[[1]][1]
  
  strand<-str_split(resta,' ')[[1]][2]
  stran2<-substring(strand,2,2)
  strand2list<-append(strand2list,stran2)
  
  
  start<-str_split(coordinates,'-')[[1]][1]
  startlist<-append(startlist,start)
  
  end<-str_split(coordinates,'-')[[1]][2]
  endlist<-append(endlist,end)
  
  len<-as.numeric(end)-as.numeric(start)
  lenlist<-append(lenlist,len)
  
}

geneListAnnot<-data.frame(Geneid=geneList,
                          Chr=unlist(chrlist),
                          start=unlist(startlist),
                          End=unlist(endlist),
                          strand=unlist(strand2list),
                          length=unlist(lenlist))
rownames(geneListAnnot)<-geneList


##########################################################################
############  GENERATION OF THE EXCEL FOR EACH COMPARISON     ############
##########################################################################

### BRCA1.SA vs NOMUT.SA

# exprMat=exprs(brca_smk_sa1)

# Fit the model 
exprMat_brca1_sa=exprs(brca_smk_sa1)
fit <- lmFit(exprMat_brca1_sa, design) 

# --- Contrasts ---
contrast.matrix <- makeContrasts(cond1=BRCA1.SA-NOMUT.SA,
                                 levels = design)

# --- Fit contrasts and Bayesian adjustment ---
fit2 <- contrasts.fit(fit, contrast.matrix)
fite <- eBayes(fit2)

fit.main=fite
cond<-as.factor(pData(brca_smk_sa1)$Phenotype) #Phenotype
contrast<-list(c('BRCA1.SA','NOMUT.SA'))
BasicFunctions::RNAseq.resAnnot(exprMat_brca1_sa,geneListAnnot,cond,fitMain = fit.main, contrast=contrast,resultsDir = getwd(),resAnnotFilename = "resultsAnnot_brca1_sa_", Excel = TRUE,excelFilename = 'excelfile_brca1_sa_',GO=FALSE,pvalue=0.05,logFC = 1)

#------------------------------------------------------------------------

### BRCA1.AF vs NOMUT.AF

# Fit the model 
exprMat_brca1_af=exprs(brca_smk_af1)
fit <- lmFit(exprMat_brca1_af, design) 

# --- Contrasts ---
contrast.matrix <- makeContrasts(cond1=BRCA1.AF-NOMUT.AF,
                                 levels = design)

# --- Fit contrasts and Bayesian adjustment ---
fit2 <- contrasts.fit(fit, contrast.matrix)
fite <- eBayes(fit2)

fit.main=fite
cond<-as.factor(pData(brca_smk_af1)$Phenotype) #Phenotype
contrast<-list(c('BRCA1.AF','NOMUT.AF'))
BasicFunctions::RNAseq.resAnnot(exprMat_brca1_af,geneListAnnot,cond,fitMain = fit.main, contrast=contrast,resultsDir = getwd(),resAnnotFilename = "resultsAnnot_brca1_af", Excel = TRUE,excelFilename = 'excelfile_brca1_af',GO=FALSE)

#-----------------------------------------------------------------------------

### BRCA2.SA vs NOMUT.SA

# Fit the model 
exprMat_brca2_sa=exprs(brca_smk_sa2)
fit <- lmFit(exprMat_brca2_sa, design) 

# --- Contrasts ---
contrast.matrix <- makeContrasts(cond1=BRCA2.SA-NOMUT.SA,
                                 levels = design)

# --- Fit contrasts and Bayesian adjustment ---
fit2 <- contrasts.fit(fit, contrast.matrix)
fite <- eBayes(fit2)

fit.main=fite
cond<-as.factor(pData(brca_smk_sa2)$Phenotype) #Phenotype
contrast<-list(c('BRCA2.SA','NOMUT.SA'))
BasicFunctions::RNAseq.resAnnot(exprMat_brca2_sa,geneListAnnot,cond,fitMain = fit.main, contrast=contrast,resultsDir = getwd(),resAnnotFilename = "resultsAnnot_brca2_sa", Excel = TRUE,excelFilename = 'excelfile_brca2_sa',GO=FALSE)

#------------------------------------------------------------------------

### BRCA2.AF vs NOMUT.AF

# Fit the model 
exprMat_brca2_af=exprs(brca_smk_brca2.af)
fit <- lmFit(exprMat_brca2_af, design) 

# --- Contrasts ---
contrast.matrix <- makeContrasts(cond1=BRCA2.AF-NOMUT.AF,
                                 levels = design)

# --- Fit contrasts and Bayesian adjustment ---
fit2 <- contrasts.fit(fit, contrast.matrix)
fite <- eBayes(fit2)

fit.main=fite
cond<-as.factor(pData(brca_smk_brca2.af)$Phenotype) #Phenotype
contrast<-list(c('BRCA2.AF','NOMUT.AF'))
BasicFunctions::RNAseq.resAnnot(exprMat_brca2_af,geneListAnnot,cond,fitMain = fit.main, contrast=contrast,resultsDir = getwd(),resAnnotFilename = "resultsAnnot_brca2_af", Excel = TRUE,excelFilename = 'excelfile_brca2_af',GO=FALSE)
















