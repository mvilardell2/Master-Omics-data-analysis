##################################################################
############# VIRAL PRACTICE STATISTICAL METHODS #################
##################################################################

setwd("~/MASTER/BIOINFORMATICS/STATISTICS")
viraldata<-read.table(file='viral.28.txt',header=T)
#head(viraldata)
#str(viraldata)
library(ggpubr)
library(kableExtra)
win.graph()
###########################################################################################

#1.Describe the main characteristics of the dataset: perform a univariate descriptive 
# analysis of the first 6 variables

#Transform variables to factor (infection, sind,gender, hosp)
viraldata$infection<-factor(viraldata$infection, levels=c(0,1),labels=c('bacterial infection','viral infection'))
viraldata$sind<-factor(viraldata$sind,levels = c(0,1),labels=c('symptoms remain','symptoms finished'))
viraldata$gender<-factor(viraldata$gender,levels = c(0,1),labels=c('female','male'))
viraldata$hosp<-factor(viraldata$hosp,levels = c(1,0),labels=c('hospitalization','no hospitalization'))

str(viraldata[,1:6])
summary(viraldata)
#### CATEGORICAL DATA

##INFECTION 
infection_abs_freq<-table(viraldata$infection) #Absolute frequency
infection_rel_freq<-prop.table(infection_abs_freq) #Relative frequency
freqtable_infection<-cbind(infection_abs_freq,infection_rel_freq) #Both frequencies
colnames(freqtable_infection)<-c('Absolute frequency','Relative frequency') #change name of colnames 

#Creation of a kable to build a stylish table. 
kable_infection<-kable(freqtable_infection, caption = ' Table X')%>%kable_styling(bootstrap_options = c('responsive'),full_width = F, html_font = "Cambria")%>%
  row_spec(0, color = "black", background = "gainsboro") 
kable_infection

# Infection data is equal distributed in each kind of infection, so in this study there is the same amount of bacterial and viral infected patients, concretely 70.


#plot

plot_inf<-ggbarplot(data.frame(infection_rel_freq),x='Var1',y='Freq',xlab='Kind of infection',
          ylab='Relative frequency',lab.nb.digits = 2, label = T, title = ' Figure 1. Representation of relatives frequencies',
          fill='Var1',legend='none',palette = c('skyblue1','royalblue3'),ylim=c(0,0.7),subtitle='A) Infection')



##SIND
sind_abs_freq<-table(viraldata$sind)
sind_rel_freq<-prop.table(sind_abs_freq)
freqtable_sind<-cbind(sind_abs_freq,sind_rel_freq)

colnames(freqtable_sind)<-c('Absolute frequency','Relative frequency')
kable_sind<-kable(round(freqtable_sind,2), caption = ' Table X')%>%kable_styling(bootstrap_options = c('responsive'),full_width = F, html_font = "Cambria")%>%
  row_spec(0, color = "black", background = "gainsboro")

# There are more patients that have still symptoms than patients that have no longer symptoms. Around 66% of the patients remain with symptoms and 34% of them finished with the symptoms. 


#plot
plot_sind<-ggbarplot(data.frame(sind_rel_freq),x='Var1',y='Freq',xlab=' ',
          ylab='Relative frequency',lab.nb.digits = 2, label = T, subtitle = 'B) Sind(symptoms)',
          fill='Var1',legend='none',palette = c('lemonchiffon1','gold1'),ylim=c(0,0.7))



#GENDER
gender_abs_freq<-table(viraldata$gender)
gender_rel_freq<-round(prop.table(gender_abs_freq),2)
freqtable_gender<-cbind(gender_abs_freq,gender_rel_freq)

colnames(freqtable_gender)<-c('Absolute frequency','Relative frequency')
kable_infection<-kable(freqtable_gender, caption = ' Table X')%>%kable_styling(bootstrap_options = c('responsive'),full_width = F, html_font = "Cambria")%>%
  row_spec(0, color = "black", background = "gainsboro")

#More males where included in the studies than females. 56% of the patients were man and 44% were woman. 


#plot
plot_gender<-ggbarplot(data.frame(gender_rel_freq),x='Var1',y='Freq',xlab='Gender',
          ylab='Relative frequency',lab.nb.digits = 2, label = T, subtitle = 'C) Gender',
          fill='Var1',legend='none',palette = c('darkolivegreen1','limegreen'),ylim=c(0,0.7))



#HOSP
hosp_abs_freq<-table(viraldata$hosp)
hosp_rel_freq<-prop.table(hosp_abs_freq)
freqtable_hosp<-cbind(hosp_abs_freq,hosp_rel_freq)

colnames(freqtable_hosp)<-c('Absolute frequency','Relative frequency')
kable_hosp<-kable(round(freqtable_hosp,2), caption = ' Table X')%>%kable_styling(bootstrap_options = c('responsive'),full_width = F, html_font = "Cambria")%>%
  row_spec(0, color = "black", background = "gainsboro")
kable_hosp

# 54% of patients required hospitalization whereas 65% didn't.

plot_hosp<-ggbarplot(data.frame(hosp_rel_freq),x='Var1',y='Freq',xlab=' ',
          ylab='Relative frequency',lab.nb.digits = 2,ylim=c(0,0.7),label = T, subtitle = 'D) Hospitalization',
          fill='Var1',legend='none',palette = c('pink1','indianred1'))


##PLOT OF CATEGORICAL DATA

ggarrange(plot_inf,plot_sind,plot_gender,plot_hosp,ncol=2,nrow=2)



#### NUMERICAL DATA
#STIME

stime_sum<-t(as.array(summary(viraldata$stime)))
rownames(stime_sum)<-'Stime'

sd_dev<-sd(viraldata$stime)
stime_<-cbind(stime_sum,sd_dev)
IQR<-IQR(viraldata$stime)
stime_<-cbind(stime_,IQR)

hist_stime<-gghistogram(viraldata, x='stime',
                  add='mean', fill='#FF6347', color='black',bins=15,title='Figure X',
                  add_density = T,xlab='Time with symptoms (days)')

box_stime<-ggboxplot(viraldata, y='stime',ylab='Time with symptoms(days)',fill='salmon',title=' ',
                add='jitter')+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggarrange(hist_stime,box_stime)

#AGE

age_sum<-t(as.array(summary(viraldata$age)))
rownames(age_sum)<-'Age'


summ_age_stime<-round(rbind(stime_sum,age_sum),2)
kable(summ_age_stime, caption = ' Table X')%>%kable_styling(bootstrap_options = c('responsive'),full_width = F, html_font = "Cambria")%>%
  row_spec(0, color = "black", background = "lightskyblue") 
  
hist<-gghistogram(viraldata, x='age',
            add='mean', fill='#E7B800', color='black',bins=25,title='Figure X',
            add_density = T,xlab='Age')


boxp<-ggboxplot(viraldata, y='age',ylab='Age',fill='lightgoldenrod',title=' ',
          add='jitter')+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggarrange(hist,boxp)
boxplot(viraldata$age)

#######################################################################################################
#######################################################################################################

#2. Perform hierarchical clustering of (scaled) gene expression levels and explore
# possible relationships between genes. How many gene clusters are observed?

#normalize gene expression levels:

genes<-viraldata[,8:57] #get gene expression values
genes_scaled<-scale(genes)
# To ensure the function has scaled the values properly, the mean should be 0 and the standard deviation shold be 1.
# we can check this using the following:
apply(genes_scaled,2,mean)
apply(genes_scaled,2,sd)

#Compute the distance between genes. 
genes_in_rows<-t(genes_scaled) #We need the genes in rows, so we transpose the matrix.
dist_genes<-dist(genes_in_rows,method='euclidian')
distmat_genes<-as.matrix(dist_genes)#contains the distance for each pair of genes.
dim(distmat_genes)

#cluster method
hcgenes<-hclust(dist_genes,method = 'average')

#Simple visualization
plot(hcgenes)
abline(h = 17, col = 'blue') #definir un treshold
#2 clusters
rect.hclust(hcgenes , k = 2, border = 2:6)
cut_avg <- cutree(hcgenes, k = 4)
cut_avg
cl1<-rownames(data.frame(cut_avg[which(cut_avg==1)]))

ggcorr(viraldata[cl1], nbreaks=8, palette='RdGy', label=TRUE, label_size=5, label_color='white')
# It is difficult to determine the exact real number of clusters, and the biological implications of each gene should 
#be studied to see in which pathway they act, thus, to get closer to the real number of clusters. 


#Better visualization:
#install.packages(c("factoextra", "dendextend"))

library(factoextra)
fviz_dend(hcgenes)
fviz_dend(hcgenes, cex = 0.75, 
          main = "Figure 4.Hierarchical clustering Dendrogram ",
          xlab='Gene Name',ylab = "Distance", sub = "",
          k=2, k_colors = c('#00AFBB','#FF6A48'),
          color_labels_by_k = TRUE,lwd=0.8)+
  geom_hline(yintercept = 17, linetype = "dashed")


fviz_dend(hcgenes, cex = 0.75, 
          main = "Figure 5.Hierarchical clustering Dendrogram ",
          xlab='Gene Name',ylab = "Distance", sub = "",
          k=4,
          color_labels_by_k = TRUE,lwd=0.8)+
  geom_hline(yintercept = 16, linetype = "dashed")


# A treshold of 17 is defined in order to create the number of clusters. In this case, and with this specific treshold, 2 clusters can be clearly observated in Figure 4.

#In cluster 1 (colored by orange), for example genes RUNDC1 and SLC2A3 seem to be associated, as the distance between them is shorter compared to other genes. The same for genes GPR180 and
# Contig63649_RC in cluster 2 (colored by blue). The distance between this two genes is short, and thus may indicate this genes are similar or involved in a same biological pathway. 
#On the contrary, for instance, gene DLT of cluster 1 and gene MELK of cluster 2, are far from each other, hence there is not a relationship between them. 


##########################################################################################################
###########################################################################################################

#3. Perform hierarchical clustering of individuals according to their (scaled) gene 
# expression levels and explore possible relationships between them. How many
# clusters of individuals are observed? Check visually whether the clustering is 
# related to infection, gender, hospitalization or ancestry

#distance between individuals
dist_individuals<-dist(genes_scaled,method = 'euclidian')
distmat_indiv<-as.matrix(dist_individuals)
dim(distmat_indiv)

#cluster method and plot
hcindividuals<-hclust(dist_individuals,method = 'average')
plot(hcindividuals)

fviz_dend(hcindividuals)

fviz_dend(hcindividuals, cex = 0.75, 
          main = "Figure 4.Hierarchical clustering Dendrogram ",
          xlab='Individual number',ylab = "Distance", sub = "",
          k=9,
          color_labels_by_k = TRUE,lwd=0.8)+
  geom_hline(yintercept = 11, linetype = "dashed")



#compare the classification with some caracteristics of the individuals:
#infection
#change the lables to make it easy to visuallize. 

xdata<-factor(viraldata$infection,labels=c('*','V'))
hcindividuals$labels=xdata #change the labels of the hclust object 

plot(hcindividuals,labels=xdata)

fviz_dend(hcindividuals, cex = 0.75, 
          main = "Figure 7. Cluster Dendrogram associated to infection ",
          xlab='Individual number',ylab = "Distance", sub = "",lwd=0.8, k=20, color_labels_by_k = TRUE)

#The gene expression levels do not provide a classification according to the type of infection. 



#gene expression levels not provide classification according to infection category. 
#So the classification of gene expression is not related with the classsfication according to infection. 

pch=viraldata$infection
pch
#gender
xdata<-factor(viraldata$gender,labels=c('*','Male'))
plot(hcindividuals,labels=xdata)

#ancestry
type_ancestry<-factor(viraldata$ancestry,labels=c('A','B','c'))
hcindividuals$labels=type_ancestry #change the labels of the hclust object 
fviz_dend(hcindividuals, cex = 0.75, 
          main = "Figure 9. Cluster Dendrogram associated to ancestry ",
          xlab='Individual number',ylab = "Distance", sub = "",lwd=0.8, k=9, color_labels_by_k = TRUE)

#hospitalization
xdata<-factor(viraldata$hosp,labels=c('NO','*'))
plot(hcindividuals,labels=xdata)

#########################################################################################################
#########################################################################################################

# 4. Perform K-means clustering with k=2 and test whether the clustering is associated 
# to (a) the kind of infection and (b) the risk of hospitalization. Interpret the results.

#we want the clustering of individuals
k<-2 
kmeans<-kmeans(genes_scaled,k)#dataset(individuals in the rows and genes scaled) and number of clusters
kmeans

kmeans$cluster


plot(genes_scaled, col=color[kmeans$cluster]) +  points(kmeans$centers,col = color, pch = 8, cex=2)
#pch=viraldata$infection

color<-c(c("#DC143C", "darkgreen"))
#posar llegenda: color _-->cluster 1   

#although the centers are quite separated (which means the mean gene expression of the two groups of clusters are a bit separated), there is not enough good separation, because some individuals overlap 
#there is a not god separation, because the classified individuals overlap each other in some points
#So, the 2 variables used (AYTL2 and TGFB3) are not able to separate the 2 groups of clusters correctly, and the 

data2<-data.frame(genes_scaled)
data2$cluster<-kmeans$cluster
data2$cluster<-as.factor(data2$cluster)

ggscatter(data2, x='AYTL2',y='TGFB3',color='cluster', title='Figure 11. k-means cluster plot ', mean.point = TRUE,mean.point.size = 7, shape='cluster')


#a)cluster-infection
cluster_infection<-prop.table(table(viraldata$infection,kmeans$cluster),2)
cluster_infection
c_i<-data.frame(cluster_infection)
c_i$Frequency<-round(c_i$Frequency,3)
colnames(c_i)<-c('Infection','Cluster','Frequency')
ggbarplot(c_i,x='Cluster',y='Frequency',fill='Infection',label=T,legend='right',title='Figure X',palette = get_palette(palette='Dark2',2))

#H0:cluster and infection are independent
#H1:cluster and infection are associated
chisq.test(cluster_infection)
#p.value<0.05-->Reject H0 and accept H1.---Semblaria que hi ha relació entre el cluster i infecció
#cluster 2--mes risc de bacterial infection
#cluster 1--mes risc de viral infection


#b)cluster-hospitalization
cluster_hosp<-table(viraldata$hosp,kmeans$cluster)
cluster_hosp
c_h<-data.frame(cluster_hosp)
colnames(c_h)<-c('Hospitalization','Cluster','Frequency')

#H0:cluster and Hosp are independent
#H1:cluster and Hosp are associated
fisher.test(cluster_hosp)
#p.value>0.05---so accept H0. cluster and Hospitalization are independent.
ggbarplot(c_h,x='Cluster',y='Frequency',fill='Hospitalization',label=T,legend='right',title='Figure X',palette = get_palette(palette='Dark2',3))


###########################################################################################################
#############################################################################################################

# 5. Perform PCA for exploring possible relationships between individuals according to 
# their (scaled) gene expression levels. Provide the variance explained plot. How 
# much variability is explained by the first two principal components? Which is the 
# eigen-value of PC1 and how can be interpreted? Check, using concentration 
# ellipses, whether PCA projections of individuals are associated to infection, 
# gender, hospitalization or ancestry. Which are the 10 genes that most contribute 
# to PC1 and PC2? (follow similar steps as in section 1.5.8 in "Solutions Exercises 
# section 2"). Discuss the results.

#PCA
pcdata<-prcomp(genes_scaled)
summary(pcdata) #explain variation of each variable. For the first PCA is 0.252

#As we have 50 variables (genes), we will obtain 50 principal components, and the first will be the one that explain the most variability of our data. 
#In this summary

#So, with the first PC we explain 25.2% of variablity. 
#By the first 2 components 33,9% of variability is explained.


fviz_cluster(kmeans, data = genes_scaled)

###############################################################################################################
###########################################################################################################

# 6. Perform a nice heatmap with dendrograms for genes and individuals, individuals 
# divided in two groups according to k-means (k=2), and annotations for infection 
# and hospitalization (similar to the one proposed in section 1.4 in "Solutions 
# Exercises section 2").

library(ComplexHeatmap)
library(circlize)
library(dendextend)

row_dend = as.dendrogram(hclust(dist(genes_scaled)))

row_dend = color_branches(row_dend, k = 2, col = c('blue','red'))


mycols <- colorRamp2(breaks = c(-4, 0, 4),
                     colors = c("limegreen", "white", "red"))

Heatmap(genes_scaled,
        name = 'gene expression',
        column_title = 'genes', row_title = 'individuals',
        cluster_rows = row_dend,
        row_split=2,
        column_names_gp = gpar(fontsize = 10),
        col=mycols)+Heatmap(viraldata$infection, name = "infection", width = unit(5, "mm"), col=c('orange','darkred'))+ 
  Heatmap(viraldata$hosp, name = "hospitalization", width = unit(5, "mm"), col=c('turquoise1','darkblue'))

# The cluster colored by red seems to be quite associated to viral infection, and it is characterized by 
# a high gene expression levels of the gens on the left and in the middle, and a low expression levels 
# of the genes in the right. On the other hand, it cannot be concluded that the cluster colored by blue is 
# associated to one kind of infection (neither bacterial nor viral infection). This cluster is characterized 
#by a low expression levels of the genes on the left and middle, and a high expression of the genes in the right.
#Regarding the hospitalization status, it seems that there is no relation with the clusters. Both of them have
# approximately a similar amount of individuals that have been hospitalizated and the ones that haven't required.
#So hospitalization can't be associated to the clusters. 


########################################################################################
##################################################################################

# 7. Test if the mean expression levels of the first gene are different between viral and 
# bacterial infections

###If the p-value is lower the significance level (usually 0.05) then we can say that we have statistically significant evidences to reject the null hypothesis, and thus to accept that the data are different in your case.
library(dplyr)
library(kableExtra)

#First of all data has to be prepared. The column which contains the infection information and the first gene are selected. In the following table, the gene expression means between the 2 infection groups can be observed:
firstgene<-viraldata %>% select(1, 8)
mean_inf<-tapply(firstgene$AYTL2,firstgene$infection,mean)
data<-t(data.frame(mean_inf))
rownames(data)<-'Mean'

kable(data, caption = 'Gene expression mean')%>%kable_classic(full_width=F)%>%row_spec(0, background = 'gainsboro')

#First of all--> shapiro test.
#H0: Data follow a normal distribution
#H1: Data do not follow a normal distribution.

shapiro.test(firstgene$AYTL2)
#p-value < 0.05. We cannot accept the H0 hypotesis, so that means that data do not follow a normal distribution. 


#As data do not follow a normal distribution-->Wilcoxon test
#H0: mean expression of viral infected patients = mean expression of bacterial infected patients
#H1: mean expression of viral infected patients is different from the mean expression of bacterial infected patients. 

wilcox.test(firstgene$AYTL2~firstgene$infection)

ggboxplot(firstgene, x='infection',y='AYTL2',
          fill='infection', palette = get_palette(palette = 'Dark2',2))+stat_compare_means(label.y=0.5,label.x = 1.3)
#p-value < 0.05, so we cannot accept the H0 hypotesis (it is rejected). Per tant, H1 is accepted, and we can conclude that
# there are diferencies significatives entre la mitjana de gene expression levels de viral infected patients and bacterial infected. 

#per tant, les mitjanes son diferents. 

############################################################################################################################
#######################################################################################################################

# 8. Test if the mean expression levels of the first gene are different among ancestry 
# groups

firstgene_ancestry<-viraldata[,7:8]
data16<-data.frame(round(tapply(firstgene_ancestry$AYTL2,firstgene_ancestry$ancestry,mean),4))
colnames(data16)<-'Mean'

kable(t(data16), caption = 'Gene expression mean')%>%kable_classic(full_width=F)%>%row_spec(0, background = 'gainsboro')
#First of all--> shapiro test.
#H0: Data follow a normal distribution
#H1: Data do not follow a normal distribution.
shapiro.test(firstgene_ancestry$AYTL2)
#mateix que l'anterior, per tant, es pot posar directament que les dades no segueixen una distribució normal.

###MORE THAN 2 MEANS-->Kruskal-Wallis test. 
#H0: all means are equal. Mean1=mean2=mean3
#H1: at lesat one of the means is different. 
kruskal.test(firstgene_ancestry$AYTL2~firstgene_ancestry$ancestry)
#p-value > 0.05. We don't have enough statistically evidences to reject the H0, so we accept that there is no stadistically significant
# diferences between the mean of the 3 ancestry groups. 
my_comparisons <- list( c("A", "B"), c("B", "C"), c("A", "C") )
ggboxplot(firstgene_ancestry, x='ancestry',y='AYTL2',
          fill='ancestry', palette = c('#7FFFD4','#ADFF2F','#32CD32'))+stat_compare_means(comparisons = my_comparisons)+ stat_compare_means(label.y = 0.9)

###############################################################################################################################
##########################################################################################################################

# 9. Test whether mean expression levels of the first and second genes are equal for 
# viral infections.

#Prepare the data
first_second_gene<-viraldata%>% filter(., infection == 'viral infection')%>%select(8,9)
data_9<-data.frame(round(apply(first_second_gene,2,mean),4)) #see the means
colnames(data_9)<-'Gene expression mean'

kable(t(data_9), caption = 'Table X')%>%kable_classic(full_width=F)%>%row_spec(0, background = 'gainsboro')
#As we have selected only those viral infected patients, we have to check again if the data is normally distributed or not

#H0: data is normally distributed
#H1: data is not normally distributed

shapiro.test(first_second_gene$AYTL2)
shapiro.test(first_second_gene$TGFB3)
#both p-values are >0.05. So data follow a normal distribution. 

#T-TEST.
#H0: mean expression of gene AYTL2 is equal to the mean expression of gene TGFB3 
#H1: means are different.


d<-first_second_gene$AYTL2-first_second_gene$TGFB3
t.test(d,mu=0)

#p-value < 0.05. Rebutjem H0, i acceptem H1, per tant les mitjanes son differents entre els pacients virics del gen 1 i del gen 2
library(tidyr)
a<-pivot_longer(first_second_gene, AYTL2:TGFB3)
ggboxplot(a, x='name',y='value',legend.title='Gene', legend='right', font.legend=12, xlab='Gene', ylab='Expression',
          fill='name',palette=c('#FF7F50','#D53E4F') ) + stat_compare_means(method='t.test',paired = T,label.x = 1.3)


##############################################################################################
######################################################################################

# 10. Perform a nonparametric test for association of the kind of infection (viral or bacterial) and the risk of hospitalization. 
#Provide the OR of the risk of hospitalization for viral vs bacterial infections.

#Table of frequencies
data_10<-table(viraldata$infection,viraldata$hosp)

kable(data_10, caption = 'Table X')%>% kable_classic(full_width = F)%>%row_spec(0, background = 'gainsboro')
# H0: Infection and hospitalization are independent variables
# H1: Infection and hospitalization are related variables. 
fisher.test(data_10)


#p-value > 0.05. ACcEPTEM H0. --SOn independents---ods ratio equal to 1.

#The odds ratio turns out to be 0.6699788. 
#The odds of being hospitalizated for a bacetrial infected patient is 0.6699 times the odds that a viral patient is hospitalizated. 
# in other word, we can say that the odds that a bacterial infected patient is hospializated are lowered by about 33%


#The fisher test also provides us the 95% confidnece interval for the odds ratio. In this case is....


# Since the confidence interval contains the odds ratio value of 1, then the calculated odds ratio is not considered statistically significant,
# because it means that the expected true odds ratio can be below or above 1, so we don't really know if being ifectd with bacterial or viral, 
#incrases or decreases the odds og being hospitalizated. 


mosaicplot(data_10,color = T)
oddsratio(data_10)

# To visually represent the data, it can be drawn a <span style="color:blue"> mosaic plot</span>:

library(ggmosaic)
ggplot(data = viraldata) +
  geom_mosaic(aes(x = product(infection), fill=hosp)) +
  theme_bw()+xlab('Infection')+ylab('Hospitalization')+theme(legend.position = 'none')+ggtitle('Figure X')
###########################################################################################
##########################################################################################

# 11. Test the normality of expression levels of the 50 genes (use function apply). How 
# many genes are not normally distributed and which are their names?

genes<-viraldata[,8:57]
test<-apply(genes,2,function(x) shapiro.test(x)$p.value)
sum(test<0.05)

a<-rbind(genes,test)#add the last row--p-value

p_valeminor<-data.frame(test[which(test<0.05)])
normal_genes<-data.frame(rownames(p_valeminor))
colnames(normal_genes)<-'Gene name'
kable(normal_genes, caption = 'Table X')%>%kable_styling(bootstrap_options = c("striped", "hover", "condensed",'bordered'),full_width = F)%>%row_spec(0,2, background = 'lightskyblue')


###############################################################################################
###############################################################################################

# #12. Identify those genes that are differentially expressed between viral and bacterial 
# infections (use function apply). Create a function that checks whether the gene 
# expression levels are normally distributed or not and, accordingly, applies the 
# most appropriate test for comparing gene expression levels between viral and 
# bacterial infections. Adjust the p-values for multiple testing according to an fdr 
# threshold equal to 0.1. Interpret the results


#SHapiro
#distribució no normal--->wilcoxon
#wilcoxon not equal <0.05---dif
#wilcoxon>0.05---no dif

#distribució normal-->ttest
#ttest equal--->no son dif
##ttest not equal-->dif


deg<-function(x){
  if (shapiro.test(x)$p.value<0.05){
    wilcox.test(x~viraldata$infection)$p.value
  }
  else{
    if (var.test(x~viraldata$infection)$p.value<0.05){
      t.test(x~viraldata$infection,var.equal=F)$p.value
    }
    else{
      t.test(x~viraldata$infection,var.equal=T)$p.value
    }
  }
}


pvalues<-apply(genes, 2, deg)
adjusted_pvalues<-p.adjust(pvalues, method = "fdr", n = length(pvalues))
sum(adjusted_pvalues<0.01)
data_pvaladj<-data.frame(adjusted_pvalues[which(adjusted_pvalues<0.01)])
rownames(data_pvaladj)


#############################################################################################################################
############################################################################################################################

#13. Consider a regression model for the kind of infection as a function of gender, age 
# and ancestry and the first 10 genes (scaled). Use stepwise variable selection and 
# denote the selected model as “best.model”. Interpret the obtained model.  

#infection--dichotomous variable-->logistic regression


model1<-glm(viraldata$infection~viraldata$gender+viraldata$age+viraldata$ancestry+genes_scaled[,1:10], family = binomial())
summary(model1)
coef(model1)

#Moltes variables--volem les que estan relacionades nomes. 

backward<-step(model1,direction = 'backward')
summary(backward)
backward$coefficients
formula(backward)

best.model<-glm(viraldata$infection ~ viraldata$gender + viraldata$ancestry + 
                  genes_scaled[, 1:10],family=binomial())


###############################################################################################################################
#############################################################################################################################

# 14. Analyze the classification ability of “best.model” (ROC curve and AUC) 
# according to the following schemes:
#   a. Apparent validation of “best.model” using the same data that was used 
# for model building.
# b. Cross-validation with k = 5 for “best.model”.
# c. Though the cv-classification is better than the apparent classification, it still
# is over-estimating the real classification of "best-model". Discuss why and 
# how to obtain a more accurate classification estimation (slides 262:264).

lp<-best.model$linear.predictors


library(ROCR)

pred<-prediction(lp,viraldata$infection)
perf<-performance(pred,'tpr','fpr')
plot(perf)
abline(a=0,b=1)

auc<-slot(performance(pred,"auc"), "y.values")[[1]]
auc
######################

#BiocManager::install("ROCR")
library(ROCR)

k<-5
n<-nrow(viraldata)
fold<-sample(as.numeric(cut((1:n),breaks = k)))

fold

pred <- NULL # vector of predictions

for(i in 1:k){
  indTest <- which(fold==i)   # Test indices 
  indTrain <- which(fold!=i)  # Train indices
  model.i <- glm(viraldata$infection~viraldata$gender + viraldata$ancestry + genes_scaled[, 1:10], data=viraldata[indTrain,], family=binomial())  # Adjust the model with training data
  pred.i <- predict(model.i, newdata=viraldata[indTest, ])   # Predicts test data at step i
  pred[indTest] <- pred.i   # Store predicted values for test data at step i 
}  

prediction <- prediction(pred, viraldata$infection) 
perf <- performance(prediction, "tpr", "fpr" )
plot(perf)
abline(a=0, b= 1)

auc<-slot(performance(prediction,"auc"), "y.values")[[1]]
auc
#################################################################################################################################
####################################################################################################################################
# 15. Consider a regression model for the kind of infection as a function of all 50 genes 
# (scaled) and adjusted by age. Perform variable selection with LASSO and 
# interpret the results.

library(glmnet)
#all genes scaled + age
age<-viraldata[,6]
genes_age<-cbind(age,genes_scaled)
x<-genes_age

#
y<-viraldata[,1]

mlasso <- glmnet(x, y, standardize=TRUE, alpha=1,family='binomial')

plot(mlasso)
#######################################################################################################################################
####################################################################################################################################
# 16. Obtain Kaplan-Meier survival curves for the time of symptoms as a function of the 
# kind of infection and test for the significance of the difference in duration of 
# symptoms. Discuss the results

#install.packages("KMsurv")
library(survival)
library(KMsurv)

kmbysex<-survfit(Surv(viraldata$stime,viraldata$sind)~ viraldata$infection)
summary(kmbysex)
plot(kmbysex)
plot(kmbysex, main='FIgure x',xlab='time',ylab='survival',col=2:3)
legend('topright',col=2:3,legend=c('bacterial','viral'),lty=1)

#Individuals with a viral infection take less time with the symptoms.
#At 15 days, only aprox 40% of the viral infected patients remain with symptoms,
#and aprox 60% of the bacterial patients remain with symptoms

#Are these diferences significative?

#H0: no diference in the survival for the 2 groups
#H1: the 2 groups have different survival distributions

survdiff(Surv(viraldata$stime,viraldata$sind)~ viraldata$infection)


###############################################################################################################################
################################################################################################################################

#17. Perform a Cox regression model for duration symptoms as a function of the 
#covariates (ignore gene expression levels). Discuss the results

coxmodel<-coxph(Surv(viraldata$stime,viraldata$sind)~viraldata$infection+viraldata$gender+viraldata$hosp+viraldata$age+viraldata$ancestry)
coxmodel
