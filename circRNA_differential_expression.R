### R project

#setwd("E:/circRNA data")
setwd("C:/Users/owner/Downloads")

##loading packages----
library("ggplot2")
library("edgeR")
library("gridExtra")
library("rafalib")

options(stringsAsFactors = F)


## data filtering and cleaning -----

#read files
liner_genes=readRDS('htseq_counted_dedup.RDS')
load('good_samps_no_noncontrol.rdata')
ens_table=read.csv('ens_genes_with_name',sep = '\t')
load('mata_data_info.rdata')

#remove samples not in matadata
keep_samples=match(mata_data_orderd$name,colnames(liner_genes))
liner_genes=liner_genes[,keep_samples]

#add gene name liner genes
rownames(liner_genes)<- sapply(strsplit(rownames(liner_genes),'.',fixed = T),function(x) x[1])
liner_genes<- liner_genes[rownames(liner_genes) %in% ens_table$hg19.ensGene.name2 ,]
m=match(rownames(liner_genes),ens_table$hg19.ensGene.name2)
gene_symb=ens_table$hg19.ensemblToGeneName.value[m]
liner_genes=cbind.data.frame(gene_symb,liner_genes)

#add gene name to circ and liner junction
gene_names=sapply(strsplit(circ$name,':',fixed = T),function(x) x[2])
m=match(gene_names,ens_table$X.hg19.ensGene.name)
gene_symb=ens_table$hg19.ensemblToGeneName.value[m]
circ=cbind.data.frame(gene_symb,circ)


circ_e=circ[,-c(1:5)]
liner_genes_e<- liner_genes[,-1]
rownames(circ_e)=circ$name

#match liner and circ
m=match(colnames(circ_e),colnames(liner_genes_e))
liner_genes_e=liner_genes_e[,m]

## create circ to linear ratio----
tot_circ<- as.vector(apply(data.matrix(circ_e), 2, sum))
tot_linear<- apply(data.matrix(liner_genes_e), 2, sum)

circ_lin_ratio<- tot_circ/tot_linear


## boxplot circ to linear and t-test -----
mata_data_orderd$ASD.CTL<- droplevels(mata_data_orderd$ASD.CTL, "NCTL") 
mypar(1,2)

#cortex only
boxplot(log(circ_lin_ratio)[mata_data_orderd$RegionID!="vermis"]~ mata_data_orderd$ASD.CTL[mata_data_orderd$RegionID!="vermis"],
        ylab="circRNA to linear RNA ratio", xlab= "", main = "Cortex")

t.test(log(circ_lin_ratio)[mata_data_orderd$RegionID!="vermis"]~ mata_data_orderd$ASD.CTL[mata_data_orderd$RegionID!="vermis"] )

#varmis only
boxplot((log(circ_lin_ratio)[mata_data_orderd$RegionID=="vermis"])~ mata_data_orderd$ASD.CTL[mata_data_orderd$RegionID=="vermis"],
        ylab="circRNA to linear RNA ratio", xlab = "", main ="Vermis")

t.test(log(circ_lin_ratio)[mata_data_orderd$RegionID=="vermis"]~ mata_data_orderd$ASD.CTL[mata_data_orderd$RegionID=="vermis"] )





## detiled diagnosis (chromosome15 duplication)-----

brain_region<- as.vector(mata_data_orderd$RegionID)
brain_region[brain_region != "vermis"] <- "cortex"

ggplot(data.frame(log(circ_lin_ratio),ncol=1), aes(y=log(circ_lin_ratio), x=brain_region, fill=mata_data_rel$Detailed.Diagnosis))+
  geom_boxplot()+ ylab("circular to linear reads (log)")+ xlab("")+ ggtitle("Chromosome 15q duplication")+
  scale_fill_discrete(labels = c("Control", "ASD without Chr15q duplication", "Chr15q duplication"))+labs(fill = "Detailed diagnosis") 

#anova cortex
circ_lin_ratio_cortex<- circ_lin_ratio[mata_data_orderd$RegionID != "vermis"]
mata_data_orderd_cortex<- mata_data_orderd[mata_data_orderd$RegionID!= "vermis",]

summary(aov(circ_lin_ratio_cortex~mata_data_orderd_cortex$Age+mata_data_orderd_cortex$RIN+mata_data_orderd_cortex$SeqBatch+mata_data_orderd_cortex$Detailed.Diagnosis))

#anova vermis
circ_lin_ratio_vermis<- circ_lin_ratio[mata_data_orderd$RegionID == "vermis"]
mata_data_orderd_vermis<- mata_data_orderd[mata_data_orderd$RegionID== "vermis",]

summary(aov(circ_lin_ratio_vermis~mata_data_orderd_vermis$Age+mata_data_orderd_vermis$RIN+mata_data_orderd_vermis$SeqBatch+mata_data_orderd_vermis$Detailed.Diagnosis))


## Seizures------
m<- !is.na(mata_data_missing$Seizures)
Seizures_info<- mata_data_missing$Seizures[m]
circ_lin_ratio_seiz<- circ_lin_ratio[m]
brain_region_seiz<- brain_region[m]

ggplot(data.frame(log(circ_lin_ratio_seiz),ncol=1), aes(y=log(circ_lin_ratio_seiz), x=brain_region_seiz, fill=Seizures_info))+
  geom_boxplot()+ ylab("circular to linear reads (log)")+ xlab("")+ ggtitle("Seizures")+
  scale_fill_discrete(labels = c("Control", "ASD with Seizures", "ASD with Seizures"))+labs(fill = "Seizures") 


#anova cortex
summary(aov(circ_lin_ratio_cortex~mata_data_orderd_cortex$Age+mata_data_orderd_cortex$RIN+mata_data_orderd_cortex$SeqBatch+
              mata_data_missing$Seizures[mata_data_orderd$RegionID!= "vermis"]))

#anova vermis
summary(aov(circ_lin_ratio_vermis~mata_data_orderd_vermis$Age+mata_data_orderd_vermis$RIN+mata_data_orderd_vermis$SeqBatch+
              mata_data_missing$Seizures[mata_data_orderd$RegionID== "vermis"]))


## Differential expression ----

#filter only cortex samples
circ_cortex=circ_e[,mata_data_rel$RegionID=='cortex']
liner_genes_cortex=liner_genes_e[,mata_data_rel$RegionID=='cortex']
mata_data_rel_cortex=mata_data_rel[mata_data_rel$RegionID=='cortex',]
mata_data_orderd_cortex<-mata_data_orderd[mata_data_orderd$RegionID!='vermis',] 
mata_data_missing_cortex=mata_data_missing[mata_data_rel$RegionID=='cortex',]
mata_data_rel_cortex<- cbind.data.frame(mata_data_rel_cortex, mata_data_missing_cortex$Seizures)
mata_data_rel_cortex$ASD.CTL<- relevel(mata_data_rel_cortex$ASD.CTL, ref = "CTL")

de_liner=DGEList(counts=liner_genes_cortex,group = mata_data_rel_cortex$ASD.CTL)
de_liner <- calcNormFactors(de_liner)

##dif exp of linear reads
de_liner<- DGEList(counts=liner_genes_cortex,group = mata_data_rel_cortex$ASD.CTL)
mata_data_rel_cortex$ASD.CTL<- relevel(mata_data_rel_cortex$ASD.CTL, ref = "CTL")
de_liner <- calcNormFactors(de_liner)

#filter out non expressing genes
saf=55
exp_liner=apply(de_liner$counts,1,function(x) sum(x>1)>saf)
de_liner=de_liner[exp_liner,]
mata_data_rel_cortex$Detailed.Diagnosis <- relevel(mata_data_rel_cortex$Detailed.Diagnosis, ref="NO_15")
design_liner <- model.matrix(~ASD.CTL+Age+Sex+RIN+SeqBatch,data=mata_data_rel_cortex)
rownames(design_liner) <- colnames(de_liner)

de_liner <- estimateDisp(de_liner, design_liner, robust=TRUE)
plotBCV(de_liner)
fit_liner <- glmFit(de_liner, design_liner)
lrt_liner <- glmLRT(fit_liner,coef = 3)
topTags(lrt_liner)

res_liner<- data.frame(topTags(lrt_liner,n=nrow(de_liner)))

## dif exp of circ reads
de_circ<- DGEList(counts=data.matrix(circ_cortex),group = mata_data_rel_cortex$ASD.CTL,lib.size = de_liner$samples$lib.size)
de_circ <- calcNormFactors(de_circ)

de_circ$samples$norm.factors<-de_liner$samples$norm.factors

#filter out non expressing circs
a<- data.frame(de_circ$counts)
saf=50
exp_circ=apply(a,1,function(x) sum(x>1)>saf)
de_circ$counts=de_circ$counts[exp_circ,]

design <- model.matrix(~Age+Sex+RIN+SeqBatch+ASD.CTL,data=mata_data_rel_cortex)
rownames(design) <- colnames(de_circ)

de_circ <- estimateDisp(de_circ, design, robust=TRUE)

fit <- glmFit(de_circ, design)
lrt <- glmLRT(fit,coef = 2)
plotMD(lrt)

res<- data.frame(topTags(lrt,n=nrow(de_circ)))

plot(res$logFC, -log10(res$PValue), cex = 0.5)

## volcano plot of logFC vs FDR ----
library("RColorBrewer")
library("ggrepel")
cols <-brewer.pal(11,name =  'Spectral')
cols[11]='grey'

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#plot for circs
circRNA_Density <- get_density(res$logFC, -log(res$FDR))

p1 <- ggplot(res, aes(y = -log(res$PValue), x = logFC))
p2<- p1 + geom_point(aes(col = circRNA_Density))+scale_color_gradientn(guide = 'colorbar',colors=rev(cols))+ xlim(-0.06, 0.06)+ ylim(-0.01, 12)+
  geom_hline(yintercept = 6.352876,color = "Red",linetype="dashed")+labs(title="circRNAs", x="log (Fold Change)", y= "-log p-Value")+theme_bw()+theme(legend.position = c(0.18, 0.18))+
  theme(plot.title = element_text(size = 18), axis.title.x = element_text(size = 13), axis.title.y = element_text( size = 13))

p2

sum(res$logFC>0)/nrow(res)

#plot for linear reads
linear_Density <- get_density(res_liner$logFC, -log(res_liner$FDR))
l1 <- ggplot(res_liner, aes(y = -log(PValue), x = logFC))
l2<- l1 + geom_point(aes(col = linear_Density))+scale_color_gradientn(guide = 'colorbar',colors=rev(cols))+ ylim(-0.01, 26)+ xlim(-0.11,0.11)+
  geom_hline(yintercept = -log(0.007926753), color = "Red", linetype="dashed")+labs(title="Linear Genes", x="log (Fold Change)")+ theme_bw()+theme(legend.position = c(0.16, 0.19))+
  theme(plot.title = element_text(size = 18), axis.title.x = element_text(size = 13), axis.title.y = element_text( size = 13))
l2

sum(res_liner$logFC>0)/nrow(res_liner)


#creat volvano plot for linear genes that express high level of circs
high_exp_circ <- sapply(strsplit(rownames(de_circ$counts),':',fixed = T),function(x) x[2])
m<- match(high_exp_circ, ens_table$X.hg19.ensGene.name)
high_exp_gene<- ens_table$hg19.ensGene.name2[m]
m<-  match(high_exp_gene, rownames(res_liner))
m<- m[!is.na(m)]
dup<- !duplicated(m)
m<- m[dup]

res_liner_exp<- res_liner[m,]
linear_Density <- get_density(res_liner_exp$logFC[-8], -log(res_liner_exp$FDR[-8]))
e1 <- ggplot(res_liner_exp[-8,], aes(y = -log(PValue), x = logFC))
e2<- e1 + geom_point(aes(col = linear_Density))+scale_color_gradientn(guide = 'colorbar',colors=rev(cols))+ 
  ylim(-0.01, 17)+ xlim(-0.06,0.06)+ geom_hline(yintercept = -log(0.007926753), color = "Red", linetype="dashed")+
  labs(title="Genes that Express circRNAs", x="log (Fold Change)")+ theme_bw()+theme(legend.position = c(0.15, 0.18))+
  theme(plot.title = element_text(size = 18), axis.title.x = element_text(size = 13), axis.title.y = element_text( size = 13))
e2

sum(res_liner_exp$logFC>0)/nrow(res_liner_exp)

## circRNA validation ----
#from article: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0571-3#Sec2

load("val circs.csv")
val_circs<- val_circs[1:24,1:6]
colnames(val_circs)<- c("gene", "chr", "start", "end", "region in genome", "junction reads")

all_matching<- NULL

for ( k in 1:nrow(val_circs)){
  matching_circs<- NULL
  for ( i in c(as.numeric(val_circs$start[k]-1):as.numeric(val_circs$start[k]+1))){
    matching_circs <- rbind.data.frame(matching_circs, circ[i == circ$s,])
  }
  for ( j in c(as.numeric(val_circs$end[k]-1):as.numeric(val_circs$end[k]+1))){
    all_matching<- rbind.data.frame(all_matching, matching_circs[matching_circs$e== j,])
  }
}


circ_index<- c()
for (i in 1:nrow(all_matching)){
  circ_name<- paste(all_matching$name[i])
  circ_index<- c(circ_index, which(rownames(res)== circ_name))
}

valid_circs<- NULL
for (i in 1:nrow(res)){
  if (i %in% circ_index){valid_circs[i]<-T}
  else {valid_circs[i]<-F}
}

a1<- ggplot(data = res,aes(x=logFC))+geom_density(fill="blue",alpha=0.4)+ xlim(-0.03,0.03)+labs(title="logFC all circRNAs")+ geom_vline(xintercept = mean(res$logFC))
a2<- ggplot(data = res[valid_circs,],aes(x=logFC))+geom_density(fill="red",alpha=0.4) + xlim(-0.03,0.03)+labs(title="logFC validated circRNAs")+geom_vline(xintercept = mean(res$logFC[valid_circs]))
logFC_density<- grid.arrange(a1, a2, ncol=1)


#permutation test
mean_logFC_per=c()
for(i in 1:1000) {
  shuff=sample(res$logFC)
  mean_logFC_per[i]=(mean(shuff[res$valid_circs==T]))
}
density(mean_logFC_per, main = "Permutation test for validated circs", xlab = "Permutated mean")
abline(v=mean(res$logFC[valid_circs]),col='red')

p_val_per<- (sum(mean_logFC_per<=mean(res$logFC[valid_circs]))+1)/(1000+1)
p_val_per

## correlation of RNA editing levels and circ to lin ratio ----

#read RNA editing file
RNA_edit_ba9 <- read.csv("RNA_edit_ba9.csv")
brain_data_edit <- read.csv("brain_data_edit.csv")

#choose only samples that have RNA editind data and match all tables. 
table(paste(mata_data_orderd$BrainID, mata_data_orderd$RegionID) %in%
        paste(brain_data_edit$BrainID, brain_data_edit$RegionID))
sample_name<- paste("Sample","_",mata_data_orderd$BrainID,"_",mata_data_orderd$RegionID, sep = "")
mata_data_orderd<- cbind.data.frame(sample_name, mata_data_orderd)

m<- brain_data_edit$sample_name %in% mata_data_orderd$sample_name 
brain_data_edit<- brain_data_edit[m,]

m<- mata_data_orderd$sample_name %in% brain_data_edit$sample_name
mata_data_orderd<- mata_data_orderd[m,]

mata_data_orderd<- mata_data_orderd[!duplicated(mata_data_orderd$sample_name),]

m<- match(brain_data_edit$sample_name, mata_data_orderd$sample_name)
mata_data_orderd<- mata_data_orderd[m,]

circ_e<- circ_e[,which((colnames(circ_e) %in% mata_data_orderd$name))]
m<- match(mata_data_orderd$name, colnames(circ_e))
circ_e<- circ_e[,m]

liner_genes_e<- liner_genes_e[,which((colnames(liner_genes_e) %in% mata_data_orderd$name))]
m<- match(mata_data_orderd$name, colnames(liner_genes_e))
liner_genes_e<- liner_genes_e[,m]

#only ba9 data
RNA_ba9_cut<- RNA_edit_ba9[,-c(1:26)]
RNA_ba9_cut<- RNA_ba9_cut[,which((colnames(RNA_ba9_cut) %in% brain_data_edit$sample_name))]

#convert empty cells to NAs
convert_to_NA<- function(x) {
  for (i in 1:length(x)){
    if (x[i]== "") {x[i]<- NA}
  }
  return(x)
}

RNA_ba9_cut<- apply(RNA_ba9_cut, 1, convert_to_NA)
RNA_ba9_cut<- t(RNA_ba9_cut)

#count how many rows contain NA 

keep_row<- function(x){
  keep_list<- rep(T, times = nrow(x))
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)){
      if (is.na(x[i,j])== T) {keep_list[i]<- F
      break} 
    }
  }
  return(keep_list)
}
table(keep_row(RNA_ba9_cut))
# FALSE  TRUE 
# 2089  1225 
#most genes have at least 1 NA value

#create table of numbers
for (i in 1:nrow(RNA_ba9_cut)){
  for (j in 1:ncol(RNA_ba9_cut)){
    if (is.na(RNA_ba9_cut[i,j])==T){RNA_ba9_cut[i,j]<- NA} else {RNA_ba9_cut[i,j]<-as.numeric(unlist(strsplit(RNA_ba9_cut[i,j], "^",fixed = T))[1])}
  }
}

#create mean without NAs
mean_editing <-apply(RNA_ba9_cut, 2, function(x) mean(as.numeric(x[!is.na(x)])))

#match samples in mean editing and circ to linear ratio
names(circ_lin_ratio)<- mata_data_orderd$sample_name
shared_samples<- match(names(mean_editing), names(circ_lin_ratio))
circ_to_linear_ba9<- circ_lin_ratio[shared_samples]

cor_data<- cbind.data.frame(circ_to_linear_ba9,mean_editing)
ggplot(data=cor_data, aes(circ_to_linear_ba9,mean_editing))+ geom_point()+ 
  ggtitle("circular to linear ratio vs. RNA editing in ba9")+ xlab("circular to linear RNA ratio")+
  ylab(" mean of RNA editing rates")+ 
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text( size=14)
  )

cor(circ_to_linear_ba9, mean_editing)
#-0.200583

cor.test(circ_to_linear_ba9, mean_editing)
# Pearson's product-moment correlation
# 
# data:  circ_to_linear_ba9 and mean_editing
# t = -1.5322, df = 56, p-value = 0.1311
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.43627546  0.06086623
# sample estimates:
#       cor 
# -0.200583 


#only ba41-42-22 data
RNA_edit_ba42 <- read.csv("RNA_edit_ba42.csv")

convert_to_NA<- function(x) {
  for (i in 1:length(x)){
    if (x[i]== "") {x[i]<- NA}
  }
  return(x)
}
RNA_edit_ba42<- apply(RNA_edit_ba42, 1, convert_to_NA)
RNA_edit_ba42<- t(RNA_edit_ba42)

keep_row(RNA_edit_ba42)
# FALSE  TRUE 
# 2026   386 
# only 386 row dont contain NAs

#create table of numbers
for (i in 1:nrow(RNA_edit_ba42)){
  for (j in 1:ncol(RNA_edit_ba42)){
    if (is.na(RNA_edit_ba42[i,j])==T){RNA_edit_ba42[i,j]<- NA} else {RNA_edit_ba42[i,j]<-as.numeric(unlist(strsplit(RNA_edit_ba42[i,j], "^",fixed = T))[1])}
  }
}


#create mean without NAs
mean_editing_ba42 <-apply(RNA_edit_ba42, 2, function(x) mean(as.numeric(x[!is.na(x)])))

mata_data_ba42<- mata_data_orderd[which(mata_data_orderd$RegionID== "ba41-42-22"),]

samples_ba42<- vector()
for (i in 1:length(mean_editing_ba42)) {samples_ba42[i]<-unlist(strsplit(names(mean_editing_ba42)[i], "_"))[2]}
m<- samples_ba42 %in% mata_data_ba42$BrainID
samples_ba42<- samples_ba42[m] 
mean_editing_ba42<- mean_editing_ba42[m]

m<- mata_data_ba42$BrainID %in% samples_ba42 
mata_data_ba42<- mata_data_ba42[m,]

#match samples in mean editing and circ to linear ratio
names(circ_lin_ratio)<- mata_data_orderd$sample_name
shared_samples<- match(mata_data_ba42$sample_name, names(circ_lin_ratio))
circ_to_linear_ba42<- circ_lin_ratio[shared_samples]

cor_data<- cbind.data.frame(circ_to_linear_ba42,mean_editing_ba42)
ggplot(data=cor_data, aes(circ_to_linear_ba42,mean_editing_ba42))+ geom_point()+ 
  ggtitle("circular to linear ratio vs. RNA editing in ba41_42_22")+ xlab("circular to linear RNA ratio")+
  ylab(" mean of RNA editing rates")+ 
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text( size=14)
  )

cor(circ_to_linear_ba42, mean_editing_ba42)
# -0.07218531

cor.test(circ_to_linear_ba42, mean_editing_ba42)
# Pearson's product-moment correlation
# 
# data:  circ_to_linear_ba42 and mean_editing_ba42
# t = -0.50142, df = 48, p-value = 0.6184
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3436286  0.2103897
# sample estimates:
#         cor 
# -0.07218531 

#only vermis data
RNA_edit_vermis <- read.csv("RNA_edit_vermis.csv")

RNA_edit_vermis<- apply(RNA_edit_vermis, 1, convert_to_NA)
RNA_edit_vermis<- t(RNA_edit_vermis)

table(keep_row(RNA_edit_vermis))
# FALSE  TRUE 
# 3558   782
# only 782 row dont contain NAs

#create table of numbers
for (i in 1:nrow(RNA_edit_vermis)){
  for (j in 1:ncol(RNA_edit_vermis)){
    if (is.na(RNA_edit_vermis[i,j])==T){RNA_edit_vermis[i,j]<- NA} 
    else {RNA_edit_vermis[i,j]<-as.numeric(unlist(strsplit(RNA_edit_vermis[i,j], "^",fixed = T))[1])}
  }
}


#create mean without NAs
mean_edit_vermis <-apply(RNA_edit_vermis, 2, function(x) mean(as.numeric(x[!is.na(x)])))

#match RNA edit table to mata data
mata_data_vermis<- mata_data_orderd[which(mata_data_orderd$RegionID== "vermis"),]
m<- names(mean_edit_vermis) %in% mata_data_vermis$sample_name

mean_edit_vermis<- mean_edit_vermis[m]

m<- match(names(mean_edit_vermis), mata_data_vermis$sample_name)
mata_data_vermis<- mata_data_vermis[m,]

#match samples in mean editing and circ to linear ratio
shared_samples<- match(mata_data_vermis$sample_name, names(circ_lin_ratio))
circ_to_linear_vermis<- circ_lin_ratio[shared_samples]

cor_data<- cbind.data.frame(circ_to_linear_vermis,mean_edit_vermis)
ggplot(data=cor_data, aes(circ_to_linear_vermis,mean_edit_vermis))+ geom_point()+ 
  ggtitle("circular to linear ratio vs. RNA editing in vermis")+ xlab("circular to linear RNA ratio")+
  ylab(" mean of RNA editing rates")+ 
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text( size=14)
  )

cor(circ_to_linear_vermis, mean_edit_vermis)
# 0.234254

cor.test(circ_to_linear_vermis, mean_edit_vermis)
# Pearson's product-moment correlation
# 
# data:  circ_to_linear_vermis and mean_edit_vermis
# t = 1.7208, df = 51, p-value = 0.09135
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.03847597  0.47450352
# sample estimates:
#      cor 
# 0.234254 


