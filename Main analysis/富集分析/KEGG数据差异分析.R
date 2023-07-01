# 清空环境 ----------------------------------------------------------------
library("pacman")
p_unload(p_loaded(), character.only = TRUE)
rm(list = ls())
gc()

# 加载需要的R包 -----------------------------------------------------------------
library("dplyr")
library("vegan")
library("clusterProfiler")
library("tidyverse")
library("ggplot2")
library("pathview")
library("org.Hs.eg.db")
library("DOSE")

# 设置工作路径并读取文件 ------------------------------------------------------------------
setwd("C:/Users/Dell/Desktop/每日记录/20230629/")
KO_pathway <- read.csv("KEGG_pathway_extracted.csv",check.names = F)
KEGG_module <- read.csv("KEGG_module_extracted.csv",check.names = F)
species <- read.csv("species.csv",check.names = F,row.names = 1)
metabolites <- read.csv("metabolites.csv",check.names = F,row.names = 1)
pathway <- read.csv("function.csv",check.names = F,row.names = 1)
RC_metadata <- read.csv("RC_metadata.csv",check.names = F,row.names = 1)
#4个特殊样本名的重命名
colnames(KO_pathway)[which(colnames(KO_pathway)=="PWV0298-1PWV0676")] <- "PWV0298"
colnames(KO_pathway)[which(colnames(KO_pathway)=="PWV0157-1PWV0673")] <- "PWV0157"
colnames(KO_pathway)[which(colnames(KO_pathway)=="PWV0255-1PWV0699")] <- "PWV0255"
colnames(KO_pathway)[which(colnames(KO_pathway)=="PWV0159-1PWV0743")] <- "PWV0159"
colnames(KEGG_module)[which(colnames(KEGG_module)=="PWV0298-1PWV0676")] <- "PWV0298"
colnames(KEGG_module)[which(colnames(KEGG_module)=="PWV0157-1PWV0673")] <- "PWV0157"
colnames(KEGG_module)[which(colnames(KEGG_module)=="PWV0255-1PWV0699")] <- "PWV0255"
colnames(KEGG_module)[which(colnames(KEGG_module)=="PWV0159-1PWV0743")] <- "PWV0159"

# 数据预处理 -------------------------------------------------------------------
#菌群与metacyc通路信息log处理
species <- log(species+1,2)
pathway <- log(pathway+1,2)
#metabolites信息Zscore处理
metabolites <- scale(metabolites,center = F)
species <- as.data.frame(t(species))
pathway <- as.data.frame(t(pathway))
metabolites <- as.data.frame(t(metabolites))
#提取目标菌群/metacyc通路/代谢物
metabolites <- metabolites[,which(colnames(metabolites)%in%c("b21","b283","b133","b333"))]
species <- species[,which(colnames(species)%in%c("Adlercreutzia_equolifaciens","Collinsella_aerofaciens","Roseburia_hominis","Fusicatenibacter_saccharivorans",
                                                 "Dorea_longicatena","Alistipes_putredinis"))]
pathway <- pathway[,which(colnames(pathway)%in%c("a205","a117"))]
#KEGG通路计算相对丰度
KEGG_pathway <- aggregate(KO_pathway[,6:ncol(KO_pathway)], list(KO_pathway$Pathway),sum)
KEGG_pathway_id <- KEGG_pathway[,1]
KEGG_pathway <- decostand(KEGG_pathway[,2:ncol(KEGG_pathway)],method="total",2)
rownames(KEGG_pathway) <- KEGG_pathway_id
KEGG_pathway <- as.data.frame(t(KEGG_pathway))
#KEGG module计算相对丰度
KEGG_module <- aggregate(KEGG_module[,4:ncol(KEGG_module)], list(KEGG_module$module),sum)
KEGG_module_id <- KEGG_module[,1]
KEGG_module <- decostand(KEGG_module[,2:ncol(KEGG_module)],method="total",2)
rownames(KEGG_module) <- KEGG_module_id
KEGG_module <- as.data.frame(t(KEGG_module))
#KEGG KO计算相对丰度
KEGG_gene <- aggregate(KO_pathway[,6:ncol(KO_pathway)], list(KO_pathway$KO),sum)
KEGG_gene_id <- KEGG_gene[,1]
KEGG_gene <- decostand(KEGG_gene[,2:ncol(KEGG_gene)],method="total",2)
rownames(KEGG_gene) <- KEGG_gene_id
KEGG_gene <- as.data.frame(t(KEGG_gene))

# 数据合并-------------------------------------------------------------------
#KEGG通路合并
KEGG_pathway$ID <- rownames(KEGG_pathway)
RC_metadata$ID <- rownames(RC_metadata)
harmonized_KEGG_pathway <- inner_join(KEGG_pathway,RC_metadata,by="ID")
harmonized_KEGG_pathway <- dplyr::select(harmonized_KEGG_pathway,c("ID","RC","age","gender","BMI","smoking","drinking","diabetes","hypertension","MET","DD",colnames(KEGG_pathway)))
#KEGG module合并
KEGG_module$ID <- rownames(KEGG_module)
RC_metadata$ID <- rownames(RC_metadata)
harmonized_KEGG_module <- inner_join(KEGG_module,RC_metadata,by="ID")
harmonized_KEGG_module <- dplyr::select(harmonized_KEGG_module,c("ID","RC","age","gender","BMI","smoking","drinking","diabetes","hypertension","MET","DD",colnames(KEGG_module)))
#KEGG KO合并
KEGG_gene$ID <- rownames(KEGG_gene)
RC_metadata$ID <- rownames(RC_metadata)
harmonized_KEGG_gene <- inner_join(KEGG_gene,RC_metadata,by="ID")
harmonized_KEGG_gene <- dplyr::select(harmonized_KEGG_gene,c("ID","RC","age","gender","BMI","smoking","drinking","diabetes","hypertension","MET","DD",colnames(KEGG_gene)))

# 广义线性模型差异分析 --------------------------------------------------------------------
source("C:/Users/Dell/Desktop/每日记录/20230629/广义线性模型变量筛选.R")
#KEGG通路
outcome <- "RC"
exposure <- colnames(harmonized_KEGG_pathway)[3:ncol(harmonized_KEGG_pathway)]
cov <- c("age","gender","BMI","smoking","drinking","diabetes","hypertension","MET","DD")
KEGG_pathway_glm_summary <- linear_model_feature_selection(harmonized_KEGG_pathway,outcome,exposure,cov)
#KEGG module
outcome <- "RC"
exposure <- colnames(harmonized_KEGG_module)[3:ncol(harmonized_KEGG_module)]
cov <- c("age","gender","BMI","smoking","drinking","diabetes","hypertension","MET","DD")
KEGG_module_glm_summary <- linear_model_feature_selection(harmonized_KEGG_module,outcome,exposure,cov)
#KEGG KO
outcome <- "RC"
exposure <- colnames(harmonized_KEGG_gene)[3:ncol(harmonized_KEGG_gene)]
cov <- c("age","gender","BMI","smoking","drinking","diabetes","hypertension","MET","DD")
KEGG_gene_glm_summary <- linear_model_feature_selection(harmonized_KEGG_gene,outcome,exposure,cov)

# KEGG gene富集分析 --------------------------------------------------------------
Diff_gene_list <- KEGG_gene_glm_summary%>%.[which(.[,"FDR"]<0.05),"ID"]
options(clusterProfiler.download.method = "wininet")
#以Gene Ratio为坐标
Enrichment_analysis <- enrichKEGG(gene = Diff_gene_list,
                       organism = "ko",
                       keyType = "kegg",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05)
Enrichment_table <- Enrichment_analysis@result
Enrichment_table <- Enrichment_table[Reduce(intersect,list(c(which(Enrichment_table$p.adjust<0.05),which(Enrichment_table$qvalue<0.05)))),]
Enrichment_table <- separate(data=Enrichment_table, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
Enrichment_table <- mutate(Enrichment_table, GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
#pdf("Enrichment_analysis.pdf",width = 8,height = 8)
ggplot(Enrichment_table,aes(GeneRatio,fct_reorder(factor(Description), GeneRatio))) + 
  geom_point(aes(size=Count,color=qvalue)) +
  scale_color_gradient(low="purple",high = "red") + 
  labs(color="pvalue",size="Count",x="GeneRatio",y="KEGG term",title="KEGG enrichment in Pathway")+
  theme_bw()
#dev.off()

#以Enrichment factor为坐标
Enrichment_table <- Enrichment_analysis@result
Enrichment_table <- separate(data=Enrichment_table, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
Enrichment_table <- separate(data=Enrichment_table, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
Enrichment_table <- mutate(Enrichment_table, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
Enrichment_table <- Enrichment_table[Reduce(intersect,list(c(which(Enrichment_table$p.adjust<0.05),which(Enrichment_table$qvalue<0.05)))),]
#pdf("Enrichment_analysis_EF.pdf",width = 6,height =8)
ggplot(Enrichment_table,aes(enrichment_factor,fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=qvalue)) +
  scale_color_gradient(low="purple",high = "red") + 
  labs(color="pvalue",size="Count",x="Enrichment Factor",y="KEGG term",title="KEGG enrichment in Pathway")+
  theme_bw()
#dev.off()

# KEGG module富集分析 ----------------------------------------------------------------
options(clusterProfiler.download.method = "wininet")
Diff_gene_list <- KEGG_gene_glm_summary%>%.[which(.[,"FDR"]<0.05),"ID"]

#以Gene Ratio为坐标
Enrichment_analysis <- enrichMKEGG(gene = Diff_gene_list,
                                  organism = "ko",
                                  keyType = "kegg",
                                  pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05)
Enrichment_table <- Enrichment_analysis@result
Enrichment_table <- Enrichment_table[Reduce(intersect,list(c(which(Enrichment_table$p.adjust<0.05),which(Enrichment_table$qvalue<0.05)))),]
Enrichment_table <- separate(data=Enrichment_table, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
Enrichment_table <- mutate(Enrichment_table, GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
#pdf("Enrichment_module_analysis.pdf",width = 8,height = 8)
ggplot(Enrichment_table,aes(GeneRatio,fct_reorder(factor(Description), GeneRatio))) + 
  geom_point(aes(size=Count,color=qvalue)) +
  scale_color_gradient(low="purple",high = "red") + 
  labs(color="pvalue",size="Count",x="GeneRatio",y="KEGG term",title="KEGG enrichment in Module")+
  theme_bw()
#dev.off()

#以Enrichment factor为坐标
Enrichment_table <- Enrichment_analysis@result
Enrichment_table <- separate(data=Enrichment_table, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
Enrichment_table <- separate(data=Enrichment_table, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
Enrichment_table <- mutate(Enrichment_table, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
Enrichment_table <- Enrichment_table[Reduce(intersect,list(c(which(Enrichment_table$p.adjust<0.05),which(Enrichment_table$qvalue<0.05)))),]
#pdf("Enrichment_module_analysis_EF.pdf",width = 6,height =8)
ggplot(Enrichment_table,aes(enrichment_factor,fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=qvalue)) +
  scale_color_gradient(low="purple",high = "red") + 
  labs(color="pvalue",size="Count",x="Enrichment Factor",y="KEGG term",title="KEGG enrichment in Module")+
  theme_bw()
#dev.off()

# 展示特定的KEGG通路 ------------------------------------------------------------
ko.data <- KEGG_gene_glm_summary[which(KEGG_gene_glm_summary$FDR<1),]
ko.data <- ko.data[,c(3,5)]
rownames(ko.data) <- ko.data$ID
ko.data$ID <- NA
ko.data <- t(ko.data)
ko.data <- ko.data[1,]
p <- pathview(gene.data = ko.data, pathway.id = "00270", species = "ko", out.suffix = "ko.data", kegg.native = T)

# 关联分析 -------------------------------------------------
metabolites$ID <- rownames(metabolites)
temp <- left_join(metabolites,harmonized_KEGG_module,by="ID")
temp <- temp[which(temp$ID%in%RC_metadata$ID),]
summary(lm(b333~temp$M00017,data = temp))
cor.test(temp$b333,temp$M00010,method = "spearman")

