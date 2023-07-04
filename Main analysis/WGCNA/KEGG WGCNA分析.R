# 清空环境 ----------------------------------------------------------------
library("pacman")
p_unload(p_loaded(), character.only = TRUE)
rm(list = ls())
gc()

# 加载需要的R包 -----------------------------------------------------------------
library("dplyr")
library("vegan")
library("tidyverse")
library("ggplot2")
library("pathview")
library("DOSE")
library("WGCNA")
library("pheatmap")
library("psych")
library(doParallel) 
library(caret)
library(caretEnsemble)
library(MLmetrics)
library(data.table)
library(Boruta) 
library(tidyverse)
library(cvAUC)
library(ggplot2)
library(stringr)
library(pROC)
library(ggExtra)

# 设置工作路径并读取文件 ------------------------------------------------------------------
setwd("C:/Users/Dell/Desktop/每日记录/20230701/")
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

# 去除变量全为0值的样本 -------------------------------------------------------------
harmonized_KEGG_module <- harmonized_KEGG_module[,colSums(abs(harmonized_KEGG_module[,12:ncol(harmonized_KEGG_module)])) !=0]

# WGCNA分析 -----------------------------------------------------------------
sampleTree <- hclust(dist(harmonized_KEGG_module[,12:ncol(harmonized_KEGG_module)]), method ="average")
plot(sampleTree)
traitColors <- numbers2colors(harmonized_KEGG_module[,2], signed = FALSE)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels ="RC",
                    main ="Sample dendrogram and trait heatmap")
#选择软阈值
powers =c(c(1:10),seq(from = 12, to=20,by=2))
sft = pickSoftThreshold(harmonized_KEGG_module[,12:ncol(harmonized_KEGG_module)], powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
cex1 = 0.9
#绘图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="SoftThreshold(power)",ylab="ScaleFreeTopologyModelFit,signedR^2",type="n",
     main =paste("Scaleindependence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="SoftThreshold(power)",ylab="MeanConnectivity", type="n",
     main =paste("Meanconnectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],labels=powers, cex=cex1,col="red")
abline(h=100,col="red")
#WGCNA正式分析
net = blockwiseModules(harmonized_KEGG_module[,12:ncol(harmonized_KEGG_module)],power= 3,
                       TOMType ="signed", minModuleSize = 30,	
                       reassignThreshold = 0, mergeCutHeight = 0.25,	
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase ="RCTOM",
                       verbose = 3,randomSeed = 9999)
table(net$colors)
mergedColors <- labels2colors(net$colors)
# 绘制树状图和模块颜色图
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 每个KO所属的模块与颜色
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
# 每个样本对应的模块评分
MEs <- net$MEs
rownames(MEs) <- harmonized_KEGG_module$ID
MEs$ID <- rownames(MEs)

# 整合数据框，使用线性模型进行关联分析确定关键Module --------------------------------------------------------
RC_metadata_module <- inner_join(RC_metadata,MEs,by="ID")
source("C:/Users/Dell/Desktop/每日记录/20230629/广义线性模型变量筛选.R")
outcome <- "RC"
exposure <- colnames(RC_metadata_module)[38:ncol(RC_metadata_module)]
cov <- c("age","gender","BMI","smoking","drinking","diabetes","hypertension","MET","DD")
RC_metadata_module_glm_summary <- linear_model_feature_selection(RC_metadata_module,outcome,exposure,cov)
RC_metadata_module <- dplyr::select(RC_metadata_module,c(colnames(RC_metadata),RC_metadata_module_glm_summary$ID[which(RC_metadata_module_glm_summary$FDR<0.05)]))

# 根据得到的差异Module挑出相应的KEGG KO -----------------------------------------------
module_KO_list <- as.data.frame(moduleLabels)
module_KO_list$ME <- "ME"
module_KO_list$Type <- str_c(module_KO_list$ME,module_KO_list$moduleLabels)
module_KO_list$KO <- rownames(module_KO_list)
diff_module_KO_list <- dplyr::select(module_KO_list,c(KO,Type))
diff_module_KO_list <- diff_module_KO_list[which(diff_module_KO_list$Type%in%RC_metadata_module_glm_summary$ID[which(RC_metadata_module_glm_summary$FDR<0.05)]),]
diff_module_KO_list <- diff_module_KO_list[which(diff_module_KO_list$Type!="ME0"),]

# 对筛选得到的KO根据模块分别进行关联分析 ------------------------------------------------------
harmonized_KEGG_module <- dplyr::select(harmonized_KEGG_module,c(colnames(harmonized_KEGG_module)[which(colnames(harmonized_KEGG_module)%in%colnames(RC_metadata))],diff_module_KO_list$KO[which(diff_module_KO_list$Type=="ME2")]))
spearman_test <- corr.test(harmonized_KEGG_module[,12:103],harmonized_KEGG_module[,12:103],method = "spearman",adjust = "BH")
p <- pheatmap(spearman_test$r,fontsize_number=4,fontsize = 7,
         cellwidth = 7, cellheight = 7,
         cluster_rows = T,
         display_numbers = matrix(ifelse(spearman_test$p < 0.001, 
                                         "***",
                                         ifelse(spearman_test$p < 0.01 ,"**",
                                                ifelse(spearman_test$p<0.05 , "*"," "))
         ), nrow(spearman_test$p)),
         cluster_cols = T)
row_cluster <- cutree(p$tree_row, k=2)
newOrder <- harmonized_KEGG_module[,12:103][p$tree_row$order,]
cluster_results <- data.frame(ID=colnames(newOrder),cluster=row_cluster)
cluster1 <- cluster_results$ID[which(cluster_results$cluster==1)]
harmonized_KEGG_module_cluster1 <- dplyr::select(harmonized_KEGG_module,c(RC,cluster1,ID))

# 计算 in sclico score ------------------------------------------------------
source("interative_linear_modeling_with_continuous_outcome_without_feature_selection.R")
source("readin_interation_results_continuous_without_feature_selection.R")

interative_linear_modeling_with_continuous_outcome_without_feature_selection(pathway ="C:/Users/Dell/Desktop/每日记录/20230701/Output/",
                                                                             dataframe = harmonized_KEGG_module_cluster1,ID="ID",outcome = "RC",
                                                                             ninteration = 20,mproportion = 0.7,keep_list = cluster1)

readin_interation_results_continuous_without_feature_selection(pathway="C:/Users/Dell/Desktop/每日记录/20230701/Output/",
                                                               dataframe=harmonized_KEGG_module_cluster1, outcome="RC")

harmonized_KEGG_module_cluster1 <- read.table("C:/Users/Dell/Desktop/每日记录/20230701/Output/weighted_score.tsv",header = T)
harmonized_KEGG_module_cluster1$weighted_score <- -harmonized_KEGG_module_cluster1$weighted_score
summary(lm(RC~harmonized_KEGG_module_cluster1$weighted_score,data = harmonized_KEGG_module))


