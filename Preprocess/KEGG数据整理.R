# 清空环境 ----------------------------------------------------------------
library("pacman")
p_unload(p_loaded(), character.only = TRUE)
rm(list = ls())
gc()

# 加载需要的R包 -----------------------------------------------------------------
library("dplyr")

# 设置工作路径并读取文件 ------------------------------------------------------------------
setwd("C:/Users/Dell/Desktop/每日记录/20230629/KEGG数据的整理/")
KEGG_fecal_id <- read.csv("KEGG_fecal_id.csv",check.names = F)
standard_id <- read.csv("standard_id.csv",check.names = F)
id_match <- inner_join(KEGG_fecal_id,standard_id,by=c("KEGG_fecal_id"="standard_fecal_id"))
#write.csv(id_match,"id_match.csv",row.names = F)

# 进行id转换 ------------------------------------------------------------
id_match <- read.csv("id_match.csv",check.names = F)
KEGG_pathway <- read.csv("huanan_KEGG.csv",check.names = F)
KEGG_module <- read.csv("huanan_MODULE.csv",check.names = F)
for (i in c(1:ncol(KEGG_pathway))) {
  if(colnames(KEGG_pathway)[i]%in%id_match$KEGG_fecal_id){
    colnames(KEGG_pathway)[i] <- id_match$standard_population_id[which(id_match$KEGG_fecal_id==colnames(KEGG_pathway)[i])]
  }
}
for (i in c(1:ncol(KEGG_module))) {
  if(colnames(KEGG_module)[i]%in%id_match$KEGG_fecal_id){
    colnames(KEGG_module)[i] <- id_match$standard_population_id[which(id_match$KEGG_fecal_id==colnames(KEGG_module)[i])]
  }
}
KEGG_pathway <- KEGG_pathway[,c(1:5,which(colnames(KEGG_pathway)%in%id_match$standard_population_id))]
KEGG_module <- KEGG_module[,c(1:3,which(colnames(KEGG_module)%in%id_match$standard_population_id))]
#write.csv(KEGG_pathway,"KEGG_pathway_extracted.csv",row.names = F)
#write.csv(KEGG_module,"KEGG_module_extracted.csv",row.names = F)
