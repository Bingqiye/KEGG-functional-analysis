linear_model_feature_selection <- function(dataframe,outcome,exposure,cov){
  temp <- dplyr::select(dataframe,c(outcome,cov,exposure))
  results_summary <- NULL
  for (i in c((sum(length(outcome),length(cov))+1):ncol(temp))) {
    temp_glm <- dplyr::select(dataframe,c(outcome,cov,colnames(temp)[i]))
    fomula <- paste(outcome,"~",".",sep = "")
    model <- glm(fomula,data = temp_glm)
    if((length(coef(summary(model)))/4)>(sum(length(outcome),length(cov)))){
      model_coef <- as.data.frame(t(coef(summary(model))[sum(length(outcome),length(cov))+1,]))
      model_coef$ID <- colnames(temp)[i]
      results_summary <- rbind(results_summary,model_coef)
    }else {
      model_coef <- data.frame("Estimate"=NA,"Std. Error"=NA,"t value"=NA,"Pr(>|t|)"=NA,"ID"=colnames(temp)[i],check.names = F)
      results_summary <- rbind(results_summary,model_coef)
    }
  }
  results_summary <- results_summary[order(results_summary$`Pr(>|t|)`,decreasing = F),]
  results_summary$FDR <- p.adjust(results_summary$`Pr(>|t|)`,method = "BH")
  return(results_summary)
}
