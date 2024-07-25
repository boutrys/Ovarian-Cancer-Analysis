#/*****************************************************************************************
#  *
#  * Pipeline for Ovarian cancer analysis - Copyright (C) <2017-2023> <Université catholique de Louvain (UCLouvain)>
#  * 	
#  * List of the contributors to the development of Excalibur simulation: see LICENSE file.
#* Description and complete License: see LICENSE file.
#* 	
#  * This program (Excalibur simulation) is free software: 
#  * you can redistribute it and/or modify it under the terms of the 
#* GNU General Public License as published by the Free Software Foundation, 
#* either version 3 of the License, or (at your option) any later version.
#* 
#  * This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#* 
#  * You should have received a copy of the GNU General Public License
#* along with this program (see COPYING file).  If not, 
#* see <http://www.gnu.org/licenses/>.
#* 
#  *****************************************************************************************/
#  
#  /**
#  *
#  * @author Simon Boutry
#*
#  */

#' survival_analysis
#'
#' given output of data_pretreatment function, will perform a survival analysis
#' 
#' @param data tibble return by data_pretreatment function
#' @variable_for_survival variable of the form c(event, time) (Must be in that order!!) in order to perform survival analysis based on these criteria
#' @to_investigate_survival clinical annotation to investigate for the survival analysis
#' @param output where to store the results
#'
#' @return pvalue survival analysis for variable_for_survival based on all to_investigate_survival criteria
#'

survival_analysis <- function(data_survival,
                              variable_for_survival = c(),
                              to_investigate_survival = c(),
                              to_remove = NA,
                              output = getwd()){
  ###Load needed libraries
  require("survminer")
  require("survival")
  library(writexl)
  library(tidyverse)
  library(gtsummary)
  
  
  #clean data
  if(!is.na(to_remove)){
    data_survival <- data_survival[-which(data_survival$Patient %in% to_remove),]  
  }
  
  
  #TODO median follow-up time 
  data_survival_new <- data_survival
  data_survival_new$Group <- "BigDFS"
  data_survival_new$Group[which(data_survival_new$DFS <= 12)] <- "smallDFS"
  new_data <- data_survival_new
  new_data$Group <- "All"
  data_survival_new <- rbind(data_survival_new, new_data)
  data_survival_new$surv_object <- Surv(time = eval(parse(text=paste("data_survival_new$", variable_for_survival[2], sep = ""))),
                                        event = eval(parse(text=paste("data_survival_new$", variable_for_survival[1], sep = ""))))
  surv_object <- "surv_object"
  tmp_model <- paste0(surv_object, "~", "Group")
  fit0 <- survminer::surv_fit(as.formula(tmp_model), data = data_survival_new)
  median_followup <- tbl_survfit(fit0,
                                 probs = 0.5,
                                 label_header = "**Median survival (95% CI)**"
  )
  median_followup <- median_followup$table_body[2:4,5:6]
  colnames(median_followup)[2] <- "Median survival (95% CI)"
  writexl::write_xlsx(median_followup, paste(output, "median_followUp_originalData.xlsx", sep = ""))
  g <- ggsurvplot(fit0, pval = FALSE, risk.table = TRUE, conf.int = FALSE,
                  title=paste("Kaplan-Meier Curve based on ", variable_for_survival[2], sep = ""),
                  data = data_survival_new) +
    ylab("Percent followed") +
    xlab("Time (in month)")
  png(paste(output,"Median_followUp_originalData_", variable_for_survival[2], ".png", sep = ""),
      width = 880, height = 880)
  print(g, newpage = FALSE)
  dev.off()
  
  
  #Change DFS because package surminor use the opposit 
  new_DFS <- data_survival$censored_DFS
  data_survival$censored_DFS[which(new_DFS == 1)] <- 0
  data_survival$censored_DFS[which(new_DFS == 0)] <- 1
  
  
  #Create a survival object based on variable for survival
  data_survival$surv_object <- Surv(time = eval(parse(text=paste("data_survival$", variable_for_survival[2], sep = ""))),
                                    event = eval(parse(text=paste("data_survival$", variable_for_survival[1], sep = ""))))
  surv_object <- "surv_object"
  
  
  #first model 
  tmp_model <- paste0(surv_object, "~", 1)
  fit0 <- survminer::surv_fit(as.formula(tmp_model), data = data_survival)
  g <- ggsurvplot(fit0, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                  title=paste("Kaplan-Meier Curve based on ", variable_for_survival[2], " with Median: 16, 95% CI: 13, 29", sep = ""),
                  data = data_survival) +
    ylab("Disease Free Survival") +
    xlab("Time (in month)")
  png(paste(output,"survplot_for_", variable_for_survival[2], ".png", sep = ""),
      width = 880, height = 880)
  print(g, newpage = FALSE)
  dev.off()
  
  
  #General model survival with all parameter to_investigate_survival
  #TODO while not converge, keep removing smallest pvalue
  tmp_model <- paste0(surv_object, "~", paste(to_investigate_survival, collapse = "+"))
  fit_coxph <- coxph(as.formula(tmp_model), data = data_survival)
  forest_plot <- ggforest(fit_coxph, data = as.data.frame(data_survival)) #Hazar ratio forest plot
  pdf(paste(output,"forest_plot_hazard_ratio_of_general_model.pdf", sep = ""))
  print(forest_plot, newpage = FALSE)
  dev.off()
  
  tmp_summary <- summary(fit_coxph)
  coef <- as.data.frame(tmp_summary$coefficients)
  coef <- cbind(data.frame(Variable = rownames(coef)), coef)
  write_xlsx(coef, paste(output, "coefficient_full_model_coxRegression.xlsx", sep = ""))
  name_test <- c("Likelihood ratio", "Wald", "Score (logrank)")
  statistics <- data.frame(Test = name_test,
                           Statistic = c(tmp_summary$logtest[1], tmp_summary$sctest[1],  tmp_summary$waldtest[1]),
                           Degree_of_freedom = c(tmp_summary$logtest[2], tmp_summary$sctest[2],  tmp_summary$waldtest[2]),
                           Pvalue = c(tmp_summary$logtest[3], tmp_summary$sctest[3],  tmp_summary$waldtest[3]))
  write_xlsx(statistics, paste(output, "statistics_full_model_coxRegression.xlsx", sep = ""))
  
  
  #Generate plots and pvalues for each vairable to investigate for survival
  for (i in 1:length(to_investigate_survival)) {
    tmp_model <- paste0(surv_object, "~", to_investigate_survival[i])
    fit <- survminer::surv_fit(as.formula(tmp_model), data = data_survival)
    method_2 <- coxph(as.formula(tmp_model), data = data_survival)
    tmp <- summary(method_2)
    #store results
    if(i == 1){
      result <- data.frame(Variable = to_investigate_survival,
                           Likelihood_ratio_pvalue = NA,
                           Wald_pvalue = NA,
                           Score_logrank_pvalue = NA)
    }
    result[i,-1] <- c(tmp$logtest[3], tmp$sctest[3],  tmp$waldtest[3])
    
    
    tmp_label <- c(unique(eval(parse(text=paste("data_survival$", to_investigate_survival[i], sep = "")))))
    tmp_label <- tmp_label[which(!is.na(tmp_label))]
    #only for categorical variable with less than 5 categories
    if(length(tmp_label) < 5){
      g <- ggsurvplot(fit, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                      legend.title=gsub("_", " ", to_investigate_survival[i]),
                      title=paste("Kaplan-Meier Curve for ", gsub("_", " ", to_investigate_survival[i]), sep = ""),
                      data = data_survival) +
        ylab("Disease Free Survival") +
        xlab("Time (in month)")
      
      #tables <- ggsurvtable(fit, data =  data_survival, color = "strata",
      #                      y.text = FALSE)
      
      png(paste(output,"survplot_", to_investigate_survival[i], "_for_", variable_for_survival[2], ".png", sep = ""),
          width = 880, height = 880)
      print(g, newpage = FALSE)
      dev.off() 
    }
    
  }#end of for loop
  result$surv_object <- paste(variable_for_survival, collapse = " and ")
  
  #save results for pvalues
  write_xlsx(result, paste(output, "pvalue_survival_analysis_for_", variable_for_survival[2], ".xlsx", sep = ""))
  
  
  return(result)
}
