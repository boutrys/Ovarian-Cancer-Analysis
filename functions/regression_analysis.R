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

#' regression_analysis
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

regression_analysis <- function(data_regression,
                                patient_to_exclude = c(),
                                dependent_variable = c(),
                                restricted_model = c(),
                                full_model = c(),
                                output = getwd()){
  #libraries
  require(broom)
  library(rcompanion)
  library(writexl)
  require(nnet)
  library(regclass)
  #library(bestglm)
  
  ###Data pre treatment
  #remove 
  if(length(patient_to_exclude) > 0){
    data_regression <- data_regression[-which(data_regression$Patient %in% patient_to_exclude),]  
  }
  
  #Change DFS into categorical variable
  data_regression$DFS[which(data_regression$DFS < 13)] <- 0
  data_regression$DFS[which(data_regression$DFS > 12)] <- 1
  data_regression$DFS <- as.factor(data_regression$DFS)
  
  ###Choosing the family distribution
  for (i in 1:length(dependent_variable)) {
    tmp_dependent_value <- dependent_variable[i]
    nbr_level <- unique(eval(parse(text=paste("data_regression$", dependent_variable[i], sep = ""))))
    if(length(nbr_level) == 2){
      #Binary, best choice Binomial distribution and logistic function
      data_full <- data_regression[,which(colnames(data_regression) %in% c(full_model,tmp_dependent_value))]
      #Null model
      model <- paste0(tmp_dependent_value, "~", 1)
      model.null <- glm(as.formula(model), data = data_full, family = binomial(link="logit"))
      
      #Full model
      model <- paste0(tmp_dependent_value, "~", paste(full_model, collapse = "+"))
      model.full <- glm(as.formula(model), data = data_full, family = binomial(link="logit"))
      tmp_summary <- summary(model.full)
      coef <- as.data.frame(tmp_summary$coefficients)
      coef <- cbind(data.frame(Variable = rownames(coef)), coef)
      colnames(coef)[dim(coef)[2]] <- "P-value"
      write_xlsx(coef, paste(output, tmp_dependent_value, "_coefficient_full_model.xlsx", sep = ""))
      stat_full <- glance(model.full)
      stat_full <- stat_full[,-c(2,7)]
      stat_full$converged <- model.full$converged
      write_xlsx(stat_full, paste(output, tmp_dependent_value, "_statistics_full_model.xlsx", sep = ""))
      colinearity <- as.data.frame(VIF(model.full))
      write_xlsx(colinearity, paste(output, tmp_dependent_value, "_colinearity_full_model.xlsx", sep = ""))
      
      #restricted model
      model_restricted <- paste0(tmp_dependent_value, "~", paste(restricted_model, collapse = "+"))
      model.restricted <- glm(as.formula(model_restricted), data = data_full, family = binomial(link="logit"))
      tmp_summary_restricted <- summary(model.restricted)
      coef_restricted <- as.data.frame(tmp_summary_restricted$coefficients)
      coef_restricted <- cbind(data.frame(Variable = rownames(coef_restricted)), coef_restricted)
      colnames(coef_restricted)[dim(coef_restricted)[2]] <- "P-value"
      write_xlsx(coef_restricted, paste(output, tmp_dependent_value, "_coefficient_restricted_model.xlsx", sep = ""))
      stat_full <- glance(model.restricted)
      stat_full <- stat_full[,-c(2,7)]
      stat_full$converged <- model.restricted$converged
      write_xlsx(stat_full, paste(output, tmp_dependent_value, "_statistics_restricted_model.xlsx", sep = ""))
      colinearity <- as.data.frame(VIF(model.restricted))
      write_xlsx(colinearity, paste(output, tmp_dependent_value, "_colinearity_restricted_model.xlsx", sep = ""))
      
      
      #potential Best full model stepwise approach
      if(FALSE){
        step_model <- step(model.null, 
                           scope = list(upper = model.full),
                           test = "Chisq",
                           data = data_full)
        #potential test = "Rao","F", "LRT", "Chisq"
        tmp_best <- as.character(step_model$formula)
        potential <- strsplit(tmp_best[3], split = " \\+ ")[[1]]
        potential <- potential[-length(potential)]
        if(length(potential) > 0){
          model <- paste0(tmp_dependent_value, "~", paste(potential, collapse = "+"))
          model.full_best <- glm(as.formula(model), data = data_full, family = binomial(link="logit"))
          tmp_summary <- summary(model.full_best)
          coef <- as.data.frame(tmp_summary$coefficients)
          coef <- cbind(data.frame(Variable = rownames(coef)), coef)
          colnames(coef)[dim(coef)[2]] <- "P-value"
          write_xlsx(coef, paste(output, tmp_dependent_value, "_coefficient_best_full_model.xlsx", sep = ""))
          stat_full <- glance(model.full_best)
          stat_full <- stat_full[,-c(2,7)]
          stat_full$converged <- model.full_best$converged
          write_xlsx(stat_full, paste(output, tmp_dependent_value, "_statistics_best_full_model.xlsx", sep = ""))
        }
        
        #potential Best restricted model stepwise approach
        step_model <- step(model.null, 
                           scope = list(upper = model.restricted),
                           test = "Chisq",
                           data = data_full)
        #potential test = "Rao","F", "LRT", "Chisq"
        tmp_best <- as.character(step_model$formula)
        potential <- strsplit(tmp_best[3], split = " \\+ ")[[1]]
        #potential <- potential[-length(potential) TO AUTOMATIZE : if no pvalue ok then model stops at the best, not the one just before
        if(length(potential) > 0){
          model <- paste0(tmp_dependent_value, "~", paste(potential, collapse = "+"))
          model.restricted_best <- glm(as.formula(model), data = data_full, family = binomial(link="logit"))
          tmp_summary <- summary(model.restricted_best)
          coef <- as.data.frame(tmp_summary$coefficients)
          coef <- cbind(data.frame(Variable = rownames(coef)), coef)
          colnames(coef)[dim(coef)[2]] <- "P-value"
          write_xlsx(coef, paste(output, tmp_dependent_value, "_coefficient_best_restricted_model.xlsx", sep = ""))
          stat_full <- glance(model.restricted_best)
          stat_full <- stat_full[,-c(2,7)]
          stat_full$converged <- model.restricted_best$converged
          write_xlsx(stat_full, paste(output, tmp_dependent_value, "_statistics_best_restricted_model.xlsx", sep = ""))
        }
      }
      

      #TO DO carefull if one of the model do not exists it will crash
      #comparing models to assess fit statistics TO DO this works, but not within a function, why? 
      #res_compare <- compareGLM(model.null,model.full, model.restricted, model.full_best,model.restricted_best)
      #res_compare <- cbind(res_compare$Fit.criteria,res_compare$Models)
      #write_xlsx(res_compare, paste(output, "summary_model_comparisson.xlsx", sep = ""))
      
      #compare each model to the previous one using anova
      #anova(model.null,model.restricted,model.restricted_best, model.full_best)
      #anova(model.null,model.restricted_best, model.full_best,model.restricted)
      
      #bestglm package
      #colnames(data_full)[which(colnames(data_full) == dependent_variable[2])] <- "y"
      #data_full <- data.frame(data_full)
      #res <- bestglm(Xy = data_full,
      #               family = binomial,
      #               IC = "AIC", 
      #               method = "exhaustive")
      
      #Predicting values
      #library(caret)
      #data.omit <- na.omit(data_full)
      #to_train <- createDataPartition(data.omit$censored_DFS, p=0.8, list = FALSE)
      #training <- data.omit[to_train,]
      #testing <- data.omit[-to_train,]
      #model <- glm(censored_DFS ~ BRCA+FIGO_stage+AINS_per_op+reponse_per_op+Residual_disease+Anapath_response, data = training, family = "binomial")
      #probabilities <- model %>% predict(testing, type = "response")
      #predicted.classes <- ifelse(probabilities > 0.5, "1", "0")
      #mean(predicted.classes == testing$censored_DFS)
      
    }else if (length(nbr_level) == 3){
      #Categorical with more than two categories, mulitnomial regression
      data_full <- data_regression[,which(colnames(data_regression) %in% c(full_model,tmp_dependent_value))]
      #Null model
      model <- paste0(tmp_dependent_value, "~", 1)
      model.null <- multinom(as.formula(model), data = data_full, family = binomial(link="logit"))
      
      #Full model
      model <- paste0(tmp_dependent_value, "~", paste(full_model, collapse = "+"))
      model.full <- multinom(as.formula(model), data = data_full, family = binomial(link="logit"))
      tmp_summary <- summary(model.full)
      coef <- as.data.frame(tmp_summary$coefficients)
      coef <- cbind(data.frame(Variable = rownames(coef)), coef)
      write_xlsx(coef, paste(output, tmp_dependent_value, "_coefficient_full_model.xlsx", sep = ""))
      stat_full <- glance(model.full)
      write_xlsx(stat_full, paste(output, tmp_dependent_value, "_statistics_full_model.xlsx", sep = ""))
      #colinearity <- as.data.frame(VIF(model.full))
      #write_xlsx(colinearity, paste(output, tmp_dependent_value, "_colinearity_full_model.xlsx", sep = ""))
      #TODO colinearity do not work here
      
      #Restricted model
      model <- paste0(tmp_dependent_value, "~", paste(restricted_model, collapse = "+"))
      model.restricted <- multinom(as.formula(model), data = data_full, family = binomial(link="logit"))
      tmp_summary <- summary(model.restricted)
      coef <- as.data.frame(tmp_summary$coefficients)
      coef <- cbind(data.frame(Variable = rownames(coef)), coef)
      write_xlsx(coef, paste(output, tmp_dependent_value, "_coefficient_restricted_model.xlsx", sep = ""))
      stat_full <- glance(model.restricted)
      write_xlsx(stat_full, paste(output, tmp_dependent_value, "_statistics_restricted_model.xlsx", sep = ""))
      #colinearity <- as.data.frame(VIF(model.restricted))
      #write_xlsx(colinearity, paste(output, tmp_dependent_value, "_colinearity_restricted_model.xlsx", sep = ""))
      
    }else{
      #TO DO make difference in between continuous / continuous non-negative / discrete (counts)
      #For the moment we assume continuous non-negative
      #Linear regression with Gamma distribution (could also use inverse Gaussian)
    }
  }#endo of for loops 

  
  #return(result)
}