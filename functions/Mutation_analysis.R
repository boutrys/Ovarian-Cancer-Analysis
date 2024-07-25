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

#' Mutation_analysis
#'
#' TODO description
#' Given results for general Facets from Facets_analysis() and clinical data and mutaiton data (mandatory column : sample, variant_type)
#' allow to analysis of regions of LOH/gain
#' plot these results and computes statistics and annotate clinical data with the different cni score 
#' 
#' @param group_name which group of cni score to be compared (TuPo TuPr) 
#' @param path_input where to retrieve clinical data and the cni score from different software
#' @param path_output where to store results
#' @param path_analysis where to store the new annotated clinical data
#'
#'

Mutation_analysis <- function(path_mutation_data = c(),
                              no_mutation = c(),
                              clinical_data_path = c(),
                              to_change = c(),
                              new_val = c(),
                              to_remove = c(),
                              path_output = getwd()){
  #libraries
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(ggpubr)
  library(ggplot2)
  library(cowplot)
  library(car)
  library(plotly)
  library(RColorBrewer)
  
  
  dir.create(path_output)
  #Load and pre treatment clinical data
  clinical_data <- read_xlsx(clinical_data_path)
  clinical_data$DFS[which(clinical_data$DFS < 13)] <- "=<12"
  clinical_data$DFS[which(clinical_data$DFS > 12)] <- ">12"
  
  #Load and pre treatment of mutation data from Highlander
  mutation_data <- read_tsv(path_mutation_data)
  mutation_data$sample <- gsub("VDB-", "", mutation_data$sample)
  if(length(to_change) > 0){
    for (i in 1:length(to_change)) {
      mutation_data$sample[which(mutation_data$sample == to_change[i])] <- new_val[i]
    } 
  }
  for (i in 1:length(to_remove)) {
    tmp_to_rm <- which(mutation_data$sample == to_remove[i])
    if(length(tmp_to_rm) > 0){
      mutation_data <- mutation_data[-tmp_to_rm,]
    }
  }
  to_rm_mnv <- which(mutation_data$variant_type == "MNV")
  if(length(to_rm_mnv) > 0){
    mutation_data <- mutation_data[-to_rm_mnv,]
  }
  mutation_data$variant_type[mutation_data$variant_type == "INS"] <- "INDEL"
  mutation_data$variant_type[mutation_data$variant_type == "DEL"] <- "INDEL"
  write_tsv(mutation_data, paste(path_output, "mutation_data_cleaned.tsv", sep = ""))
  
  
  #data pretreatment
  to_plot_mutation <- data.frame(patient = unique(mutation_data$sample))
  to_plot_mutation$SNVs <- NA
  to_plot_mutation$INDELs <- NA
  to_plot_mutation$group <- NA
  to_plot_mutation$Anapath_response <- NA
  to_plot_mutation$DFS <- NA
  to_plot_mutation$PFI <- NA
  for (i in 1:dim(to_plot_mutation)[1]) {
    tmp <- strsplit(to_plot_mutation$patient[i], split = "-")[[1]]
    to_plot_mutation$group[i] <- tmp[2]
    to_plot_mutation$Anapath_response[i] <- clinical_data$Anapath_response[which(clinical_data$Patient == tmp[1])]
    to_plot_mutation$DFS[i] <- clinical_data$DFS[which(clinical_data$Patient == tmp[1])]
    to_plot_mutation$PFI[i] <- clinical_data$PFI[which(clinical_data$Patient == tmp[1])]
    to_plot_mutation$SNVs[i] <- length(which(mutation_data$variant_type[which(mutation_data$sample == to_plot_mutation$patient[i])] == "SNV"))
    to_plot_mutation$INDELs[i] <- length(which(mutation_data$variant_type[which(mutation_data$sample == to_plot_mutation$patient[i])] == "INDEL"))
    to_plot_mutation$patient[i] <- tmp[1]
  }
  to_plot_mutation <- to_plot_mutation[-which(to_plot_mutation$group == "ASC"),] #TODO here we remove ASC, but we could include them
  
  #add rows for patient with no mutation 
  for (i in 1:length(no_mutation)) {
    tmp <- strsplit(no_mutation[i], split = "-")[[1]]
    if(tmp[1] %in% c("FETT", "KDAU")){
      new_row <- data.frame(patient = tmp[1], SNVs = 0, INDELs = 0, group = tmp[2], 
                            Anapath_response = clinical_data$Anapath_response[which(clinical_data$Patient == tmp[1])],
                            DFS = clinical_data$DFS[which(clinical_data$Patient == tmp[1])],
                            PFI = clinical_data$PFI[which(clinical_data$Patient == tmp[1])])
    }else{
      new_row <- to_plot_mutation[which(to_plot_mutation$patient == tmp[1]),]
    }
    new_row$group <- tmp[2]
    new_row$SNVs <- 0
    new_row$INDELs <- 0
    to_plot_mutation <- rbind(to_plot_mutation, new_row)
  }
  
  #change to Tumor Mutation Burden
  #TMB <- Ask raph
  TMB <- 33.755675
  to_plot_mutation$SNVs <- to_plot_mutation$SNVs / TMB
  to_plot_mutation$INDELs <- to_plot_mutation$INDELs / TMB
  
  
  ###Start loop for statistics and plot for each predictors
  to_plot_mutation$Anapath_response <- ordered(as.factor(to_plot_mutation$Anapath_response), c("Bad", "Intermediate", "Good"))
  to_plot_mutation$DFS <- ordered(as.factor(to_plot_mutation$DFS), c("=<12", ">12"))
  to_plot_mutation$PFI <- ordered(as.factor(to_plot_mutation$PFI), c("platinum resistant", "semi sensitive", "platinum sensitive"))
  to_plot_mutation$group <- ordered(as.factor(to_plot_mutation$group), c("TuPr", "TuPo"))
  var_predictor <- c("SNVs", "INDELs")
  group_name <- c("TuPr", "TuPo")
  group_name <- ordered(as.factor(group_name), c("TuPr", "TuPo"))
  nbr_group <- length(group_name)
  var_of_interest <- c("Anapath_response", "DFS", "PFI")
  y_axis_max <- max(c(to_plot_mutation$SNVs, to_plot_mutation$INDELs), na.rm = TRUE)
  y_axis_max <- y_axis_max * 0.5 + y_axis_max
  list_plot <- vector(mode = "list", length = length(var_predictor)*length(var_of_interest))
  count <- 1
  for (m in 1:length(var_predictor)) {
    data_plot <- to_plot_mutation
    colnames(data_plot)[which(colnames(data_plot) == var_predictor[m])] <- "predictor"
    #Loop go through all var of interest
    for (i in 1:length(var_of_interest)) {
      
      ###STATISTICS 
      value_var_of_interest <- unique(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))))
      value_var_of_interest <- value_var_of_interest[order(value_var_of_interest)]
      nbr_value_var_of_interest <- length(value_var_of_interest)
      res_stat <- vector("list", nbr_group)
      name_stat_test <- c()
      to_plot_y_max <- 0
      for (k in 1:nbr_group) {
        tmp_data <- data_plot[which(data_plot$group == group_name[k]),]
        colnames(tmp_data)[which(colnames(tmp_data) == var_of_interest[i])] <- "variable"
        #check normality
        normality_table <- data.frame(matrix(NA, ncol = 2, nrow = nbr_value_var_of_interest))
        colnames(normality_table) <- c("condition", "pvalue")
        for (j in 1:nbr_value_var_of_interest) {
          err <- try(tmp <- shapiro.test(tmp_data$predictor[which(tmp_data$variable == value_var_of_interest[j])])$p.value,
                     silent = TRUE)
          if(!is.character(err)){
            normality_table$pvalue[j] <- err
          }
          normality_table$condition[j] <- as.character(value_var_of_interest[j]) 
        }
        normality_to_add_name <- paste("normality", normality_table$condition, sep = "_")
        normality_to_add <- data.frame(matrix(normality_table$pvalue, ncol = length(normality_to_add_name), nrow = 1))
        colnames(normality_to_add) <- normality_to_add_name
        #reject normality assumptions if at least one pvalue < alpha = 0.05 or NA
        reject_normality <- (length(which(normality_table$pvalue < 0.05)) > 0 || length(which(is.na(normality_table$pvalue))) > 0)
        #Check homogeneity of variance
        if(reject_normality){
          if(nbr_value_var_of_interest == 2){
            #Fligner-Killeen test
            homo_var <- fligner.test(predictor ~ variable, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor , na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Rank-Sum test")
          }else if(nbr_value_var_of_interest > 2){
            #Levene’s test
            homo_var <- leveneTest(predictor ~ variable, data = tmp_data)$`Pr(>F)`[1]
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_value_var_of_interest,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(group_name[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Rank-Sum test")
            
          }
        }else{
          #comparing two groups normaly distributed
          if(nbr_value_var_of_interest == 2){
            #F-test
            homo_var <- var.test(predictor ~ variable, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is true 
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "T-test")
          }else if(nbr_value_var_of_interest > 2){
            #Barlett's test
            homo_var <- bartlett.test(predictor ~ variable, data = tmp_data)$p.value
            
            #TODO
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_value_var_of_interest,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(group_name[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "T-test")
          }
        }
        to_plot_y_max <- max(c(to_plot_y_max, y_max))
        res_stat[[k]] <- cbind(res_stat[[k]], normality_to_add, data.frame(homogeneity_variance = homo_var))
        
      }
      to_plot_stat <- res_stat[[1]]
      for (k in 2:nbr_group) {
        to_plot_stat <- rbind(to_plot_stat, res_stat[[k]])
      }
      to_plot_stat$name_test <- name_stat_test
      write_xlsx(to_plot_stat, paste(path_output, "Stat_", var_predictor[m], "_", var_of_interest[i], ".xlsx", sep = ""))
      to_plot_stat$p.adj <- round(to_plot_stat$p.adj, digits = 4)
      name_stat_test <- unique(name_stat_test)
      if(length(name_stat_test) > 1){
        print("different stat test used")
        name_stat_test <- ""
      }
      to_plot_stat$group <- ordered(as.factor(to_plot_stat$group), c("TuPr", "TuPo"))
      
      
      #Box plot for each variable of interest
      tmp_x_lab <- table(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))))
      g <- ggplot(data_plot, aes(y = predictor, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
        geom_boxplot() +
        geom_point() +
        xlab(var_of_interest[i]) +
        ylab(paste("TMB", var_predictor[m])) +
        scale_x_discrete(labels=names(tmp_x_lab)) +
        stat_pvalue_manual(to_plot_stat, label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_axis_max) +
        theme_classic()
      p <- g + facet_wrap(~group, nrow = 1)
      list_plot[[count]] <- p
      names(list_plot)[[count]] <- paste(var_predictor[m], var_of_interest[i], sep = "_")
      count <- count + 1
      fig <- ggplotly(p)
      #ggsave(p, filename = paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], ".png", sep = ""), height = 8, width = 12)
      #htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], ".html", sep = ""))
      
      data_plot_tupr <- data_plot[which(data_plot$group == "TuPr"),]
      tmp_x_lab_tupr <- table(eval(parse(text=paste("data_plot_tupr$",var_of_interest[i], sep = ""))))
      data_plot_tupo <- data_plot[which(data_plot$group == "TuPo"),]
      tmp_x_lab_tupo <- table(eval(parse(text=paste("data_plot_tupo$",var_of_interest[i], sep = ""))))
      plot1 <- ggplot(data_plot_tupr,aes(y = predictor, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
        geom_boxplot() +
        geom_point() +
        xlab(var_of_interest[i]) +
        ylab(paste("TMB", var_predictor[m])) +
        stat_pvalue_manual(to_plot_stat[which(to_plot_stat$group == "TuPr"),], label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_axis_max) +
        scale_x_discrete(labels=paste(names(tmp_x_lab_tupr), "=",tmp_x_lab_tupr)) +
        labs("TuPr") +
        theme_classic()
      plot2 <- ggplot(data_plot_tupo,aes(y = predictor, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
        geom_boxplot() +
        geom_point() +
        xlab(var_of_interest[i]) +
        ylab(paste("TMB", var_predictor[m])) +
        stat_pvalue_manual(to_plot_stat[which(to_plot_stat$group == "TuPo"),], label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_axis_max) +
        scale_x_discrete(labels=paste(names(tmp_x_lab_tupo), "=",tmp_x_lab_tupo)) +
        labs("TuPo") +
        theme_classic()
      new <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)
      ggsave(new, filename = paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], ".png", sep = ""), height = 8, width = 12)
      
      
    }#end Loop go through all var of interest
  }#end loop go through all predictors 
  saveRDS(list_plot, paste(path_output, "list_plot.rds", sep = ""))
  
  
  count_table <- table(to_plot_mutation$patient)
  to_plot_mutation$statut <- NA
  for (i in 1:length(count_table)) {
    if(count_table[i] == length(group_name)){
      to_plot_mutation$statut[which(to_plot_mutation$patient == names(count_table[i]))] <- "paired"
    }else{
      to_plot_mutation$statut[which(to_plot_mutation$patient == names(count_table[i]))] <- "not_paired"
    }
  }
  not_present_in_mutation_patient <- to_plot_mutation[which(to_plot_mutation$statut == "not_paired"),]
  write_xlsx(not_present_in_mutation_patient, paste(path_output, "Not_present_in_mutation_TuPrVSTuPo_analysis.xlsx", sep = ""))
  
  
  list_plot <- vector(mode = "list", length = length(var_predictor)*length(var_of_interest))
  count <- 1
  for (m in 1:length(var_predictor)) {
    data_plot <- to_plot_mutation[which(to_plot_mutation$statut == "paired"),]
    colnames(data_plot)[which(colnames(data_plot) == var_predictor[m])] <- "predictor"
    #Loop go through all var of interest
    for (i in 1:length(var_of_interest)) {
      
      ###STATISTICS 
      value_var_of_interest <- unique(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))))
      value_var_of_interest <- value_var_of_interest[order(value_var_of_interest)]
      nbr_value_var_of_interest <- length(value_var_of_interest)
      res_stat <- vector("list", nbr_value_var_of_interest)
      name_stat_test <- c()
      to_plot_y_max <- 0
      for (k in 1:nbr_value_var_of_interest) {
        tmp_data <- data_plot[which(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))) == value_var_of_interest[k]),]
        #check normality
        normality_table <- data.frame(matrix(NA, ncol = 2, nrow = nbr_group))
        colnames(normality_table) <- c("condition", "pvalue")
        for (j in 1:nbr_group) {
          err <- try(tmp <- shapiro.test(tmp_data$predictor[which(tmp_data$group == group_name[j])])$p.value,
                     silent = TRUE)
          if(!is.character(err)){
            normality_table$pvalue[j] <- err
          }
          normality_table$condition[j] <- as.character(group_name[j]) 
        }
        normality_to_add_name <- paste("normality", normality_table$condition, sep = "_")
        normality_to_add <- data.frame(matrix(normality_table$pvalue, ncol = length(normality_to_add_name), nrow = 1))
        colnames(normality_to_add) <- normality_to_add_name
        #reject normality assumptions if at least one pvalue < alpha = 0.05 or NA
        reject_normality <- (length(which(normality_table$pvalue < 0.05)) > 0 || length(which(is.na(normality_table$pvalue))) > 0)
        #Check homogeneity of variance
        if(reject_normality){
          if(nbr_group == 2){
            #Fligner-Killeen test
            homo_var <- fligner.test(predictor ~ group, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor , na.rm = TRUE)
            y_max <- y_max * 0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ group, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = TRUE) %>%
              mutate(y.position = y_max)
            if(dim(tmp_res_stat)[1] < 1){
              tmp_res_stat <- data.frame(.y. = "predictor", group1 = "TuPr", group2 = "TuPo", p = "Na", p.adj = 0, p.format = "Na", p.signif = "Na", method = "Na", y.position =  y_max)
            }
            res_stat[[k]] <- cbind(data.frame(group = value_var_of_interest[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Signed-rank test")
          }else if(nbr_group > 2){
            #Levene’s test
            homo_var <- leveneTest(predictor ~ group, data = tmp_data)$`Pr(>F)`[1]
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_group) + y_max
            tmp_res_stat <- compare_means(predictor ~ group, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = TRUE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_group,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(value_var_of_interest[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Signed-rank test")
            
          }
        }else{
          #comparing two groups normaly distributed
          if(nbr_group == 2){
            #F-test
            homo_var <- var.test(predictor ~ group, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is true 
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ group, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = TRUE) %>%
              mutate(y.position = y_max)
            if(dim(tmp_res_stat)[1] < 1){
              tmp_res_stat <- data.frame(.y. = "predictor", group1 = "TuPr", group2 = "TuPo", p = "Na", p.adj = 0, p.format = "Na", p.signif = "Na", method = "Na", y.position =  y_max)
            }
            res_stat[[k]] <- cbind(data.frame(group = value_var_of_interest[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Paired t-test")
          }else if(nbr_group > 2){
            #Barlett's test
            homo_var <- bartlett.test(predictor ~ group, data = tmp_data)$p.value
            
            #TODO
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_group) + y_max
            tmp_res_stat <- compare_means(predictor ~ group, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = TRUE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_group,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(value_var_of_interest[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Paired T-test")
          }
        }
        to_plot_y_max <- max(c(to_plot_y_max, y_max))
        res_stat[[k]] <- cbind(res_stat[[k]], normality_to_add, data.frame(homogeneity_variance = homo_var))
        
      }
      to_plot_stat <- res_stat[[1]]
      for (k in 2:nbr_value_var_of_interest) {
        to_plot_stat <- rbind(to_plot_stat, res_stat[[k]])
      }
      to_plot_stat$name_test <- name_stat_test
      write_xlsx(to_plot_stat, paste(path_output, "Stat_", var_predictor[m], "_", var_of_interest[i], "_TuPrVSTuPo.xlsx", sep = ""))
      to_plot_stat$p.adj <- round(to_plot_stat$p.adj, digits = 4)
      name_stat_test <- unique(name_stat_test)
      if(length(name_stat_test) > 1){
        print("different stat test used")
        name_stat_test <- ""
      }
      if(i == 1){
        to_plot_stat$group <- ordered(as.factor(to_plot_stat$group), c("Bad", "Intermediate", "Good"))
        colnames(to_plot_stat)[1] <- "Anapath_response"
      }
      if(i == 2){
        to_plot_stat$group <- ordered(as.factor(to_plot_stat$group), c("=<12", ">12"))
        colnames(to_plot_stat)[1] <- "DFS"
      }
      if(i == 3){
        to_plot_stat$group <- ordered(as.factor(to_plot_stat$group), c("platinum resistant", "semi sensitive", "platinum sensitive"))
        colnames(to_plot_stat)[1] <- "PFI"
      }
      
      
      #Box plot for each variable of interest
      tmp_x_lab <- table(data_plot$group)
      g <- ggplot(data_plot, aes(y = predictor, x = group, label = patient)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 2, width = 0.1) +
        stat_pvalue_manual(to_plot_stat, label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylab(paste("TMB", var_predictor[m])) +
        scale_x_discrete(labels=names(tmp_x_lab)) +
        ylim(0,y_axis_max) +
        theme_classic()
      p <- g + facet_wrap(~eval(parse(text=paste(var_of_interest[i]))), nrow = 1) +
        labs(title = var_of_interest[i]) 
      list_plot[[count]] <- p
      names(list_plot)[[count]] <- paste(var_predictor[m], var_of_interest[i], sep = "_")
      count <- count + 1
      fig <- ggplotly(p)
      ggsave(p, filename = paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], "_TuPrVSTuPo.png", sep = ""), height = 8, width = 12)
      htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], "_TuPrVSTuPo.html", sep = ""))
      
    }#end Loop go through all var of interest
  }#end loop go through all predictors 
  saveRDS(list_plot, paste(path_output, "list_plot_TuPrVSTuPo.rds", sep = ""))
  
  
  
  #Ploting the allelic_depth_proportion_alt TuPr vs TuPo
  af_table <- mutation_data[-grep(mutation_data$sample, pattern = "ASC"),] 
  only_name <- gsub("-TuPr", "", af_table$sample)
  only_name <- gsub("-TuPo", "", only_name)
  af_table$id_variant <- paste(af_table$chr, af_table$pos, af_table$reference, af_table$alternative, af_table$gene_symbol, only_name, sep = "_")
  
  to_plot_af <- data.frame(id = unique(af_table$id_variant))
  to_plot_af$TuPr_allelic_depth_proportion_alt <- 0
  to_plot_af$TuPr_allelic_depth_proportion_ref <- 1
  to_plot_af$TuPo_allelic_depth_proportion_alt <- 0
  to_plot_af$TuPo_allelic_depth_proportion_ref <- 1
  to_plot_af$Enrich_proportion <- NA
  to_plot_af$TuPr_allelic_depth_alt <- 0
  to_plot_af$TuPr_allelic_depth_ref <- 0
  to_plot_af$TuPo_allelic_depth_alt <- 0
  to_plot_af$TuPo_allelic_depth_ref <- 0
  to_plot_af$Enrich <- NA
  count <- 0
  for (i in 1:dim(to_plot_af)[1]) {
    tmp <- af_table[which(af_table$id_variant == to_plot_af$id[i]),]
    #allelic_depth_proportion_alt
    tmp_TuPr <- tmp$allelic_depth_proportion_alt[grep(tmp$sample, pattern = "TuPr")]
    tmp_TuPo <- tmp$allelic_depth_proportion_alt[grep(tmp$sample, pattern = "TuPo")]
    if(length(tmp_TuPr) > 0){
      to_plot_af$TuPr_allelic_depth_proportion_alt[i] <- tmp_TuPr
    }
    if(length(tmp_TuPo) > 0){
      to_plot_af$TuPo_allelic_depth_proportion_alt[i] <- tmp_TuPo
    }
    if(dim(tmp)[1] > 2){
      count <- count + 1
      print(i)
    }
    if(to_plot_af$TuPo_allelic_depth_proportion_alt[i] > to_plot_af$TuPr_allelic_depth_proportion_alt[i]){
      to_plot_af$Enrich_proportion[i] <- "Yes"
    }else{
      to_plot_af$Enrich_proportion[i] <- "No"
    }
    #allelic_depth_proportion_ref
    tmp_TuPr <- tmp$allelic_depth_proportion_ref[grep(tmp$sample, pattern = "TuPr")]
    tmp_TuPo <- tmp$allelic_depth_proportion_ref[grep(tmp$sample, pattern = "TuPo")]
    if(length(tmp_TuPr) > 0){
      to_plot_af$TuPr_allelic_depth_proportion_ref[i] <- tmp_TuPr
    }
    if(length(tmp_TuPo) > 0){
      to_plot_af$TuPo_allelic_depth_proportion_ref[i] <- tmp_TuPo
    }
    if(dim(tmp)[1] > 2){
      count <- count + 1
      print(i)
    }
    
    #Same but for allelic depth alt and ref
    tmp_TuPr <- tmp$allelic_depth_alt[grep(tmp$sample, pattern = "TuPr")]
    tmp_TuPo <- tmp$allelic_depth_alt[grep(tmp$sample, pattern = "TuPo")]
    if(length(tmp_TuPr) > 0){
      to_plot_af$TuPr_allelic_depth_alt[i] <- tmp_TuPr
    }
    if(length(tmp_TuPo) > 0){
      to_plot_af$TuPo_allelic_depth_alt[i] <- tmp_TuPo
    }
    if(dim(tmp)[1] > 2){
      count <- count + 1
      print(i)
    }
    if(to_plot_af$TuPo_allelic_depth_proportion_alt[i] > to_plot_af$TuPr_allelic_depth_proportion_alt[i]){
      to_plot_af$Enrich[i] <- "Yes"
    }else{
      to_plot_af$Enrich[i] <- "No"
    }
    #allelic_depth_ref
    tmp_TuPr <- tmp$allelic_depth_ref[grep(tmp$sample, pattern = "TuPr")]
    tmp_TuPo <- tmp$allelic_depth_ref[grep(tmp$sample, pattern = "TuPo")]
    if(length(tmp_TuPr) > 0){
      to_plot_af$TuPr_allelic_depth_ref[i] <- tmp_TuPr
    }
    if(length(tmp_TuPo) > 0){
      to_plot_af$TuPo_allelic_depth_ref[i] <- tmp_TuPo
    }
    if(dim(tmp)[1] > 2){
      count <- count + 1
      print(i)
    }
    
  }
  
  
  af_table$pvalue_enrich <- NA
  af_table$enrich <- NA
  to_plot_af$pvalue_enrich <- NA
  #Chisquar to compute enrichissment and annotate mutation data (af_table = mutation data without ASC)
  for(i in 1:dim(to_plot_af)[1]){
    contingency_matrix = matrix(c(to_plot_af$TuPr_allelic_depth_alt[i], to_plot_af$TuPr_allelic_depth_ref[i], 
                                  to_plot_af$TuPo_allelic_depth_alt[i], to_plot_af$TuPo_allelic_depth_ref[i]),
                                nrow=2, ncol=2,  byrow = TRUE)
    result <- chisq.test(contingency_matrix, simulate.p.value = FALSE)
    pvalue <- result$p.value
    #pvalue <- format(pvalue, scientific=FALSE)
    if(pvalue == "NaN"){
      if((to_plot_af$TuPr_allelic_depth_proportion_alt[i] == 0 && to_plot_af$TuPo_allelic_depth_proportion_alt[i] > 0.1) || 
         (to_plot_af$TuPr_allelic_depth_proportion_alt[i] > 0.1 && to_plot_af$TuPo_allelic_depth_proportion_alt[i] == 0)){
        to_plot_af$pvalue_enrich[i] <- "Significant"
        pvalue <- 0.0000000000000000009
      }else{
        to_plot_af$pvalue_enrich[i] <- "No"
        pvalue <- 1
      }
      
    }else{
      if(pvalue > 0.05){
        to_plot_af$pvalue_enrich[i] <- "No"
      }else if(pvalue <= 0.05){
        to_plot_af$pvalue_enrich[i] <- "Significant"
      }
    }
    af_table$enrich[which(af_table$id_variant == to_plot_af$id[i])] <- to_plot_af$Enrich[i]
    af_table$pvalue_enrich[which(af_table$id_variant == to_plot_af$id[i])] <- pvalue
  }
  
  to_plot_af$Statut <- paste(to_plot_af$Enrich, to_plot_af$pvalue_enrich, sep = "_")
  to_plot_af$Statut[which(to_plot_af$Statut == "No_No")] <- "not significant"
  to_plot_af$Statut[which(to_plot_af$Statut == "No_Significant")] <- "Depleted"
  to_plot_af$Statut[which(to_plot_af$Statut == "Yes_No")] <- "not significant"
  to_plot_af$Statut[which(to_plot_af$Statut == "Yes_Significant")] <- "Enriched"
  g <- ggplot(to_plot_af, aes(x=TuPo_allelic_depth_proportion_alt, y=TuPr_allelic_depth_proportion_alt), col=c("#FF0000", "#0033FF", "#666666")) +
    geom_jitter(aes(color = Statut)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = "Allelic depth proportion ALT") +
    theme_classic()
  ggsave(g, filename = paste(path_output, "allelic_depth_proportion_alt_TuPr_vs_TuPo.png", sep = ""), 
         height = 8, width = 8)
  
  
  #TODO same plot as above but coloring based on DFS
  new_to_plot <- to_plot_af
  new_to_plot$DFS <- NA
  for (i in 1:dim(new_to_plot)[1]) {
    new_to_plot$DFS[i] <- clinical_data$DFS[which(clinical_data$Patient == strsplit(new_to_plot$id[i], split = "_")[[1]][6])]
  }
  new_to_plot$DFS <- ordered(as.factor(new_to_plot$DFS), c("=<12", ">12"))
  g <- ggplot(new_to_plot, aes(x=TuPo_allelic_depth_proportion_alt, y=TuPr_allelic_depth_proportion_alt)) +
    geom_jitter(aes(color = DFS)) +
    scale_color_manual(values = c("#33A02C", "#FB9A99")) +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = "Allelic depth proportion ALT") +
    theme_classic()
  ggsave(g, filename = paste(path_output, "allelic_depth_proportion_alt_TuPr_vs_TuPo_DFS_color.png", sep = ""), 
         height = 8, width = 8)
  
  
  af_table <- af_table %>%
    mutate(Statut = case_when(
      pvalue_enrich <= 0.05 & enrich == "No" ~ "Depleted",
      pvalue_enrich <= 0.05 & enrich == "Yes" ~ "Enriched",
      TRUE ~ NA_character_
    ))
  write_tsv(af_table, paste(path_output, "mutation_data_withoutASC_annotEnrich.tsv", sep = ""))
  
  sign_pval <- filter(af_table, pvalue_enrich <= 0.05)
  if(dim(sign_pval)[1] > 0){
    write_tsv(sign_pval, paste(path_output, "Only_sig_enrich_mutation_data.tsv", sep = ""))
  }
  
  write_tsv(to_plot_af, paste(path_output, "table_plot_allelic_depth_prop.tsv", sep = ""))
  
  
  
  ### enrichment /depleted plots 
  #data pre treatment
  enrich_mutation_data <- af_table[which(!af_table$pvalue_enrich == "NaN"),]
  #enrich_mutation_data$qvalue_enrich <- p.adjust(enrich_mutation_data$pvalue_enrich, method = "BH")
  enrich_mutation_data$qvalue_enrich <- enrich_mutation_data$pvalue_enrich
  
  to_plot_mutation$Nbr_enrich <- 0
  to_plot_mutation$Nbr_depleted <- 0
  to_plot_mutation$Nbr_total <- 0
  to_plot_mutation$Nbr_variant_unique <- 0
  for (i in 1:dim(to_plot_mutation)[1]) {
    tmp_id <- grep(to_plot_af$id, pattern = to_plot_mutation$patient[i])
    to_plot_mutation$Nbr_enrich[i] <- length(which(to_plot_af$Statut[tmp_id] == "Enriched"))
    to_plot_mutation$Nbr_depleted[i] <- length(which(to_plot_af$Statut[tmp_id] == "Depleted"))
    to_plot_mutation$Nbr_variant_unique[i] <- length(grep(to_plot_mutation$patient[i], to_plot_af$id))
  }
  to_plot_mutation$Nbr_total <- to_plot_mutation$Nbr_enrich + to_plot_mutation$Nbr_depleted
  to_plot_mutation$fraction_nbr_enrich <- to_plot_mutation$Nbr_enrich / to_plot_mutation$Nbr_variant_unique
  to_plot_mutation$fraction_nbr_depleted <- to_plot_mutation$Nbr_depleted / to_plot_mutation$Nbr_variant_unique
  to_plot_mutation$fraction_nbr_total <- to_plot_mutation$Nbr_total / to_plot_mutation$Nbr_variant_unique
  write_xlsx(to_plot_mutation, paste(path_output, "summary_mutation_analysis.xlsx", sep = ""))
  
  
  var_predictor <- c("Nbr_enrich", "Nbr_depleted", "Nbr_total")
  group_name <- c("TuPr")
  nbr_group <- length(group_name)
  var_of_interest <- c("Anapath_response", "DFS", "PFI")
  y_axis_max <- max(c(to_plot_mutation$Nbr_enrich, to_plot_mutation$Nbr_depleted), na.rm = TRUE)
  y_axis_max <- y_axis_max * 0.5 + y_axis_max
  list_plot <- vector(mode = "list", length = length(var_predictor)*length(var_of_interest))
  count <- 1
  for (m in 1:length(var_predictor)) {
    data_plot <- to_plot_mutation[which(to_plot_mutation$group == "TuPr"),] #TODO only take paired sample, we took TuPr because NACH do not have it
    colnames(data_plot)[which(colnames(data_plot) == var_predictor[m])] <- "predictor"
    #Loop go through all var of interest
    for (i in 1:length(var_of_interest)) {
      
      ###STATISTICS 
      value_var_of_interest <- unique(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))))
      value_var_of_interest <- value_var_of_interest[order(value_var_of_interest)]
      nbr_value_var_of_interest <- length(value_var_of_interest)
      res_stat <- vector("list", nbr_group)
      name_stat_test <- c()
      to_plot_y_max <- 0
      for (k in 1:nbr_group) {
        tmp_data <- data_plot[which(data_plot$group == group_name[k]),]
        colnames(tmp_data)[which(colnames(tmp_data) == var_of_interest[i])] <- "variable"
        #check normality
        normality_table <- data.frame(matrix(NA, ncol = 2, nrow = nbr_value_var_of_interest))
        colnames(normality_table) <- c("condition", "pvalue")
        for (j in 1:nbr_value_var_of_interest) {
          err <- try(tmp <- shapiro.test(tmp_data$predictor[which(tmp_data$variable == value_var_of_interest[j])])$p.value,
                     silent = TRUE)
          if(!is.character(err)){
            normality_table$pvalue[j] <- err
          }
          normality_table$condition[j] <- as.character(value_var_of_interest[j]) 
        }
        normality_to_add_name <- paste("normality", normality_table$condition, sep = "_")
        normality_to_add <- data.frame(matrix(normality_table$pvalue, ncol = length(normality_to_add_name), nrow = 1))
        colnames(normality_to_add) <- normality_to_add_name
        #reject normality assumptions if at least one pvalue < alpha = 0.05 or NA
        reject_normality <- (length(which(normality_table$pvalue < 0.05)) > 0 || length(which(is.na(normality_table$pvalue))) > 0)
        #Check homogeneity of variance
        if(reject_normality){
          if(nbr_value_var_of_interest == 2){
            #Fligner-Killeen test
            homo_var <- fligner.test(predictor ~ variable, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor , na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Rank-Sum test")
          }else if(nbr_value_var_of_interest > 2){
            #Levene’s test
            homo_var <- leveneTest(predictor ~ variable, data = tmp_data)$`Pr(>F)`[1]
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_value_var_of_interest,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(group_name[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Rank-Sum test")
            
          }
        }else{
          #comparing two groups normaly distributed
          if(nbr_value_var_of_interest == 2){
            #F-test
            homo_var <- var.test(predictor ~ variable, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is true 
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "T-test")
          }else if(nbr_value_var_of_interest > 2){
            #Barlett's test
            homo_var <- bartlett.test(predictor ~ variable, data = tmp_data)$p.value
            
            #TODO
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_value_var_of_interest,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(group_name[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "T-test")
          }
        }
        to_plot_y_max <- max(c(to_plot_y_max, y_max))
        res_stat[[k]] <- cbind(res_stat[[k]], normality_to_add, data.frame(homogeneity_variance = homo_var))
        
      }
      to_plot_stat <- res_stat[[1]]
      to_plot_stat$name_test <- name_stat_test
      write_xlsx(to_plot_stat, paste(path_output, "Stat_", var_predictor[m], "_", var_of_interest[i], ".xlsx", sep = ""))
      to_plot_stat$p.adj <- round(to_plot_stat$p.adj, digits = 4)
      name_stat_test <- unique(name_stat_test)
      if(length(name_stat_test) > 1){
        print("different stat test used")
        name_stat_test <- ""
      }
      
      
      #Box plot for each variable of interest
      tmp_x_lab <- table(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))))
      g <- ggplot(data_plot, aes(y = predictor, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
        geom_boxplot() +
        geom_point() +
        xlab(var_of_interest[i]) +
        ylab(var_predictor[m]) +
        scale_x_discrete(labels=paste(names(tmp_x_lab), "=",tmp_x_lab)) +
        stat_pvalue_manual(to_plot_stat, label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_axis_max) +
        theme_classic()
      p <- g 
      list_plot[[count]] <- p
      names(list_plot)[[count]] <- paste(var_predictor[m], var_of_interest[i], sep = "_")
      count <- count + 1
      fig <- ggplotly(p)
      ggsave(p, filename = paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], ".png", sep = ""), height = 8, width = 12)
      htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], ".html", sep = ""))
      
    }#end Loop go through all var of interest
  }#end loop go through all predictors 
  saveRDS(list_plot, paste(path_output, "enrich_list_plot.rds", sep = ""))
  
  
  #same plot but for fraction nbr
  var_predictor <- c("fraction_nbr_enrich", "fraction_nbr_depleted", "fraction_nbr_total")
  group_name <- c("TuPr")
  nbr_group <- length(group_name)
  var_of_interest <- c("Anapath_response", "DFS", "PFI")
  y_axis_max <- max(to_plot_mutation$fraction_nbr_total, na.rm = TRUE)
  y_axis_max <- y_axis_max * 0.5 + y_axis_max
  list_plot <- vector(mode = "list", length = length(var_predictor)*length(var_of_interest))
  count <- 1
  for (m in 1:length(var_predictor)) {
    data_plot <- to_plot_mutation[which(to_plot_mutation$group == "TuPr"),] #TODO only take paired sample, we took TuPr because NACH do not have it
    colnames(data_plot)[which(colnames(data_plot) == var_predictor[m])] <- "predictor"
    #Loop go through all var of interest
    for (i in 1:length(var_of_interest)) {
      
      ###STATISTICS 
      value_var_of_interest <- unique(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))))
      value_var_of_interest <- value_var_of_interest[order(value_var_of_interest)]
      nbr_value_var_of_interest <- length(value_var_of_interest)
      res_stat <- vector("list", nbr_group)
      name_stat_test <- c()
      to_plot_y_max <- 0
      for (k in 1:nbr_group) {
        tmp_data <- data_plot[which(data_plot$group == group_name[k]),]
        colnames(tmp_data)[which(colnames(tmp_data) == var_of_interest[i])] <- "variable"
        #check normality
        normality_table <- data.frame(matrix(NA, ncol = 2, nrow = nbr_value_var_of_interest))
        colnames(normality_table) <- c("condition", "pvalue")
        for (j in 1:nbr_value_var_of_interest) {
          err <- try(tmp <- shapiro.test(tmp_data$predictor[which(tmp_data$variable == value_var_of_interest[j])])$p.value,
                     silent = TRUE)
          if(!is.character(err)){
            normality_table$pvalue[j] <- err
          }
          normality_table$condition[j] <- as.character(value_var_of_interest[j]) 
        }
        normality_to_add_name <- paste("normality", normality_table$condition, sep = "_")
        normality_to_add <- data.frame(matrix(normality_table$pvalue, ncol = length(normality_to_add_name), nrow = 1))
        colnames(normality_to_add) <- normality_to_add_name
        #reject normality assumptions if at least one pvalue < alpha = 0.05 or NA
        reject_normality <- (length(which(normality_table$pvalue < 0.05)) > 0 || length(which(is.na(normality_table$pvalue))) > 0)
        #Check homogeneity of variance
        if(reject_normality){
          if(nbr_value_var_of_interest == 2){
            #Fligner-Killeen test
            homo_var <- fligner.test(predictor ~ variable, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor , na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Rank-Sum test")
          }else if(nbr_value_var_of_interest > 2){
            #Levene’s test
            homo_var <- leveneTest(predictor ~ variable, data = tmp_data)$`Pr(>F)`[1]
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "wilcox.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_value_var_of_interest,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(group_name[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "Wilcoxon Rank-Sum test")
            
          }
        }else{
          #comparing two groups normaly distributed
          if(nbr_value_var_of_interest == 2){
            #F-test
            homo_var <- var.test(predictor ~ variable, data = tmp_data)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is true 
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "T-test")
          }else if(nbr_value_var_of_interest > 2){
            #Barlett's test
            homo_var <- bartlett.test(predictor ~ variable, data = tmp_data)$p.value
            
            #TODO
            y_max <- max(tmp_data$predictor, na.rm = TRUE)
            y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
            tmp_res_stat <- compare_means(predictor ~ variable, 
                                          data = tmp_data,
                                          method = "t.test",
                                          paired = FALSE) %>%
              mutate(y.position = y_max)
            #Keep only comparisson to reference group and adjust pvalue using Benjamini-Hochberg
            tmp_res_stat <- tmp_res_stat[-nbr_value_var_of_interest,-7] #we remove p.signif
            tmp_res_stat$p.adj <- p.adjust(tmp_res_stat$p, method = "BH")
            res_stat[[k]] <- cbind(data.frame(group = rep(group_name[k], dim(tmp_res_stat)[1])), tmp_res_stat)
            name_stat_test <- c(name_stat_test, "T-test")
          }
        }
        to_plot_y_max <- max(c(to_plot_y_max, y_max))
        res_stat[[k]] <- cbind(res_stat[[k]], normality_to_add, data.frame(homogeneity_variance = homo_var))
        
      }
      to_plot_stat <- res_stat[[1]]
      to_plot_stat$name_test <- name_stat_test
      write_xlsx(to_plot_stat, paste(path_output, "Stat_", var_predictor[m], "_", var_of_interest[i], ".xlsx", sep = ""))
      to_plot_stat$p.adj <- round(to_plot_stat$p.adj, digits = 4)
      name_stat_test <- unique(name_stat_test)
      if(length(name_stat_test) > 1){
        print("different stat test used")
        name_stat_test <- ""
      }
      
      
      #Box plot for each variable of interest
      tmp_x_lab <- table(eval(parse(text=paste("data_plot$",var_of_interest[i], sep = ""))))
      g <- ggplot(data_plot, aes(y = predictor, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
        geom_boxplot() +
        geom_point() +
        xlab(var_of_interest[i]) +
        ylab(var_predictor[m]) +
        scale_x_discrete(labels=paste(names(tmp_x_lab), "=",tmp_x_lab)) +
        stat_pvalue_manual(to_plot_stat, label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_axis_max) +
        theme_classic()
      p <- g 
      list_plot[[count]] <- p
      names(list_plot)[[count]] <- paste(var_predictor[m], var_of_interest[i], sep = "_")
      count <- count + 1
      fig <- ggplotly(p)
      ggsave(p, filename = paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], ".png", sep = ""), height = 8, width = 12)
      htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "Box_plot_", var_predictor[m], "_", var_of_interest[i], ".html", sep = ""))
      
    }#end Loop go through all var of interest
  }#end loop go through all predictors 
  saveRDS(list_plot, paste(path_output, "enrich_fraction_list_plot.rds", sep = ""))
  
  
  ##Gene set analysis
  data_enrich <- enrich_mutation_data[which(enrich_mutation_data$qvalue_enrich[which(enrich_mutation_data$enrich == "Yes")] <= 0.05),]
  data_depleted <- enrich_mutation_data[which(enrich_mutation_data$qvalue_enrich[which(enrich_mutation_data$enrich == "No")] <= 0.05),]
  gene_enrich_to_keep <- unique(data_enrich$gene_symbol)[-which(unique(data_enrich$gene_symbol) %in% unique(data_depleted$gene_symbol))]
  gene_depleted_to_keep <- unique(data_depleted$gene_symbol)[-which(unique(data_depleted$gene_symbol) %in% unique(data_enrich$gene_symbol))]
  data_enrich <- data_enrich[which(data_enrich$gene_symbol %in% gene_enrich_to_keep),]
  data_depleted <- data_depleted[which(data_depleted$gene_symbol %in% gene_depleted_to_keep),]
  data_combined <- rbind(data_enrich, data_depleted)
  
  res_enrich <- over_representation_analysis(geneList = unique(data_enrich$gene_symbol),
                                             go_of_interest = c("BP", "MF", "CC"),
                                             max_GO_similarity_accepted = 0.5,
                                             cutoff = 0.05,
                                             max_gene_per_item = 20,
                                             max_item_plot = 20,
                                             coocurrence = TRUE,
                                             max_coocurrence = 100,
                                             path_to_store_results = paste(path_output, "gene_set_analysis_enrich", sep = "/"))
  
  res_depleted <- over_representation_analysis(geneList = unique(data_depleted$gene_symbol),
                                               go_of_interest = c("BP", "MF", "CC"),
                                               max_GO_similarity_accepted = 0.5,
                                               cutoff = 0.05,
                                               max_gene_per_item = 20,
                                               max_item_plot = 20,
                                               coocurrence = TRUE,
                                               max_coocurrence = 100,
                                               path_to_store_results = paste(path_output, "gene_set_analysis_depleted", sep = "/"))
  
  
  #return()
}
