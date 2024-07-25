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

#' CNI_analysis
#'
#' given clinical data, and cni score
#' allow to link cni score to clinical data, plot these results and computes statistics
#' 
#' @param group_name which group of cni score to be compared (TuPo TuPr) 
#' @param software which software has been used to compute such cni score
#' @param path_input where to retrieve the cni score
#' @param path_output where to store results
#' @param clinical_data where to find the clinical data
#' @param var_of_interest which columns names of clinical data to be used to analyse the cni score
#'
#'

CNI_analysis <- function(group_name = c(),
                         software = "WisecondorX",
                         path_input = getwd(),
                         path_output = getwd(),
                         clinical_data = getwd(),
                         var_of_interest = c()){
  #library to load
  library(readxl)
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  library(writexl)
  library(ggpubr)
  library(openxlsx)
  library(car)
  library(RColorBrewer)
  library(gridExtra)
  
  dir.create(path_output)
  clinical_data <- read_xlsx(clinical_data)
  
  
  ###load data
  list_cni_per_group <- vector("list", length(group_name))
  for (i in 1:length(group_name)) {
    list_cni_per_group[[i]] <- readRDS(paste(path_input, group_name[i], "/cni_score.rds", sep = ""))
    names(list_cni_per_group)[i] <- group_name[i]
    list_cni_per_group[[i]]$group <- group_name[i]
    if(i == 1){
      cni_data <- list_cni_per_group[[i]]
    }else{
      cni_data <- rbind(cni_data, list_cni_per_group[[i]])
    }
  }
  
  ###data pre-treatment
  if(software == "QDNASeq"){
    patient_list <- unique(cni_data$name)
    nbr_patient <- length(patient_list)
    cni_data$patient <- NA
    tmp_patient <- NA
    for (i in 1:nbr_patient) {
      cni_data$patient[i] <- strsplit(cni_data$name[i], split = "-")[[1]][2]
    }
  }else if(software == "WisecondorX"){
    cni_data$patient <- cni_data$name
  }
  
  #if a patient have missing data for a group then assign "not_paired"
  count_table <- table(cni_data$patient)
  cni_data$statut <- NA
  for (i in 1:length(count_table)) {
    if(count_table[i] == length(group_name)){
      cni_data$statut[which(cni_data$patient == names(count_table[i]))] <- "paired"
    }else{
      cni_data$statut[which(cni_data$patient == names(count_table[i]))] <- "not_paired"
    }
  }
  
  
  clinical_data$DFS[which(clinical_data$DFS < 13)] <- "=<12"
  clinical_data$DFS[which(clinical_data$DFS > 12)] <- ">12"
  clinical_data$DFS <- as.character(clinical_data$DFS)
  clinical_data$censored_DFS <- as.character(clinical_data$censored_DFS)
  clinical_data$Censored_OS <- as.character(clinical_data$Censored_OS)
  
  
  #Remove cni data of patient with no clinical information available
  id_no_clinical_data <- which(!cni_data$patient %in% clinical_data$Patient)
  if(length(id_no_clinical_data) > 0){
    removed_patient_no_clinical <- unique(cni_data$patient[id_no_clinical_data])
    write_xlsx(data.frame(patient = removed_patient_no_clinical),
               paste(path_output, "patient_removed_because_no_clinical_data.xlsx", sep = ""))
    cni_data <- cni_data[-id_no_clinical_data,]
  }
  
  #combine cni and clinical data
  data_to_add <- data.frame(matrix(NA, ncol = length(var_of_interest), nrow = dim(cni_data)[1]))
  colnames(data_to_add) <- var_of_interest
  new_patient_list <- unique(cni_data$patient)
  for (i in 1:length(new_patient_list)) {
    tmp_col_clinical <- which(colnames(clinical_data) %in% var_of_interest)
    tmp_row_clinical <- which(clinical_data$Patient %in% new_patient_list[i])
    tmp_clinical <- clinical_data[tmp_row_clinical, tmp_col_clinical]
    data_to_add[which(cni_data$patient == new_patient_list[i]),] <- tmp_clinical
  }
  cni_data <- cbind(cni_data,
                    data_to_add)
  
  #Ordering of variables for the plots (can be change here)
  if(length(unique(cni_data$group)) == 2){
    cni_data$group <- ordered(as.factor(cni_data$group), c("TuPr", "TuPo"))
  }else if(length(unique(cni_data$group)) == 3){
    cni_data$group <- ordered(as.factor(cni_data$group), c("Asc", "TuPr", "TuPo"))
  }
  if(length(which(var_of_interest == "BRCA")) > 0){
    cni_data$BRCA <- ordered(as.factor(cni_data$BRCA), c("Unknown", "No", "Yes"))
  }
  if(length(which(var_of_interest == "FIGO_stage")) > 0){
    cni_data$FIGO_stage <- ordered(as.factor(cni_data$FIGO_stage), c("IIIB", "IIIc", "IVA", "IVB"))
  }
  if(length(which(var_of_interest == "CA125_post")) > 0){
    cni_data$CA125_post <- ordered(as.factor(cni_data$CA125_post), c("Bad", "Intermediate", "Good"))
  }
  if(length(which(var_of_interest == "Diff_PCI")) > 0){
    cni_data$Diff_PCI <- ordered(as.factor(cni_data$Diff_PCI), c("Bad", "Intermediate", "Good"))
  } 
  if(length(which(var_of_interest == "amount_RD")) > 0){
    cni_data$amount_RD <- ordered(as.factor(cni_data$amount_RD), c("NO residual", "Milimeter", "Supra milimeter", "+5cm"))
  }
  if(length(which(var_of_interest == "Anapath_response")) > 0){
    cni_data$Anapath_response <- ordered(as.factor(cni_data$Anapath_response), c("Bad", "Intermediate", "Good"))
  }
  if(length(which(var_of_interest == "censored_DFS")) > 0){
    cni_data$censored_DFS <- ordered(as.factor(cni_data$censored_DFS), c("0", "1"))
  }
  if(length(which(var_of_interest == "Censored_OS")) > 0){
    cni_data$Censored_OS <- ordered(as.factor(cni_data$Censored_OS), c("0", "1"))
  }
  if(length(which(var_of_interest == "DFS")) > 0){
    cni_data$DFS <- ordered(as.factor(cni_data$DFS), c("=<12", ">12"))
  }
  if(length(which(var_of_interest == "PFI")) > 0){
    cni_data$PFI <- ordered(as.factor(cni_data$PFI), c("platinum resistant", "semi sensitive", "platinum sensitive"))
  }
  
  #Loop go through all var of interest
  nbr_group <- length(group_name)
  for (i in 1:length(var_of_interest)) {
    
    #paired data plots
    paired_cni_data <- cni_data[which(cni_data$statut == "paired"),]
    tmp_x_lab <- table(paired_cni_data$group)
    #var of interest 1 must be Anapath_response
    if(i == 1 || i == 3){
      custom_palette <- c("#FB9A99", "#FFD92F", "#33A02C")
      #scale_colour_manual(values = custom_palette)
    }
    if(i == 2){
      custom_palette <- c("#FB9A99", "#33A02C")
    }
    g <- ggplot(paired_cni_data, aes(x = group, y = cni_score, label = patient)) +
      geom_point(aes(color=eval(parse(text=paste(var_of_interest[i]))))) +
      geom_line(aes(group = patient, color=eval(parse(text=paste(var_of_interest[i]))))) + 
      scale_colour_manual(values = custom_palette) + 
      labs(color = var_of_interest[i]) +
      scale_x_discrete(labels=paste(names(tmp_x_lab), "=",tmp_x_lab)) +
      theme_classic()
    fig <- ggplotly(g)
    ggsave(g, filename = paste(path_output, "paired_cni_", var_of_interest[i], ".png", sep = ""), height = 8, width = 8)
    htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "paired_cni_", var_of_interest[i], ".html", sep = ""))
    
    #not paired
    not_paired_cni_data <- cni_data[which(cni_data$statut == "not_paired"),]
    
    
    ###STATISTICS 
    value_var_of_interest <- unique(eval(parse(text=paste("cni_data$",var_of_interest[i], sep = ""))))
    value_var_of_interest <- value_var_of_interest[order(value_var_of_interest)]
    nbr_value_var_of_interest <- length(value_var_of_interest)
    res_stat <- vector("list", nbr_group)
    name_stat_test <- c()
    for (k in 1:nbr_group) {
      tmp_cni <- cni_data[which(cni_data$group == group_name[k]),-which(colnames(cni_data) %in% var_of_interest[-i])]
      colnames(tmp_cni)[which(colnames(tmp_cni) == var_of_interest[i])] <- "variable"
      #check normality
      normality_table <- data.frame(matrix(NA, ncol = 2, nrow = nbr_value_var_of_interest))
      colnames(normality_table) <- c("condition", "pvalue")
      for (j in 1:nbr_value_var_of_interest) {
        err <- try(tmp <- shapiro.test(tmp_cni$cni_score[which(tmp_cni$variable == value_var_of_interest[j])])$p.value,
                   silent = TRUE)
        if(!is.character(err)){
          normality_table$pvalue[j] <- err
        }
        normality_table$condition[j] <- as.character(value_var_of_interest[j]) 
      }
      normality_to_add_name <- paste("normality", normality_table$condition, sep = "_")
      normality_to_add <- data.frame(matrix(normality_table$pvalue, ncol = length(normality_to_add_name), nrow = 1))
      colnames(normality_to_add) <- normality_to_add_name
      #reject normality assumptions if at least one pvalue < alpha = 0.05
      reject_normality <- length(which(normality_table$pvalue < 0.05)) > 0
      #Check homogeneity of variance
      if(reject_normality){
        if(nbr_value_var_of_interest == 2){
          #Fligner-Killeen test
          homo_var <- fligner.test(cni_score ~ variable, data = tmp_cni)$p.value
          
          #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
          y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
          y_max <- y_max *0.1 + y_max
          tmp_res_stat <- compare_means(cni_score ~ variable, 
                                        data = tmp_cni,
                                        method = "wilcox.test",
                                        paired = FALSE) %>%
            mutate(y.position = y_max)
          res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
          name_stat_test <- c(name_stat_test, "Wilcoxon Rank-Sum test")
        }else if(nbr_value_var_of_interest > 2){
          #Levene’s test
          homo_var <- leveneTest(cni_score ~ variable, data = tmp_cni)$`Pr(>F)`[1]
          
          #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
          y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
          y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
          tmp_res_stat <- compare_means(cni_score ~ variable, 
                                        data = tmp_cni,
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
          homo_var <- var.test(cni_score ~ variable, data = tmp_cni)$p.value
          
          #comparing two groups : predictor is categorical, outcome is continuous, and normality is true 
          y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
          y_max <- y_max *0.1 + y_max
          tmp_res_stat <- compare_means(cni_score ~ variable, 
                                        data = tmp_cni,
                                        method = "t.test",
                                        paired = FALSE) %>%
            mutate(y.position = y_max)
          res_stat[[k]] <- cbind(data.frame(group = group_name[k]), tmp_res_stat)
          name_stat_test <- c(name_stat_test, "T-test")
        }else if(nbr_value_var_of_interest > 2){
          #Barlett's test
          homo_var <- bartlett.test(cni_score ~ variable, data = tmp_cni)$p.value
          
          #TODO
          y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
          y_max <- y_max *seq(from = 0.1, by = 0.2, length.out = nbr_value_var_of_interest) + y_max
          tmp_res_stat <- compare_means(cni_score ~ variable, 
                                        data = tmp_cni,
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
      res_stat[[k]] <- cbind(res_stat[[k]], normality_to_add, data.frame(homogeneity_variance = homo_var))
      
    }
    to_plot_stat <- res_stat[[1]]
    for (k in 2:nbr_group) {
      to_plot_stat <- rbind(to_plot_stat, res_stat[[k]])
    }
    to_plot_stat$name_test <- name_stat_test
    write_xlsx(to_plot_stat, paste(path_output, "Stat_cni_", var_of_interest[i], ".xlsx", sep = ""))
    to_plot_stat$p.adj <- round(to_plot_stat$p.adj, digits = 4)
    name_stat_test <- unique(name_stat_test)
    if(length(name_stat_test) > 1){
      print("different stat test used")
      name_stat_test <- ""
    }
    to_plot_stat$group <- ordered(as.factor(to_plot_stat$group), c("TuPr", "TuPo"))
    y_max <- max(y_max, na.rm = TRUE)
    
    
    
    #Box plot for each variable of interest
    tmp_x_lab <- table(eval(parse(text=paste("cni_data$",var_of_interest[i], sep = ""))))
    g <- ggplot(cni_data, aes(y = cni_score, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 3) +
      scale_colour_brewer(palette="Paired") + 
      xlab(var_of_interest[i]) +
      stat_pvalue_manual(to_plot_stat, label = paste(name_stat_test, "pvalue = {p.adj}")) +
      ylim(0,y_max) +
      scale_x_discrete(labels=paste(names(tmp_x_lab), "=",tmp_x_lab)) +
      theme_classic()
    p <- g + facet_wrap(~group, nrow = 1)
    #fig <- ggplotly(p)
    #ggsave(p, filename = paste(path_output, "Box_plot_cni_", var_of_interest[i], ".png", sep = ""), height = 8, width = 12)
    htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "Box_plot_cni_", var_of_interest[i], ".html", sep = ""))
    
    
    cni_data_tupr <- cni_data[which(cni_data$group == "TuPr"),]
    tmp_x_lab_tupr <- table(eval(parse(text=paste("cni_data_tupr$",var_of_interest[i], sep = ""))))
    cni_data_tupo <- cni_data[which(cni_data$group == "TuPo"),]
    tmp_x_lab_tupo <- table(eval(parse(text=paste("cni_data_tupo$",var_of_interest[i], sep = ""))))
    plot1 <- ggplot(cni_data_tupr, aes(y = cni_score, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 3) +
      scale_colour_brewer(palette="Paired") + 
      xlab(var_of_interest[i]) +
      stat_pvalue_manual(to_plot_stat[which(to_plot_stat$group == "TuPr"),], label = paste(name_stat_test, "pvalue = {p.adj}")) +
      ylim(0,y_max) +
      scale_x_discrete(labels=paste(names(tmp_x_lab_tupr), "=",tmp_x_lab_tupr)) +
      labs("TuPr") +
      theme_classic()
    plot2 <- ggplot(cni_data_tupo, aes(y = cni_score, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 3) +
      scale_colour_brewer(palette="Paired") + 
      xlab(var_of_interest[i]) +
      stat_pvalue_manual(to_plot_stat[which(to_plot_stat$group == "TuPo"),], label = paste(name_stat_test, "pvalue = {p.adj}")) +
      ylim(0,y_max) +
      scale_x_discrete(labels=paste(names(tmp_x_lab_tupo), "=",tmp_x_lab_tupo)) +
      labs("TuPo") +
      theme_classic()
    new <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)
    ggsave(new, filename = paste(path_output, "Box_plot_cni_", var_of_interest[i], ".png", sep = ""), height = 8, width = 12)
    
    
  }#end Loop go through all var of interest
  
  
  
  ###Box plots to compare TuPr VS TuPO for each level of a var of interest
  var_predictor <- c("cni_score")
  list_plot <- vector(mode = "list", length = length(var_predictor)*length(var_of_interest))
  count <- 1
  cni_data_paired <- cni_data[-which(cni_data$statut == "not_paired"),]
  y_axis_max <- max(cni_data_paired$cni_score, na.rm = TRUE)
  y_axis_max <- y_axis_max * 0.1 + y_axis_max
  for (m in 1:length(var_predictor)) {
    data_plot <- cni_data_paired
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
        ylab(paste("Cni score")) +
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
  
  #save all cni score
  write_xlsx(cni_data, paste(path_output, software, "_cni_score_summary.xlsx", sep = ""))
  
  #return()
}