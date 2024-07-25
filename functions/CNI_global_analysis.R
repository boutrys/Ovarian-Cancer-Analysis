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

#' CNI_global_analysis
#'
#' given clinical data, and several cni score computed based on different software AND/OR control groups
#' allow to compare cni score resulting from different choices (i.e. software, control group), and
#' plot these results and computes statistics and annotate clinical data with the different cni score 
#' 
#' @param group_name which group of cni score to be compared (TuPo TuPr) 
#' @param path_input where to retrieve clinical data and the cni score from different software
#' @param path_output where to store results
#' @param path_analysis where to store the new annotated clinical data
#'
#'
CNI_global_analysis <- function(group_name = c("TuPr", "TuPo"),
                                control_group_name = c("San"),
                                path_input = "",
                                path_output = getwd(),
                                path_analysis = getwd()){
  #libraries
  library(readxl)
  library(writexl)
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  library(stats)
  library(ggpubr)
  library(RColorBrewer)
  
  
  #create folder for results And load clinical data
  dir.create(path_output)
  clinical_data <- read_xlsx(paste(path_input, "clinical_data_used_for_analysis.xlsx", sep = ""))
  
  #load cni_score and combined them to annotate clinical data
  folder_program <- list.files(path_input, pattern = "score")
  myList <- vector("list", 6*length(control_group_name))
  count <- 1
  max_row <- 0
  for (program in folder_program) {
    folder_control_group <- list.files(paste(path_input, program, "/", sep = ""))
    
    for (control_group in folder_control_group) {
      myList[[count]] <- readRDS(paste(path_input, program, "/", control_group, "/Asc/cni_score.rds", sep = ""))
      if(length(grep(myList[[count]]$name, pattern = "-")) > 0){
        for (loop_name in 1:length(myList[[count]]$name)) {
          myList[[count]]$name[loop_name] <- strsplit(myList[[count]]$name[loop_name], split = "-")[[1]][2]
        }
      }
      if(max_row < dim(myList[[count]])[1]){
        max_row <- dim(myList[[count]])[1]
      }
      names(myList)[count] <- paste(program, control_group, "Asc", sep = "_")
      count <- count + 1
      myList[[count]] <- readRDS(paste(path_input, program, "/", control_group, "/TuPr/cni_score.rds", sep = ""))
      if(length(grep(myList[[count]]$name, pattern = "-")) > 0){
        for (loop_name in 1:length(myList[[count]]$name)) {
          myList[[count]]$name[loop_name] <- strsplit(myList[[count]]$name[loop_name], split = "-")[[1]][2]
        }
      }
      if(max_row < dim(myList[[count]])[1]){
        max_row <- dim(myList[[count]])[1]
      }
      names(myList)[count] <- paste(program, control_group, "TuPr", sep = "_")
      count <- count + 1
      myList[[count]] <- readRDS(paste(path_input, program, "/", control_group, "/TuPo/cni_score.rds", sep = ""))
      if(length(grep(myList[[count]]$name, pattern = "-")) > 0){
        for (loop_name in 1:length(myList[[count]]$name)) {
          myList[[count]]$name[loop_name] <- strsplit(myList[[count]]$name[loop_name], split = "-")[[1]][2]
        }
      }
      if(max_row < dim(myList[[count]])[1]){
        max_row <- dim(myList[[count]])[1]
      }
      names(myList)[count] <- paste(program, control_group, "TuPo", sep = "_")
      count <- count + 1
    }
  }
  
  cni_summary <- data.frame(matrix(NA, ncol = length(myList) + 1, nrow = length(clinical_data$Patient)))
  colnames(cni_summary)[1] <- "Patient"
  cni_summary$Patient <- clinical_data$Patient
  for (i in 1:length(myList)) {
    tmp_idx <- i + 1
    tmp_row <- which(cni_summary$Patient %in% myList[[i]]$name)
    cni_summary[tmp_row,tmp_idx] <- myList[[i]]$cni_score
    colnames(cni_summary)[tmp_idx] <- names(myList)[i]
  }
  
  #Annotate clinical data with cni scores
  idx_not_in_cni <- which(!clinical_data$Patient %in% cni_summary$Patient) 
  no_cni_data <- data.frame(matrix(NA, ncol = dim(cni_summary)[2], nrow = length(idx_not_in_cni)))
  colnames(no_cni_data) <- colnames(cni_summary)
  no_cni_data$Patient <- clinical_data$Patient[idx_not_in_cni]
  cni_to_plot <- cni_summary
  cni_summary <- rbind(cni_summary, no_cni_data)
  cni_summary <- cni_summary[match(clinical_data$Patient, cni_summary$Patient),]
  clinical_data <- cbind(clinical_data, cni_summary[,-1])
  write_xlsx(clinical_data, paste(path_analysis, "annotated_clinical_data.xlsx", sep = ""))
  
  
  ####################   Comparing software and/or controls group used + plot results
  #data pretreatment
  cni_to_plot <- pivot_longer(cni_to_plot, cols = 2:length(colnames(cni_to_plot)), names_to = "Analysis", values_to = "cni_score")
  cni_to_plot$program <- NA
  cni_to_plot$control_group <- NA
  cni_to_plot$group <- NA
  cni_to_plot$id <- NA
  for (i in 1:dim(cni_to_plot)[1]) {
    tmp <- strsplit(cni_to_plot$Analysis[i], split = "_")[[1]]
    cni_to_plot$program[i] <- tmp[1]
    cni_to_plot$control_group[i] <- tmp[4]
    cni_to_plot$group[i] <- tmp[5]
    cni_to_plot$id[i] <- paste(cni_to_plot$Patient[i], tmp[5], sep = "_")
  }
  cni_to_plot <- cni_to_plot[which(cni_to_plot$group %in% group_name),]
  cni_to_plot$abs_diff_cni_software <- NA
  cni_to_plot$abs_diff_cni_control_group <- NA
  for (i in 1:dim(cni_to_plot)[1]) {
    tmp_abs <- which(paste(cni_to_plot$control_group, cni_to_plot$id, sep = "_") == paste(cni_to_plot$control_group[i], cni_to_plot$id[i], sep = "_"))
    tmp_abs_group <- which(paste(cni_to_plot$program, cni_to_plot$id, sep = "_") == paste(cni_to_plot$program[i], cni_to_plot$id[i], sep = "_"))
    if(length(tmp_abs) == 2){
      cni_to_plot$abs_diff_cni_software[i] <-  abs(cni_to_plot$cni_score[tmp_abs[1]] - cni_to_plot$cni_score[tmp_abs[2]])
      cni_to_plot$abs_diff_cni_control_group[i] <- abs(cni_to_plot$cni_score[tmp_abs_group[1]] - cni_to_plot$cni_score[tmp_abs_group[2]])
    }else{
      stop("Not exaclty two group to compare in global cni analysis for absolute difference computation!")
    }
  }
  name_group <- unique(cni_to_plot$group)
  nbr_group <- length(name_group)
  #Add clinical data for plots
  cni_to_plot$Anapath_response <- NA
  for (i in 1:dim(cni_to_plot)[1]) {
    cni_to_plot$Anapath_response[i] <- clinical_data$Anapath_response[which(clinical_data$Patient == cni_to_plot$Patient[i])]
  }
  cni_to_plot$Anapath_response <- ordered(as.factor(cni_to_plot$Anapath_response), c("Bad", "Intermediate", "Good"))
  if(nbr_group == 2){
    cni_to_plot$group <- ordered(as.factor(cni_to_plot$group), c("TuPr", "TuPo"))
  }else if(nbr_group == 3){
    cni_to_plot$group <- ordered(as.factor(cni_to_plot$group), c("Asc", "TuPr", "TuPo"))
  }
  write_xlsx(cni_to_plot, paste(path_output, "data_plot_cni.xlsx", sep = ""))
  name_group <- unique(cni_to_plot$group)
  
  #Analysis of different software for CNI score computation
  software <- unique(cni_to_plot$program)
  nbr_software <- length(software)
  if(nbr_software > 1){
    #for each group control
    uni_group_control <- unique(cni_to_plot$control_group)
    for (i in 1:length(uni_group_control)) {
      red_cni_to_plot <- cni_to_plot[which(cni_to_plot$control_group == uni_group_control[i]),]
      
      ###STATISTICS 
      to_plot_stat <- data.frame(group = name_group)
      res_stat <- vector("list", nbr_group)
      name_stat_test <- c()
      for (k in 1:nbr_group) {
        tmp_cni <- red_cni_to_plot[which(red_cni_to_plot$group == name_group[k]),]
        #check normality
        normality_table <- data.frame(matrix(NA, ncol = 2, nrow = nbr_software))
        colnames(normality_table) <- c("condition", "pvalue")
        for (j in 1:nbr_software) {
          normality_table$pvalue[j] <- shapiro.test(tmp_cni$cni_score[which(tmp_cni$program == software[j])])$p.value
          normality_table$condition[j] <- software[j]
        }
        normality_to_add_name <- paste("normality", normality_table$condition, sep = "_")
        normality_to_add <- data.frame(matrix(normality_table$pvalue, ncol = length(normality_to_add_name), nrow = 1))
        colnames(normality_to_add) <- normality_to_add_name
        #reject normality assumptions if at least one pvalue < alpha = 0.05
        reject_normality <- length(which(normality_table$pvalue < 0.05)) > 0
        #Check homogeneity of variance
        if(reject_normality){
          #Fligner-Killeen test
          homo_var <- fligner.test(cni_score ~ program, data = tmp_cni)$p.value
          
          #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
          y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
          y_max <- y_max *0.1 + y_max
          res_stat[[k]] <- compare_means(cni_score ~ program, 
                                         data = tmp_cni,
                                         method = "wilcox.test",
                                         paired = TRUE) %>%
            mutate(y.position = y_max)
          
          name_stat_test <- c(name_stat_test, "Wilcoxon Signed-rank test")
        }else{
          #comparing two groups normaly distributed
          if(nbr_software == 2){
            #F-test
            homo_var <- var.test(cni_score ~ program, data = tmp_cni)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is true 
            y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            res_stat[[k]] <- compare_means(cni_score ~ program, 
                                           data = tmp_cni,
                                           method = "t.test",
                                           paired = TRUE) %>%
              mutate(y.position = y_max)
            name_stat_test <- c(name_stat_test, "Paired t-test")
          }else if(nbr_software > 2){
            #Barlett's test
            homo_var <- bartlett.test(cni_score ~ program, data = tmp_cni)$p.value
          }
        }
        res_stat[[k]] <- cbind(res_stat[[k]], normality_to_add, data.frame(homogeneity_variance = homo_var))
        
      }
      tmp_stat <- res_stat[[1]]
      for (k in 2:nbr_group) {
        tmp_stat <- rbind(tmp_stat, res_stat[[k]])
      }
      to_plot_stat <- cbind(to_plot_stat, tmp_stat)
      to_plot_stat$name_test <- name_stat_test
      write_xlsx(to_plot_stat, paste(path_output, "Stat_cni_different_software_with_",
                                     uni_group_control[i], "_as_control_group.xlsx", sep = ""))
      name_stat_test <- unique(name_stat_test)
      if(length(name_stat_test) > 1){
        print("different stat test used")
        name_stat_test <- name_stat_test[1]
      }
      
      #BoxPlot difference inbetween software 
      g <- ggplot(red_cni_to_plot, aes(x = program, y = cni_score, label = Patient)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(colour = Anapath_response), size = 3) +
        scale_colour_brewer(palette="Paired") + 
        stat_pvalue_manual(to_plot_stat, label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_max) +
        theme_classic()
      p <- g + facet_wrap(~group, nrow = 1)
      fig <- ggplotly(p)
      ggsave(p, 
             filename = paste(path_output, "Box_plot_cni_different_software_with_",
                              uni_group_control[i], "_as_control_group.png", sep = ""), 
             height = 8, width = 10)
      htmlwidgets::saveWidget(as_widget(fig), 
                              paste(path_output, "Box_plot_cni_different_software_with_",
                                    uni_group_control[i], "_as_control_group.html", sep = ""))
      
      #Dot plot per patient absolute difference in between software for TuPr and TuPo
      tmp_cni_to_plot <- red_cni_to_plot[which(red_cni_to_plot$program == software[1]),]
      g <- ggplot(tmp_cni_to_plot, aes(y = abs_diff_cni_software, x = Patient)) +
        geom_jitter(aes(colour = group), size = 3) +
        scale_colour_brewer(palette="Dark2") + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle=90))
      ggsave(g, filename = paste(path_output, "Box_plot_abs_diff_cni_different_software_with_",
                                 uni_group_control[i], "_as_control_group.png", sep = ""), 
             height = 8, width = 8)
      
      #Paired plot difference inbetween software
      g <- ggplot(red_cni_to_plot, aes(x = program, y = cni_score, label = Patient)) +
        geom_point(aes(color = Anapath_response)) +
        geom_line(aes(group = Patient, color= Anapath_response)) + 
        scale_colour_brewer(palette="Paired") + 
        theme_classic()
      p <- g + facet_wrap(~group, nrow = 1)
      fig <- ggplotly(p)
      ggsave(p, filename = paste(path_output, "paired_cni_different_software_with_",
                                 uni_group_control[i], "_as_control_group.png", sep = ""),
             height = 8, width = 8)
      htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "paired_cni_different_software_with_",
                                                    uni_group_control[i], "_as_control_group.html", sep = ""))
    }#end for loop going through all control groups
    
  }#end of if at least two software
  
  
  #Analysis of different controls groups for CNI score computation
  control_group_name <- unique(cni_to_plot$control_group)
  nbr_control_group <- length(control_group_name)
  if(nbr_control_group > 1){
    #for each software
    for (i in 1:nbr_software) {
       
      red_cni_to_plot <- cni_to_plot[which(cni_to_plot$program == software[i]),]
      
      ###STATISTICS 
      to_plot_stat <- data.frame(group = name_group, name_test = NA)
      res_stat <- vector("list", nbr_group)
      name_stat_test <- c()
      for (k in 1:nbr_group) {
        tmp_cni <- red_cni_to_plot[which(red_cni_to_plot$group == name_group[k]),]
        #check normality
        normality_table <- data.frame(matrix(NA, ncol = 2, nrow = nbr_control_group))
        colnames(normality_table) <- c("condition", "pvalue")
        for (j in 1:nbr_control_group) {
          normality_table$pvalue[j] <- shapiro.test(tmp_cni$cni_score[which(tmp_cni$control_group == control_group_name[j])])$p.value
          normality_table$condition[j] <- control_group_name[j]
        }
        normality_to_add_name <- paste("normality", normality_table$condition, sep = "_")
        normality_to_add <- data.frame(matrix(normality_table$pvalue, ncol = length(normality_to_add_name), nrow = 1))
        colnames(normality_to_add) <- normality_to_add_name
        #reject normality assumptions if at least one pvalue < alpha = 0.05
        reject_normality <- length(which(normality_table$pvalue < 0.05)) > 0
        #Check homogeneity of variance
        if(reject_normality){
          #Fligner-Killeen test
          homo_var <- fligner.test(cni_score ~ control_group, data = tmp_cni)$p.value
          
          #comparing two groups : predictor is categorical, outcome is continuous, and normality is false
          y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
          y_max <- y_max *0.1 + y_max
          res_stat[[k]] <- compare_means(cni_score ~ control_group, 
                                         data = tmp_cni,
                                         method = "wilcox.test",
                                         paired = TRUE) %>%
            mutate(y.position = y_max)
          
          name_stat_test <- c(name_stat_test, "Wilcoxon Signed-rank test")
        }else{
          #comparing two groups normaly distributed
          if(nbr_software == 2){
            #F-test
            homo_var <- var.test(cni_score ~ control_group, data = tmp_cni)$p.value
            
            #comparing two groups : predictor is categorical, outcome is continuous, and normality is true 
            y_max <- max(tmp_cni$cni_score, na.rm = TRUE)
            y_max <- y_max *0.1 + y_max
            res_stat[[k]] <- compare_means(cni_score ~ control_group, 
                                           data = tmp_cni,
                                           method = "t.test",
                                           paired = TRUE) %>%
              mutate(y.position = y_max)
            name_stat_test <- c(name_stat_test, "Paired t-test")
          }else if(nbr_software > 2){
            #Barlett's test
            homo_var <- bartlett.test(cni_score ~ control_group, data = tmp_cni)$p.value
          }
        }
        res_stat[[k]] <- cbind(res_stat[[k]], normality_to_add, data.frame(homogeneity_variance = homo_var))
      }
      tmp_stat <- res_stat[[1]]
      for (k in 2:nbr_group) {
        tmp_stat <- rbind(tmp_stat, res_stat[[k]])
      }
      to_plot_stat <- cbind(to_plot_stat, tmp_stat)
      to_plot_stat$name_test <- name_stat_test
      write_xlsx(to_plot_stat, paste(path_output, "Stat_cni_different_control_group_with_",
                                     software[i], ".xlsx", sep = ""))
      name_stat_test <- unique(name_stat_test)
      if(length(name_stat_test) > 1){
        print("different stat test used")
        name_stat_test <- name_stat_test[1]
      }
      
      #BoxPlot difference inbetween control groups
      g <- ggplot(red_cni_to_plot, aes(x = control_group, y = cni_score, label = Patient)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(colour = Anapath_response), size = 3) +
        scale_colour_brewer(palette="Paired") + 
        stat_pvalue_manual(to_plot_stat, label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_max) +
        theme_classic()
      p <- g + facet_wrap(~group, nrow = 1)
      fig <- ggplotly(p)
      ggsave(p, 
             filename = paste(path_output, "Box_plot_cni_different_control_group_with_",
                              software[i], ".png", sep = ""), 
             height = 8, width = 10)
      htmlwidgets::saveWidget(as_widget(fig), 
                              paste(path_output, "Box_plot_cni_different_control_group_with_",
                                    software[i], ".html", sep = ""))
      
      #Dot plot per patient absolute difference in between control groups for TuPr and TuPo
      tmp_cni_to_plot <- red_cni_to_plot[which(red_cni_to_plot$control_group == uni_group_control[1]),]
      g <- ggplot(tmp_cni_to_plot, aes(y = abs_diff_cni_control_group, x = Patient)) +
        geom_jitter(aes(colour = group), size = 3) +
        scale_colour_brewer(palette="Dark2") + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle=90))
      ggsave(g, filename = paste(path_output, "Box_plot_abs_diff_cni_different_control_group_with_",
                                 software[i], ".png", sep = ""), 
             height = 8, width = 8)
      
      #Paired plot difference inbetween control groups
      g <- ggplot(red_cni_to_plot, aes(x = control_group, y = cni_score, label = Patient)) +
        geom_point(aes(color = Anapath_response)) +
        geom_line(aes(group = Patient, color= Anapath_response)) + 
        scale_colour_brewer(palette="Paired") + 
        theme_classic()
      p <- g + facet_wrap(~group, nrow = 1)
      fig <- ggplotly(p)
      ggsave(p, filename = paste(path_output, "paired_cni_different_control_group_with_",
                                 software[i], ".png", sep = ""),
             height = 8, width = 8)
      htmlwidgets::saveWidget(as_widget(fig), paste(path_output, "paired_cni_different_control_group_with_",
                                                    software[i], ".html", sep = ""))
    }#end for loop going through all software
    
  }#end of if at least two control groups
  
  
  
  #return()
}