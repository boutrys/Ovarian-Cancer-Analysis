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

#  Calculation of the CNI score (copy number instability score) using shallow NGS data : 
#  quantitative assessment of genome-wide chromosomal instability 
#  Natasha HONORE, Simon BOUTRY and C?dric VAN MARCKE (UCLouvain, Brussels, Belgium)
#  based on doi: 10.3390/cancers14010168 
CNI_score_computation <- function(path_input = getwd(),
                                  software = c("QDNASeq", "WisecondorX"),
                                  control_name = c(),
                                  group_name =c(),
                                  path_output = getwd()){
  #library to load
  library(tidyverse)
  library(zoo)
  library(writexl)
  
  
  ### functions to process data
  process_data <- function(path = getwd(),
                           max_val = 20, 
                           min_val = -20,
                           control = TRUE,
                           output_path = getwd(),
                           germline_data,
                           wisecondor = FALSE,
                           group){
    if(wisecondor){
      #Hypothese in wisecondor : in folder controls, all sample without "San" are considered to be part of groupe controls
      if(wisecondor && group == "controls"){
        to_process <- list.files(path = path, 
                                 pattern = "San")
        to_process_2 <- list.files(path = path)
        to_process_2 <- to_process_2[-which(to_process_2 %in% to_process)]
        to_process <- to_process_2
      }else{
        to_process <- list.files(path = path, 
                                 pattern = group) 
      }
    }else{
      to_process <- list.files(path = path, 
                               pattern = "igv")
    }
    
    
    general_data <- data.frame(chr=numeric(),
                               start=numeric(),
                               end=numeric(),
                               feat=character(),
                               sample=numeric(),
                               name=character(),
                               Log2count=numeric(),
                               windows_10=numeric(),
                               mean_log2=numeric(),
                               stringsAsFactors=FALSE)
    count_under <- 0
    count_above <- 0
    for (filename in to_process) {
      if(wisecondor){
        to_read <- list.files(path = paste(path, filename, sep = ""), 
                              pattern = "bins.bed")
        shallow_data_g <- read_tsv(paste(path, filename, "/", to_read, sep = ""))
        #remove NA
        shallow_data_g <- shallow_data_g[which(!is.na(shallow_data_g$ratio)),]
      }else{
        shallow_data_g <- read_tsv(paste(path, filename, sep = ""), skip = 2)
      }
      shallow_data_g <- shallow_data_g[shallow_data_g[,5]>min_val,] #remove the extreme outliers = seq error/artifact
      count_under <- count_under + length(which(!shallow_data_g[,5]>min_val))
      shallow_data_g <- shallow_data_g[shallow_data_g[,5]<max_val,] #remove the extreme outliers = seq error/artifact
      count_above <- count_above + length(which(!shallow_data_g[,5]<max_val))
      if(wisecondor){
        name <- strsplit(to_read, split = "-")[[1]][2]
        colnames(shallow_data_g) <- c("chr", "start", "end", "feat", "sample", "wisecondor_zscore")
      }else{
        name <- colnames(shallow_data_g[,5])
        colnames(shallow_data_g) <- c("chr", "start", "end", "feat", "sample")
      }
      
      
      shallow_data_g$name <- name
      shallow_data_g['sample'][shallow_data_g['sample'] == 0] <- 0.001 #to avoid infinite log values
      
      shallow_data_g$Log2count <- log2(abs(shallow_data_g$sample))
      
      shallow_data_g_10 <- shallow_data_g %>%
        group_by(name) %>%
        mutate(mean_log2 = zoo::rollmean(Log2count, k = 10, fill = NA)) %>% # calculate means of sliding windows of 10 segments
        ungroup()
      
      general_data <- rbind(general_data, shallow_data_g_10)
    }
    
    #plot the data
    new_min <- min(general_data$mean_log2, na.rm = TRUE) - 1
    new_max <- max(general_data$mean_log2, na.rm = TRUE) + 1
    g <- ggplot(general_data, aes(mean_log2)) +
      geom_histogram(breaks = seq(new_min, new_max, 0.05)) 
    ggsave(g, filename = paste(output_path, "histogram_MeanLog2.png", sep = ""), height = 8, width = 8)
    g <- g + facet_wrap(~name)
    ggsave(g, filename = paste(output_path, "histogram_meanLog2_per_sample.png", sep = ""), height = 10, width = 10)
    
    
    #scaled and center data
    general_data$mean_log2_scaled <- as.numeric(scale(general_data$mean_log2, center = TRUE, scale = TRUE)) # center and scale the values
    general_data <- na.omit(general_data)
    
    new_min <- min(general_data$mean_log2_scaled, na.rm = TRUE) - 1
    new_max <- max(general_data$mean_log2_scaled, na.rm = TRUE) + 1
    g <- ggplot(general_data, aes(mean_log2_scaled)) +
      geom_histogram(breaks = seq(new_min, new_max, 0.05)) 
    ggsave(g, filename = paste(output_path, "histogram_ScaledMeanLog2.png", sep = ""), height = 8, width = 8)
    g <- g + facet_wrap(~name)
    ggsave(g, filename = paste(output_path, "histogram_ScaledmeanLog2_per_sample.png", sep = ""), height = 10, width = 10)
    
    
    if(control){
      general_data_by_feat = general_data %>% group_by(feat)  %>% # those means and SD per segment of chromosomes will serve to calculate the Z scores for each equivalent segment of each sample
        summarise(mean_log2_scaled_by_feat = mean(mean_log2_scaled),
                  sd_log2_scaled_by_feat = sd(mean_log2_scaled),
                  .groups = 'drop')
      
      saveRDS(general_data_by_feat, paste(output_path, "general_data_by_feat.rds", sep = ""))
      return(list(general_data_by_feat, count_under, count_above))
    }else{
      general_data$Z_score <- NA
      
      print_every_n_iteration <- 1000
      for (i in 1:dim(germline_data)[1]) { #calculate the Z score for each segment, compared to the mean and SD of that specific segment in your panel of normals
        germline_general_data_temp <- germline_data[i,]
        id <- which(general_data$feat == germline_general_data_temp$feat)
        general_data$Z_score[id] <- ((general_data$mean_log2_scaled[id] - germline_general_data_temp$mean_log2_scaled_by_feat) / germline_general_data_temp$sd_log2_scaled_by_feat)
        if (floor(i/print_every_n_iteration) * print_every_n_iteration == i) {
          msg <- sprintf("%d/%d", i, dim(germline_data)[1])
          print(msg)
        }
      }
      
      g <- ggplot(general_data, aes(Z_score)) +
        geom_histogram()
      ggsave(g, filename = paste(output_path, "histogram_Zscore.png", sep = ""), height = 8, width = 8)
      #hist(general_data$Z_score, breaks = 100)
      
      general_data_sign <- general_data %>% 
        filter(abs(Z_score) > 2.84) # keep only the Z scores that are considered extreme
      
      cni_score <- general_data_sign %>% group_by(name) %>% summarise(cni_score= sum(abs(Z_score)))
      saveRDS(cni_score, paste(output_path, "cni_score.rds", sep = ""))
      saveRDS(general_data, paste(output_path, "general_data.rds", sep = ""))
      return(list(general_data, cni_score, count_under, count_above))
    }
    
  }#end of function
  
  
  #Go through each controls groups
  for (i in 1:length(control_name)) {
    # Set the names of the files to process and where to store results
    path_output_analysis <- paste(path_output, control_name[i],  "/", sep = "")
    dir.create(path_output_analysis)
    path_output_control <- paste(path_output_analysis, control_name[i], "/", sep = "")
    dir.create(path_output_control)
    if(software == "QDNASeq"){
      path_control <- paste(path_input, control_name[i], "/", sep = "")
      resp <- process_data(path = path_control,
                           output_path = path_output_control)
      germline_control <- resp[[1]]
    }else if (software == "WisecondorX"){
      path_control <- paste(path_input, "controls/", sep = "")
      resp <- process_data(path = path_control,
                           output_path = path_output_control,
                           wisecondor = TRUE,
                           group = control_name[i])
      germline_control <- resp[[1]]
    }else{
      print("softwar must be QDNASeq or WisecondorX")
    }
    
    #Go through each groups of interest
    for (j in 1:length(group_name)) {
      path_output_group <- paste(path_output_analysis, group_name[j], "/", sep = "")
      dir.create(path_output_group)
      if(software == "QDNASeq"){
        path_input_group <- paste(path_input, group_name[j], "/", sep = "")
        resp <- process_data(path = path_input_group,
                             control = FALSE,
                             output_path = path_output_group,
                             germline_data = germline_control)
        data_group <- resp[[1]]
        cni_group <- resp[[2]]
      }else if (software == "WisecondorX"){
        path_patient <- paste(path_input, "patients/", sep = "")
        resp <- process_data(path = path_patient,
                             control = FALSE,
                             output_path = path_output_group,
                             germline_data = germline_control,
                             wisecondor = TRUE,
                             group = group_name[j])
        data_group <- resp[[1]]
        cni_group <- resp[[2]]
      }else{
        print("softwar must be QDNASeq or WisecondorX")
      }
    }#end loop go through all groups
    
  }#end loop go through all control groups
  
  #return()
}