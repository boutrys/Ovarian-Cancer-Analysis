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


qdnaseq_analysis <- function(path_input = getwd(),
                             path_output = getwd(),
                             do_segmentation = FALSE){
  #library to load
  library(writexl)
  
  
  #Qdnaseq .igv and .seg analysis
  group <- list.files(path_input)
  nbr_group <- length(group)
  data_list <- vector("list", nbr_group)
  name_col <- c("sample", "group", "average", "sd")
  data_summary <- data.frame(matrix(NA, ncol = length(name_col), nrow = 1000))
  colnames(data_summary) <- name_col
  count <- 1
  
  #for each group
  for (i in 1:nbr_group) {
    path_group <- paste(path_input, group[i], "/", sep = "")
    samples_files <- list.files(path_group)
    nbr_samples <- length(samples_files)
    data_list[[i]] <- vector("list", nbr_samples)
    
    #for each sample within a group
    for (j in 1:nbr_samples) {
      data <- read.table(paste(path_group, samples_files[j], sep = ""), header = TRUE)
      data_list[[i]][[j]] <- data
      
      if(i != 2){
        tmp <- strsplit(colnames(data)[5], split = ".", fixed = TRUE)[[1]]
        tmp_name <- tmp[2]
        tmp_group <- tmp[3]
      }else{
        tmp_name <- colnames(data)[5]
        tmp_group <- "control"
      } 
      tmp_average <-  sum(data[,5])/dim(data)[1]
      tmp_sd <- sd(data[,5])
      data_summary[count,] <- c(tmp_name, tmp_group, tmp_average, tmp_sd)
      count <- count + 1
      
      nbr_dif_feature <- length(which(!data_list[[i]][[j]]$feature == data_list[[1]][[1]]$feature))
      #check if same features
      if(nbr_dif_feature > 0){
        print(colnames(data_list[[i]][[j]]))
      }#end of if different features
      
    }#End for each sample within a group
    
    print(i)
  }#End for each group
  
  data_summary <- data_summary[-which(is.na(data_summary$sample)),]
  write_xlsx(data_summary, paste(path_output, "summary_bin.xlsx", sep = ""))
  unique_sample <- unique(data_summary$sample[which(duplicated(data_summary$sample))])
  nbr_unique_sample <- length(unique_sample)
  comparisson_per_sample <- data.frame(matrix(NA, ncol = 6, nrow = nbr_unique_sample))
  colnames(comparisson_per_sample) <- c("sample", "Asc", "San", "TuPo", "TuPr", "evol_TuPr_TuPo")
  comparisson_per_sample$sample <- unique_sample
  for (i in 1:nbr_unique_sample) {
    tmp <- data_summary[which(data_summary$sample == unique_sample[i]),]
    tmp_idx <- which(colnames(comparisson_per_sample) %in% tmp$group)
    if(length(tmp_idx) > 0){
      comparisson_per_sample[i,tmp_idx] <- tmp$average 
      comparisson_per_sample[i,6] <- as.double(tmp$average[3]) - as.double(tmp$average[4])
    }else{
      stop("problem in qdnaseq analaysis : group name do not match")
    }
  }
  write_xlsx(comparisson_per_sample, paste(path_output, "comparisson_per_sample.xlsx", sep = ""))
  
  
  
  ############################################################ SEGMENTATION ######################################################################
  if(do_segmentation){
    path <- paste(path_input, "segmentation/", sep = "")
    name_file <- list.files(path)
    nbr_file <- length(name_file)
    count <- 1
    data_list <- vector("list", nbr_file/2)
    tmp_count <- 1
    data_summary <- data.frame(matrix(NA, ncol = 7, nrow = nbr_file/2))
    colnames(data_summary) <- c("sample", "Nbr_seg_Po", "Nbr_seg_Pr", "Diff_nbr_seg", "Average_LOG2_ratio_mean_Po", "Average_LOG2_ratio_mean_Pr", "Diff_average")
    data_cnv <- data.frame()
    for (i in 1:nbr_file) {
      if(count == 1){
        data_list[[tmp_count]] <- vector("list", 2)
        count <- 2
        data <- read.table(paste(path, name_file[i], sep = ""), header = TRUE)
        data$feature <- paste(data$CHROMOSOME, data$START, data$STOP, sep = "_")
        data_list[[tmp_count]][[1]] <- data
        data_summary$sample[tmp_count] <- strsplit(data$SAMPLE_NAME[1], split = "-")[[1]][2]
        data_summary$Nbr_seg_Po[tmp_count] <- dim(data)[1]
        data_summary$Average_LOG2_ratio_mean_Po[tmp_count] <- sum(data$LOG2_RATIO_MEAN) / dim(data)[1]
      }else{
        data <- read.table(paste(path, name_file[i], sep = ""), header = TRUE)
        data$feature <- paste(data$CHROMOSOME, data$START, data$STOP, sep = "_")
        data_list[[tmp_count]][[2]] <- data
        count <- 1
        data_summary$Nbr_seg_Pr[tmp_count] <- dim(data)[1]
        data_summary$Diff_nbr_seg[tmp_count] <- data_summary$Nbr_seg_Po[tmp_count] - data_summary$Nbr_seg_Pr[tmp_count]
        data_summary$Average_LOG2_ratio_mean_Pr[tmp_count] <-  sum(data$LOG2_RATIO_MEAN) / dim(data)[1]
        data_summary$Diff_average[tmp_count] <- data_summary$Average_LOG2_ratio_mean_Po[tmp_count] - data_summary$Average_LOG2_ratio_mean_Pr[tmp_count]
        tmp_count <- tmp_count + 1
        
      }
      data_cnv <- rbind(data_cnv, data)
    }
    write_xlsx(data_summary, paste(path_output, "summary_segment.xlsx", sep = ""))
    name_dup_cnv <- unique(data_cnv$feature[which(duplicated(data_cnv$feature))])
    nbr_dup_cnv <- length(name_dup_cnv)
    for (i in 1:nbr_dup_cnv) {
      tmp <- data_cnv[which(data_cnv$feature == name_dup_cnv[i]),]
      write_xlsx(tmp, paste(path_output, i,"_dup_cnv.xlsx", sep = ""))
    }
    
  }#end of if do_segmentation
  
  
  #return()
}
