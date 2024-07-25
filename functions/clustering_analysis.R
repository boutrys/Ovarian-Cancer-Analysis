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

#' clustering_analysis
#'
#' given clinical data, and results obtained via Raphael Helaers genomicRegionIntersection java program
#' Performs a clustering analysis of the data
#' 
#' @param path where to find the clinical_data.xlsx file on your computer  
#' @param to_remove name of the column you want to get rid of
#' @param output where to store the results
#'
#' @return clinical data with categorical variable, CA 125 ans PCI pre-post computed and put into categories and irrelevent columns removed
#'

clustering_analysis <- function(input_data_cluster = getwd(),
                                clinical_data = getwd(),
                                nbr_group_cluster = c(1,2,3),
                                output = getwd()){
  ###Load needed libraries
  library(tidyverse)
  library(readxl)
  library(ComplexHeatmap)
  require(circlize)
  library(RColorBrewer)
  
  
  ###Load and pretreatment of data
  data <- read_tsv(input_data_cluster)
  data$sample <- gsub("VDB-", "", data$sample)
  data$sample <- substr(data$sample, 1, 4)
  data$id <- paste(data$chromosome, data$start, data$end, sep = "_")
  data$tcn.em <- log(data$tcn.em + 1)
  
  ###Compute the minimum and maximum indices for each item
  chr_indices <- table(data$chromosome[which(data$sample == data$sample[1])])
  chr_labels <- paste("Chr", names(chr_indices))
  tmp_chr <- data$chromosome[which(data$sample == data$sample[1])]
  min_indices <- tapply(seq_along(tmp_chr), tmp_chr, min)
  max_indices <- tapply(seq_along(tmp_chr), tmp_chr, max)
  ### Compute the difference between maximum and minimum indices
  chr_indices <- as.integer(min_indices + (max_indices - min_indices) / 2)
  ###Change label of X chromosome
  chr_labels[which(chr_labels == "Chr 23")] <- "Chr X"
  
  
  ###Reshape data for complexHeatmap
  data <- dplyr::select(data, "sample", "id", "tcn.em")
  all_sample <- unique(data$sample)
  nbr_patient <- length(all_sample)
  all_segment <- unique(data$id)
  nbr_segment <- length(all_segment)
  to_plot <- matrix(0, ncol = nbr_patient, nrow = nbr_segment)
  colnames(to_plot) <- all_sample
  rownames(to_plot) <- all_segment
  for (i in 1:nbr_patient) {
    tmp_data <- data[which(data$sample == all_sample[i]),]
    to_plot[which(all_segment %in% tmp_data$id),i] <- tmp_data$tcn.em
  }
  
  ###Load clinical data to annotate complexheatmap
  clinical_data <- read_xlsx(clinical_data)
  
  annotation_clinical <- data.frame(sample = all_sample, DFS = NA, Anapath_response = NA)
  for (i in 1:length(all_sample)) {
    annotation_clinical$DFS[i] <- clinical_data$DFS[which(clinical_data$Patient == annotation_clinical$sample[i])]
    annotation_clinical$Anapath_response[i] <- clinical_data$Anapath_response[which(clinical_data$Patient == annotation_clinical$sample[i])]
  }
  annotation_clinical$Anapath_response <- factor(annotation_clinical$Anapath_response,
                                                 levels = c("Good", "Intermediate", "Bad"))
  
  
  ###Loop for each number of cluster asked
  for (loop in 1:length(nbr_group_cluster)) {
    ###Heatmap
    ### Create a function to define the color for DFS values
    dfs_color <- function(value) {
      if (value <= 12) {
        return("red")
      } else {
        return("blue")
      }
    }
    coloring_DFS <- c()
    for (i in 1:length(annotation_clinical$DFS)) {
      coloring_DFS <- c(coloring_DFS, dfs_color(annotation_clinical$DFS[i]))
    }
    
    ha1 = HeatmapAnnotation(DFS = anno_barplot(annotation_clinical$DFS, gp = gpar(fill = coloring_DFS)), 
                            Anapath = annotation_clinical$Anapath_response,
                            col =  list(Anapath = c("Good" = "#BEBEBE", "Intermediate" = "#969696", "Bad" = "#525252"), 
                                        DFS = c("DFS <= 12" = "red", "DFS > 12" = "blue")),
                            annotation_name_side = "right")
    #col_fun <-  colorRamp2(c("<= 12", "> 12"), c("blue", "red"))#c("<= 12" = "red", "> 12" = "blue")
    #at = seq(0, 1, by = 1)
    lgd_list <- list(Legend(labels = c("<= 12", "> 12"), 
                            title = "DFS", legend_gp = gpar(fill = c("red", "blue"))))
    name_of_plot <- paste(output, "Heatmap_Facets_TuPr_with_", nbr_group_cluster[loop], "_group.png", sep = "")
    png(name_of_plot,
        width = 1080, height = 1080)
    
    h1 <- Heatmap(to_plot,
                  heatmap_legend_param = list(
                    at = c(1, 2, 3),
                    title = "Copy Number"),
                  cluster_rows = FALSE,
                  show_row_names = FALSE, 
                  clustering_distance_columns = "maximum",
                  clustering_method_columns = "ward.D2",
                  column_km = nbr_group_cluster[loop],
                  bottom_annotation = ha1,
                  left_annotation =  rowAnnotation(label = anno_mark(at = chr_indices, labels = chr_labels, side = "left")))
    draw(h1, annotation_legend_list = lgd_list)
    dev.off()
  }#end of for loop for each number of cluster asked
  
  
}#end of function