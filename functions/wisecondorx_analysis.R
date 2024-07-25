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

wisecondorx_analysis <- function(path = getwd(),
                                 path_annotsv = "",
                                 path_output = getwd()){
  #library to load
  library(writexl)
  library(tidyverse)
  library(IRanges)
  library(biomaRt)
  
  
  pattern_name_file_cnv <- "aberration" 
  path_patient <- paste(path, "patients/", sep = "")
  path_control <- paste(path, "controls/", sep = "")
  
  
  #Create results repertory
  path_output_plot <- paste(path_output, "plot_from_WisecondroX/", sep = "")
  dir.create(path_output_plot)
  
  ### Combine all cnv data patient
  name_patient_file <- list.files(path_patient)
  nbr_patient_file <- length(name_patient_file)
  for (i in 1:nbr_patient_file) {
    path_tmp_file <- paste(path_patient, name_patient_file[i], "/", sep = "")
    tmp_name_patient <- strsplit(name_patient_file[i], split = "-")[[1]][2]
    tmp_kind <- strsplit(name_patient_file[i], split = "-")[[1]][3]
    tmp_list_file <- list.files(path_tmp_file)
    tmp_plot <- paste(path_tmp_file, tmp_list_file[grep(tmp_list_file, pattern = "plots")], "/", sep = "")
    tmp_path_output_plot <- paste(paste(paste(path_output_plot, tmp_name_patient, sep = ""), tmp_kind, sep = "-"), "/", sep = "")
    dir.create(tmp_path_output_plot)
    file.copy(from = paste0(tmp_plot, list.files(tmp_plot)), 
              to = paste0(tmp_path_output_plot, list.files(tmp_plot)))
    tmp_cnv <- read.table(paste(path_tmp_file, tmp_list_file[grep(tmp_list_file, pattern = pattern_name_file_cnv)], sep = ""), header = TRUE)
    if(dim(tmp_cnv)[1] > 0){
      tmp_cnv$sample <- tmp_name_patient
      tmp_cnv$Kind_sample <- tmp_kind
      if(i ==1){
        #TO DO will crash if first patients has no cnv (file empty)
        cnv <- tmp_cnv
      }else{
        cnv <- rbind(cnv, tmp_cnv)
      }
    }
    
  }
  
  
  ### Combine all cnv data controls
  name_control_file <- list.files(path_control)
  nbr_control_file <- length(name_control_file)
  for (i in 1:nbr_control_file) {
    path_tmp_file <- paste(path_control, name_control_file[i], "/", sep = "")
    tmp_name_control <- name_control_file[i]
    if(length(grep(tmp_name_control, pattern = "VDB")) == 1){
      tmp_name_control <- strsplit(name_control_file[i], split = "-")[[1]][2]
    }
    tmp_kind <- "control"
    tmp_list_file <- list.files(path_tmp_file)
    tmp_plot <- paste(path_tmp_file, tmp_list_file[grep(tmp_list_file, pattern = "plots")], "/", sep = "")
    tmp_path_output_plot <- paste(paste(paste(path_output_plot, tmp_name_control, sep = ""), tmp_kind, sep = "-"), "/", sep = "")
    dir.create(tmp_path_output_plot)
    file.copy(from = paste0(tmp_plot, list.files(tmp_plot)), 
              to = paste0(tmp_path_output_plot, list.files(tmp_plot)))
    tmp_cnv <- read.table(paste(path_tmp_file, tmp_list_file[grep(tmp_list_file, pattern = pattern_name_file_cnv)], sep = ""), header = TRUE)
    if(dim(tmp_cnv)[1] > 0){
      tmp_cnv$sample <- tmp_name_control
      tmp_cnv$Kind_sample <- tmp_kind
      if(!exists("cnv")){
        #TO DO will crash if first patients has no cnv (file empty)
        cnv <- tmp_cnv
      }else{
        cnv <- rbind(cnv, tmp_cnv)
      }
    }
    
  }
  
  
  ### Analyze from cnv point of view
  cnv$id <- paste(cnv$chr, cnv$start, cnv$end, sep = "_")
  unique_cnv_id <- unique(cnv$id)
  nbr_unique_cnv_id <- length(unique(cnv$id))
  
  #Look at duplicated CNV(s)
  dup_id <- which(duplicated(cnv$id))
  nbr_dup_cnv <- length(dup_id)
  unique_dup_id <- unique(cnv$id[dup_id])
  nbr_unique_dup <- length(unique_dup_id)
  
  #annotate number of time a cnv is present and in which data
  cnv$Ocurrence_cnv <- 0
  cnv$List_occurence <- NA
  for (i in 1:nbr_unique_cnv_id) {
    tmp_id <- which(cnv$id == unique_cnv_id[i])
    cnv$Ocurrence_cnv[tmp_id] <- length(tmp_id)
    if(length(tmp_id) > 1){
      cnv$List_occurence[tmp_id] <- paste(paste(cnv$sample[tmp_id], cnv$Kind_sample[tmp_id], sep = "-"), collapse = "/")
    }else{
      cnv$List_occurence[tmp_id] <- paste(cnv$sample[tmp_id], cnv$Kind_sample[tmp_id], sep = "-")
    }
  }
  
  #save combined cnv and annotated results 
  write_xlsx(cnv, paste(path_output, "combined_cnv.xlsx", sep = ""))
  
  
  ###Analyze by sample
  samples <- unique(cnv$sample)
  nbr_sample <- length(samples)
  by_kind <- unique(cnv$Kind_sample)
  name_col <- c("sample",
                "Nbr_cnv",
                "Nbr_unique_cnv",
                paste("Nbr_cnv", by_kind, sep = "_"))
  data_sample <- data.frame(matrix(NA, ncol = length(name_col), nrow = nbr_sample))
  colnames(data_sample) <- name_col
  for (i in 1:nbr_sample) {
    tmp <- cnv[which(cnv$sample == samples[i]),]
    data_sample$sample[i] <- unique(tmp$sample)
    data_sample$Nbr_cnv[i] <- dim(tmp)[1]
    data_sample$Nbr_unique_cnv[i] <- length(unique(tmp$id))
    tmp_kind <- unique(tmp$Kind_sample)
    for (j in 1:length(tmp_kind)) {
      id_to_change <- grep(colnames(data_sample), pattern = tmp_kind[j])
      small_tmp <- tmp[which(tmp$Kind_sample == tmp_kind[j]),]
      data_sample[i,id_to_change[1]] <- dim(small_tmp)[1]
    }
  }
  write_xlsx(data_sample, paste(path_output, "summary_by_sample.xlsx", sep = ""))
  
  
  ###Annotsv
  if(path_annotsv != ""){
    file_annotsv <- list.files(path_annotsv)
    annot_exist <- 0
    for (i in 1:length(file_annotsv)) {
      #retrieve from file, sample name and kind of sample
      if(length(grep(file_annotsv[i], pattern = "VDB")) > 0){
        tmp_name <- strsplit(file_annotsv[i], split = "-")[[1]][2]
        if(length(grep(file_annotsv[i], pattern = "San")) > 0){
          #control sang from patient
          tmp_kind <- "control"
        }else{
          #patient
          tmp_kind <- strsplit(strsplit(file_annotsv[i], split = "-")[[1]][3], split = ".annot")[[1]][1]
        }
      }else{
        #control
        tmp_name <- strsplit(file_annotsv[i], split = ".annotated")[[1]][1]
        tmp_kind <- "control"
      }
      to_annotated <- paste(tmp_name, tmp_kind, sep = "_")
      cnv_id <- paste(cnv$sample, cnv$Kind_sample, sep = "_")
      tmp_id <- which(cnv_id == to_annotated)
      #start annotate
      if(length(tmp_id) > 0){
        tmp_annot <- read_tsv(paste(path_annotsv, file_annotsv[i], sep = ""))
        annot_exist <- annot_exist + 1
        if(annot_exist == 1){
          data_annot <- data.frame(matrix(NA, ncol = dim(tmp_annot)[2], nrow = dim(cnv)[1]))
        }
        check_if_present <- paste(tmp_annot$`SV chrom`, tmp_annot$`SV start`, tmp_annot$`SV end`, sep = "_")
        tmp_id_to_take <- which(check_if_present %in% cnv$id[tmp_id])
        a <- tmp_annot[tmp_id_to_take,]
        a <- a[which(a$`AnnotSV type` == "full"),]
        tmp_id_to_replace <- which(cnv$id[tmp_id] %in% check_if_present)
        if(dim(a)[1] > 0){
          data_annot[tmp_id[tmp_id_to_replace],] <- a
        }
        
      }
    }
    colnames(data_annot) <- colnames(tmp_annot)
    cnv_annoted <- cbind(cnv, data_annot)
    write_xlsx(cnv_annoted, paste(path_output, "combined_cnv_annotated.xlsx", sep = ""))
  }
  
  
  ### Combining CNV into larger region and see overlaps
  #remove ASc and control
  cnv <- cnv[-which(cnv$Kind_sample %in% c("Asc", "control")),]
  
  range <- split(IRanges(start = cnv$start, end = cnv$end), cnv$chr)
  #coloring <- paste(cnv[which(cnv$chr == "14"),]$sample, cnv[which(cnv$chr == "14"),]$Kind_sample, sep = "_")
  
  
  plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                         col="black", sep=0.5, ...){
    height <- 1
    if (is(xlim, "IntegerRanges")){
      xlim <- c(min(start(xlim)), max(end(xlim)))
    }
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    png(filename = paste(path_output, main, "_overlap_cnv.png", sep = ""),
        width = 1800, height = 800)
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
    title(main)
    axis(1)
    dev.off()
  }
  
  for (i in 1:length(range)) {
    coloring <- cnv[which(cnv$chr == names(range)[i]),]$Kind_sample
    palette <- rainbow(length(unique(coloring)))
    col_patient <- data.frame(matrix(NA, ncol = 2, nrow = length(coloring)))
    colnames(col_patient) <- c("patient", "color")
    col_patient$patient <- coloring
    col_patient$color[which(col_patient$patient == "TuPr")] <- palette[1]
    col_patient$color[which(col_patient$patient == "TuPo")] <- palette[2]
    #not always same order for groups on different chromosomes
    #for (j in 1:length(unique(coloring))) {
    #  col_patient$color[which(col_patient$patient == unique(coloring)[j])] <- palette[j]
    #}
    plotRanges(range[[i]], main = paste("Chr", names(range)[i]), col = col_patient$color) 
    reduce_range <- reduce(range[[i]])
  }
  #TO DO combined overlapping cnv
  
  ###Map gene to cnv
  ensembl = useEnsembl(biomart='ensembl', 
                       dataset="hsapiens_gene_ensembl")
  
  results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                   filters = c("chromosome_name", "start", "end"),
                   values = list(cnv$chr, cnv$start, cnv$end),
                   mart = ensembl)
  #remove chromosomal region without gene annotations
  results <- results[which(!results$hgnc_symbol == ""),]
  range_gene <- split(IRanges(start = results$start_position, end = results$end_position), results$chromosome_name)
  for (i in 1:length(range_gene)) {
    plotRanges(range_gene[[i]], main = paste("Gene_on_chr_", names(range_gene)[i])) 
  }
  #TODO map the gene to the cnv
  
  #return()
}