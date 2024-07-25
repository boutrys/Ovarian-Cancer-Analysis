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

#' Facets_analysis
#'
#' Given Facets data for each patients
#' allow to analysis of regions of LOH/gain
#' plot these results and computes statistics and annotate clinical data with the different cni score 
#' 
#' @param path where to retrieve FACETS data
#' @param path_output where to store results
#' @param path_analysis where to store the new annotated clinical data
#'
#'

Facets_analysis <- function(path = getwd(),
                            to_change = "",
                            new_val = "",
                            clinical_data_path = c(),
                            path_output = getwd()){
  #libraries
  library(biomaRt)
  library(tidyverse)
  library(car)
  library(readxl)
  library(writexl)
  library(plotly)
  library(cowplot)
  library(ggpubr)
  library(stats)

  
  #create output directory
  dir.create(path_output)
  
  #Load ensembl database
  ensembl = useEnsembl(biomart="ensembl", 
                       dataset="hsapiens_gene_ensembl")
  
  
  #Create the column names for the reporting files
  facets_cnv_general <- data.frame(Patient=character(),
                                   n_gains=character(),
                                   n_of_homoz_losses=character(),
                                   n_of_homoz_losses.1MB=character(),
                                   n_of_heteroz_losses=character(),
                                   n_of_heteroz_losses.1MB=character(),
                                   AMPs=character(),
                                   HOMDELs=character(),
                                   HOMDELs.1MB=character(),
                                   LOHs=character(),
                                   LOHs.1MB=character(),
                                   stringsAsFactors=FALSE)
  cols_DEL_LOH_gains <- data.frame(Patient=character(),
                                   Chromosome=character(),
                                   Location_of_segment=character(),
                                   CN_interpretation=character(),
                                   Relative_size=character(),
                                   genes_of_interest=character(),
                                   stringsAsFactors=FALSE)
  cols_individualized_genes <- data.frame(Patient=character(),
                                          Chromosome=character(),
                                          Location_of_segment=character(),
                                          CN_interpretation=character(),
                                          Relative_size=character(),
                                          genes_of_interest=character(),
                                          stringsAsFactors=FALSE)
  
  #Load clinical and explore SCORE files 
  clinical_data <- read_xlsx(clinical_data_path)
  
  full_path_list <- list.files(path = path, 
                               pattern = "scores", 
                               full.names = TRUE)
  clinical_data$ASC_ploidy <- NA
  clinical_data$TuPr_ploidy <- NA
  clinical_data$TuPo_ploidy <- NA
  clinical_data$ASC_purity <- NA
  clinical_data$TuPr_purity <- NA
  clinical_data$TuPo_purity <- NA
  clinical_data$ASC_WGD <- NA
  clinical_data$TuPr_WGD <- NA
  clinical_data$TuPo_WGD <- NA
  for (i in 1:length(full_path_list)) {
    score_file <- read_tsv(full_path_list[i], col_names = FALSE)
    tmp_patient_full <- score_file$X1[1]
    
    #Annotate clinical data
    tmp_patient <- gsub("VDB-", "", tmp_patient_full)
    need_to_be_changed <- which(to_change == tmp_patient)
    if(length(need_to_be_changed)>0){
      tmp_patient <- new_val[need_to_be_changed]
    }
    tmp_split <- strsplit(tmp_patient, split = "-")
    tmp_patient <- tmp_split[[1]][1]
    tmp_group <- tmp_split[[1]][2]
    id_col_to_fill <- which(colnames(clinical_data) %in% paste(tmp_group, c("ploidy", "purity", "WGD"), sep = "_"))
    clinical_data[which(clinical_data$Patient == tmp_patient),id_col_to_fill[1]] <- as.double(score_file$X1[which(score_file$X1 == "ploidy") + 1])
    clinical_data[which(clinical_data$Patient == tmp_patient),id_col_to_fill[2]] <- as.double(score_file$X1[which(score_file$X1 == "purity") + 1])
    if(score_file$X1[which(score_file$X1 == "wgd_boolean") + 1] == "False"){
      tmp_WGD <- "No"
    }else if(score_file$X1[which(score_file$X1 == "wgd_boolean") + 1] == "True"){
      tmp_WGD <- "Yes"
    }
    clinical_data[which(clinical_data$Patient == tmp_patient),id_col_to_fill[3]] <- tmp_WGD
  }
  to_rm <- strsplit(clinical_data_path, split = "/")[[1]]
  tmp_path <- gsub(to_rm[length(to_rm)], "", clinical_data_path)
  write_xlsx(clinical_data, paste(tmp_path, "annotated_clinical_data_Facets.xlsx", sep = ""))
  
  
  #Set the location of the files to process
  to_process <- list.files(path = path, 
                           pattern = "tsv", 
                           full.names = TRUE)
  # Set the names of the Patients
  name_patient <- list.files(path = path, 
                             pattern = "tsv")
  name_patient <- strsplit(name_patient, ".facets.tsv")
  a <- 0
  
  # For each tsv file, process every calculation
  for (filename in to_process) {
    facets_data <- filename
    a <- a + 1
    
    #use factes data per patient
    facets_data <- read_tsv(facets_data)
    facets_data <- na.omit(facets_data)
    
    #TODO problem with name patient
    facets_data <- facets_data %>% 
      mutate( size = end - start) %>% 
      mutate( Patient = name_patient[a]) %>% 
      dplyr::select(Patient, chrom, start, end, size, tcn.em, lcn.em)
    
    #General data about losses and gains
    #Gains
    facets_gains <- facets_data %>% 
      dplyr::filter(tcn.em > 5)
    n_of_gains <- dim(facets_gains) [1]
    AMPs <- sum(facets_gains$size)
    #Homozygous losses
    #total
    facets_homoz_losses <- facets_data %>% 
      dplyr::filter(tcn.em == 0)
    n_of_homoz_losses <- dim(facets_homoz_losses) [1]
    HOMDELs <- sum(facets_homoz_losses$size)
    #Larger than 1 MB
    facets_homoz_losses_MB <- facets_homoz_losses %>% 
      dplyr::filter(size > 1000000)
    n_of_homoz_losses_MB <- dim(facets_homoz_losses_MB) [1]
    HOMDELs_MB <- sum(facets_homoz_losses_MB$size)
    #Retrieve Ploidy 
    tmp_patient <- strsplit(unique(facets_data$Patient)[[1]], split = "-")[[1]]
    tmp_ploidy <- round(clinical_data[which(clinical_data$Patient == tmp_patient[2]),which(colnames(clinical_data) == paste(tmp_patient[3], "ploidy", sep = "_"))], 1)
    #Heterozygous losses
    if(tmp_ploidy < 3.5){
      facets_heteroz_losses <- facets_data %>% 
        dplyr::filter(tcn.em == 1, tcn.em != 2*lcn.em)
    }else if(tmp_ploidy >= 3.5 && tmp_ploidy < 4.5){
      facets_heteroz_losses <- facets_data %>% 
        dplyr::filter(tcn.em < 3, tcn.em != 2*lcn.em)
    }else if(tmp_ploidy >= 4.5 && tmp_ploidy < 5.5){
      facets_heteroz_losses <- facets_data %>% 
        dplyr::filter(tcn.em < 4, tcn.em != 2*lcn.em)
    }else if(tmp_ploidy >= 5.5){
      facets_heteroz_losses <- facets_data %>% 
        dplyr::filter(tcn.em < 5, tcn.em != 2*lcn.em)
    }
    n_of_heteroz_losses <- dim(facets_heteroz_losses) [1]
    LOHs <- sum(facets_heteroz_losses$size)
    #Larger than 1 MB
    facets_heteroz_losses_MB <- facets_heteroz_losses %>% 
      dplyr::filter(size > 1000000)
    n_of_heteroz_losses_MB <- dim(facets_heteroz_losses_MB) [1]
    LOHs_MB <- sum(facets_heteroz_losses_MB$size)
    
    #Agregate data in 1 table
    gains_losses_general_data <- data.frame(n_gains = n_of_gains)
    gains_losses_general_data <- gains_losses_general_data %>% 
      mutate(n_of_homoz_losses = n_of_homoz_losses) %>% 
      mutate(n_of_homoz_losses.1MB = n_of_homoz_losses_MB) %>% 
      mutate(n_of_heteroz_losses = n_of_heteroz_losses) %>% 
      mutate(n_of_heteroz_losses.1MB = n_of_heteroz_losses_MB) %>% 
      mutate(AMPs = AMPs / 1000000) %>% 
      mutate(HOMDELs = HOMDELs / 1000000) %>% 
      mutate(HOMDELs.1MB = HOMDELs_MB / 1000000) %>% 
      mutate(LOHs = LOHs / 1000000) %>% 
      mutate(LOHs.1MB = LOHs_MB / 1000000)
    gains_losses_general_data <- round(gains_losses_general_data, 0)
    gains_losses_general_data <- gains_losses_general_data %>% 
      mutate(Patient = name_patient[a]) %>% 
      dplyr::select(Patient, everything())
    
    #Consider every segment with gain or loss 
    #What is the interpretation of the CN status of the segment ?
    facets_gains <- mutate(facets_gains, CN_interpretation = "gain")
    facets_homoz_losses <- mutate(facets_homoz_losses, CN_interpretation = "del")
    facets_heteroz_losses <- mutate(facets_heteroz_losses, CN_interpretation = "LOH")
    facets_genes <- rbind(facets_gains, facets_homoz_losses, facets_heteroz_losses)

    if(dim(facets_genes)[1] > 0){
      #for each segment, search for interesting genes it contains
      all_retained_gene <- getBM(c("external_gene_name"), 
                                 filters=c("chromosome_name",
                                           "start",
                                           "end"), 
                                 values = list(facets_genes$chrom, 
                                               facets_genes$start,
                                               facets_genes$end), 
                                 attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                                 mart=ensembl,
                                 useCache = F)
      for (i in 1:dim(facets_genes)[1]) {
        tmp <- all_retained_gene[which(all_retained_gene$chromosome_name == facets_genes$chrom[i]),]
        tmp <- tmp[which(tmp$start_position >= facets_genes$start[i]),]
        the_retained_genes <- tmp[which(tmp$end_position <= facets_genes$end[i]),]
        genes <- paste(the_retained_genes$external_gene_name, collapse = " ; ")
        facets_genes[i, "genes_of_interest"] <- genes
      }
      
      #GRCH38 data about centromere position on each chromosome
      ref.dat = data.frame( chromosome_ref = seq(1:23),
                            centromere = c(122026460, 92188146, 90772459, 49708101, 46485901, 58553889, 58169654, 44033745,
                                           43236168, 39686683, 51078349, 34769408, 16000001, 16000001, 17000001, 36311159,
                                           22813680, 15460900, 24498981, 26436233, 10864561, 12954789, 58605580),
                            chr.size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 
                                         145138636, 138394717,  133797422, 135086622, 133275309, 114364328, 107043718, 
                                         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 
                                         50818468, 156040895))
      
      #Position the start and end of each segment in relation to the arm of the chromosome
      for (i in 1:dim(facets_genes)[1]) {
        ref.dat2 <- dplyr::filter(ref.dat, chromosome_ref == facets_genes$chrom[i])
        segment_start <- ref.dat2$centromere - facets_genes$start[i]
        segment_start <- ifelse(segment_start > 0, "1", "2")
        segment_end <- ref.dat2$centromere - facets_genes$end[i]
        segment_end <- ifelse(segment_end > 0, "1", "2")
        facets_genes[i, "arm_start"] <- segment_start
        facets_genes[i, "arm_end"] <- segment_end
      }
      
      #Is the segment on arm p, q, or both ?
      facets_genes <- mutate(facets_genes, segment_arm = NA)
      for (i in 1:dim(facets_genes)[1]) {
        if((facets_genes$arm_start[i] == 1) & (facets_genes$arm_end[i] == 1)) 
        {facets_genes$segment_arm[i] <- "p"
        } 
        else if( (facets_genes$arm_start[i] == 1) & (facets_genes$arm_end[i] == 2)) {
          facets_genes$segment_arm[i] <- "pq"
        } 
        else {facets_genes$segment_arm[i] <- "q"}
      }
      
      #what is the size of the segment in comparison to the chromosome ?
      facets_genes <- mutate(facets_genes, segment_rel_size = NA)
      for (i in 1:dim(facets_genes)[1]) {
        ref.dat3 <- dplyr::filter(ref.dat, chromosome_ref == facets_genes$chrom[i])
        facets_genes$segment_rel_size[i] <- (facets_genes$size[i] / ref.dat3$chr.size)
      }
      facets_genes$segment_rel_size <- round(facets_genes$segment_rel_size, 2)
      facets_genes <- facets_genes %>% 
        dplyr::select(Patient, chrom, segment_arm, CN_interpretation, segment_rel_size, `genes_of_interest`)
      colnames(facets_genes)[which(colnames(facets_genes) == "chrom")] <- "Chromosome"
      colnames(facets_genes)[which(colnames(facets_genes) == "segment_arm")] <- "Location_of_segment"
      colnames(facets_genes)[which(colnames(facets_genes) == "segment_rel_size")] <- "Relative_size"

      
      #Individualize the genes of interest for each segment
      indiv_genes <- data.frame()
      for(i in 1:dim(facets_genes)[1]){
        if(dim(as.data.frame(strsplit(facets_genes$genes_of_interest[i], ";")))[1] > 1){
          to_repeat <- facets_genes[i,]
          new_list <- unlist(strsplit(facets_genes$genes_of_interest[i], ";"))
          inter_data <- to_repeat[rep(row.names(to_repeat), length(new_list)), ]
          inter_data$genes_of_interest <- new_list
          indiv_genes <- bind_rows(indiv_genes, inter_data)
          inter_data <- data.frame()
        }else{
          indiv_genes <- bind_rows(indiv_genes, facets_genes[i,])
        } 
        indiv_genes$genes_of_interest <- gsub(" ", "", indiv_genes$genes_of_interest)
      }
      
      #summary
      #facets_cnv_general <- rbind(facets_cnv_general, gains_losses_general_data)
      cols_DEL_LOH_gains <- rbind(cols_DEL_LOH_gains, facets_genes)
      cols_individualized_genes <- rbind(cols_individualized_genes, indiv_genes)
      
    }#end if facets_gene NOT empty
    facets_cnv_general <- rbind(facets_cnv_general, gains_losses_general_data)
  }#end for loop going through all patients
  
  facets_cnv_general$Patient <- unlist(facets_cnv_general$Patient)
  cols_DEL_LOH_gains$Patient <- unlist(cols_DEL_LOH_gains$Patient)
  cols_individualized_genes$Patient <- unlist(cols_individualized_genes$Patient)
  
  #write the results in tsv files
  write_tsv(facets_cnv_general, 
            paste(path_output, "patients_general.tsv", sep = ""),
            na = "NA",
            col_names = TRUE,
            append = FALSE)
  write_tsv(cols_DEL_LOH_gains, 
            paste(path_output, "patients_genes.tsv", sep = ""),
            na = "NA",
            col_names = TRUE,
            append = FALSE)
  write_tsv(cols_individualized_genes, 
            paste(path_output, "patients_indiv genes.tsv", sep = ""),
            na = "NA",
            col_names = TRUE,
            append = FALSE)
  
  #annotate with fraction genome cna altered
  facets_cnv_general$fraction_genome_CNA_altered <- NA
  for (i in 1:length(full_path_list)) {
    score_file <- read_tsv(full_path_list[i], col_names = FALSE)
    tmp_patient_full <- score_file$X1[1]
    #add to facets data the fraction genome cna alteret for latter correlation analysis
    tmp_id <- which(facets_cnv_general$Patient %in% tmp_patient_full)
    if(length(tmp_id)> 0){
      facets_cnv_general$fraction_genome_CNA_altered[tmp_id] <- as.double(score_file$X1[which(score_file$X1 == "fraction_genome_CNA_altered") + 1])
    }
    
    #Annotate clinical data
    tmp_patient <- gsub("VDB-", "", tmp_patient_full)
    need_to_be_changed <- which(to_change == tmp_patient)
    if(length(need_to_be_changed)>0){
      tmp_patient <- new_val[need_to_be_changed]
      
      #add to facets data the fraction genome cna alteret for latter correlation analysis
      list_patient_facets <- gsub("VDB-", "", facets_cnv_general$Patient)
      tmp_id <- which(list_patient_facets == tmp_patient)
      if(length(tmp_id)> 0){
        facets_cnv_general$fraction_genome_CNA_altered[tmp_id] <- as.double(score_file$X1[which(score_file$X1 == "fraction_genome_CNA_altered") + 1])
      }
    }
  }
  
  
  #pre treatment clinical data
  clinical_data$DFS[which(clinical_data$DFS < 13)] <- "=<12"
  clinical_data$DFS[which(clinical_data$DFS > 12)] <- ">12"
  
  #Facets data treatment
  facets_cnv_general$Patient <- gsub("VDB-", "", facets_cnv_general$Patient)
  facets_cnv_general$group <- NA
  facets_cnv_general$Anapath_response <- NA
  facets_cnv_general$DFS <- NA
  facets_cnv_general$PFI <- NA
  for (i in 1:dim(facets_cnv_general)[1]) {
    tmp <- strsplit(facets_cnv_general$Patient[i], split = "-")[[1]]
    facets_cnv_general$group[i] <- tmp[2]
    facets_cnv_general$Anapath_response[i] <- clinical_data$Anapath_response[which(clinical_data$Patient == tmp[1])]
    facets_cnv_general$DFS[i] <- clinical_data$DFS[which(clinical_data$Patient == tmp[1])]
    facets_cnv_general$PFI[i] <- clinical_data$PFI[which(clinical_data$Patient == tmp[1])]
    facets_cnv_general$Patient[i] <- tmp[1]
  }
  colnames(facets_cnv_general)[which(colnames(facets_cnv_general) == "Patient")] <- "patient"
  
  #aneuploidy size computation
  facets_cnv_general$Aneuploidy <- facets_cnv_general$AMPs + facets_cnv_general$HOMDELs + facets_cnv_general$LOHs
  #put all predictor as ratio 
  facets_cnv_general$Aneuploidy <- facets_cnv_general$Aneuploidy  / 3000
  facets_cnv_general$AMPs <- facets_cnv_general$AMPs / 3000
  facets_cnv_general$HOMDELs <- facets_cnv_general$HOMDELs / 3000
  facets_cnv_general$LOHs <- facets_cnv_general$LOHs / 3000
  
  
  ### Correlation inbetween aneuploidy and fraction_genome_CNA_altered
  #The Pearson correlation
  
  tmp_cor <- cor.test(facets_cnv_general$fraction_genome_CNA_altered, facets_cnv_general$Aneuploidy, method = "spearman")
  g <- ggplot(facets_cnv_general, aes(x = fraction_genome_CNA_altered, y = Aneuploidy)) +
    geom_point(size = 2) +
    geom_smooth(method="loess", se=F) +
    labs(title = paste("Correlation full data equal", round(tmp_cor$estimate, digits = 4), "with associated pvalue =", round(tmp_cor$p.value, digits = 4))) +
    theme_classic()
  ggsave(g, filename = paste(path_output, "Correlation_Full_fraction_genome_CNA_altered_Aneuploidy.png", sep = ""), height = 8, width = 8)
  
  corr_Asc <- facets_cnv_general[which(facets_cnv_general$group == "ASC"),]
  tmp_cor <- cor.test(corr_Asc$fraction_genome_CNA_altered, corr_Asc$Aneuploidy, method = "spearman")
  g <- ggplot(corr_Asc, aes(x = fraction_genome_CNA_altered, y = Aneuploidy)) +
    geom_point(size = 2) +
    geom_smooth(method="loess", se=F) +
    labs(title = paste("Correlation ASC data equal", round(tmp_cor$estimate, digits = 4), "with associated pvalue =", round(tmp_cor$p.value, digits = 4))) +
    theme_classic()
  ggsave(g, filename = paste(path_output, "Correlation_ASC_fraction_genome_CNA_altered_Aneuploidy.png", sep = ""), height = 8, width = 8)
  
  corr_TuPr <- facets_cnv_general[which(facets_cnv_general$group == "TuPr"),]
  tmp_cor <- cor.test(corr_TuPr$fraction_genome_CNA_altered, corr_TuPr$Aneuploidy, method = "spearman")
  g <- ggplot(corr_TuPr, aes(x = fraction_genome_CNA_altered, y = Aneuploidy)) +
    geom_point(size = 2) +
    geom_smooth(method="loess", se=F) +
    labs(title = paste("Correlation TuPr data equal", round(tmp_cor$estimate, digits = 4), "with associated pvalue =", round(tmp_cor$p.value, digits = 4))) +
    theme_classic()
  ggsave(g, filename = paste(path_output, "Correlation_TuPr_fraction_genome_CNA_altered_Aneuploidy.png", sep = ""), height = 8, width = 8)
  
  corr_TuPo <- facets_cnv_general[which(facets_cnv_general$group == "TuPo"),]
  tmp_cor <- cor.test(corr_TuPo$fraction_genome_CNA_altered, corr_TuPo$Aneuploidy, method = "spearman")
  g <- ggplot(corr_TuPo, aes(x = fraction_genome_CNA_altered, y = Aneuploidy)) +
    geom_point(size = 2) +
    geom_smooth(method="loess", se=F) +
    labs(title = paste("Correlation TuPo data equal", round(tmp_cor$estimate, digits = 4), "with associated pvalue =", round(tmp_cor$p.value, digits = 4))) +
    theme_classic()
  ggsave(g, filename = paste(path_output, "Correlation_TuPo_fraction_genome_CNA_altered_Aneuploidy.png", sep = ""), height = 8, width = 8)
  
  
  #TODO here we remove ASC, but we could include them
  facets_cnv_general <- facets_cnv_general[-which(facets_cnv_general$group == "ASC"),] 
  
  ###Start loop for statistics and plot for each predictors
  facets_cnv_general$Anapath_response <- ordered(as.factor(facets_cnv_general$Anapath_response), c("Bad", "Intermediate", "Good"))
  facets_cnv_general$DFS <- ordered(as.factor(facets_cnv_general$DFS), c("=<12", ">12"))
  facets_cnv_general$PFI <- ordered(as.factor(facets_cnv_general$PFI), c("platinum resistant", "semi sensitive", "platinum sensitive"))
  facets_cnv_general$group <- ordered(as.factor(facets_cnv_general$group), c("TuPr", "TuPo"))
  var_predictor <- c("Aneuploidy", "AMPs", "HOMDELs", "LOHs")
  group_name <- c("TuPr", "TuPo")
  group_name <- ordered(as.factor(group_name), c("TuPr", "TuPo"))
  nbr_group <- length(group_name)
  var_of_interest <- c("Anapath_response", "DFS", "PFI")
  y_axis_max <- max(facets_cnv_general$Aneuploidy, na.rm = TRUE)
  y_axis_max <- y_axis_max * 0.5 + y_axis_max
  list_plot <- vector(mode = "list", length = length(var_predictor)*length(var_of_interest))
  count <- 1
  for (m in 1:length(var_predictor)) {
    data_plot <- facets_cnv_general
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
            y_max <- y_max * 0.1 + y_max
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
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 2, width = 0.1) +
        xlab(var_of_interest[i]) +
        ylab(paste("Fraction genome with", var_predictor[m])) +
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
      plot1 <- ggplot(data_plot_tupr, aes(y = predictor, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 2, width = 0.1) +
        xlab(var_of_interest[i]) +
        ylab(paste("Fraction genome with", var_predictor[m])) +
        scale_x_discrete(labels=names(tmp_x_lab)) +
        stat_pvalue_manual(to_plot_stat[which(to_plot_stat$group == "TuPr"),], label = paste(name_stat_test, "pvalue = {p.adj}")) +
        ylim(0,y_axis_max) +
        scale_x_discrete(labels=paste(names(tmp_x_lab_tupr), "=",tmp_x_lab_tupr)) +
        labs("TuPr") +
        theme_classic()
      plot2 <- ggplot(data_plot_tupo, aes(y = predictor, x = eval(parse(text=paste(var_of_interest[i]))), label = patient)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 2, width = 0.1) +
        xlab(var_of_interest[i]) +
        ylab(paste("Fraction genome with", var_predictor[m])) +
        scale_x_discrete(labels=names(tmp_x_lab)) +
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
  
  
  
  ##TODO temporary solution removed not paired
  count_table <- table(facets_cnv_general$patient)
  facets_cnv_general$statut <- NA
  for (i in 1:length(count_table)) {
    if(count_table[i] == length(group_name)){
      facets_cnv_general$statut[which(facets_cnv_general$patient == names(count_table[i]))] <- "paired"
    }else{
      facets_cnv_general$statut[which(facets_cnv_general$patient == names(count_table[i]))] <- "not_paired"
    }
  }
  not_present_in_facets_patient <- facets_cnv_general[which(facets_cnv_general$statut == "not_paired"),]
  facets_cnv_general_new <- facets_cnv_general[-which(facets_cnv_general$statut == "not_paired"),]
  write_xlsx(not_present_in_facets_patient, paste(path_output, "Not_present_in_facets_TuPrVSTuPo_analysis.xlsx", sep = ""))
  
  ###Box plots to compare TuPr VS TuPO for each level of a var of interest
  list_plot <- vector(mode = "list", length = length(var_predictor)*length(var_of_interest))
  count <- 1
  for (m in 1:length(var_predictor)) {
    data_plot <- facets_cnv_general_new
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
        ylab(paste("Fraction genome with", var_predictor[m])) +
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
  
  
  write_xlsx(facets_cnv_general, paste(path_output, "summary_facets_data.xlsx", sep = ""))
  
  #return()
}




