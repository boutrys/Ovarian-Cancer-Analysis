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
#NB
#Plots for group only show the best pvalue of the group, while all others are also significant. For Go term GO_Plot, group, only the most significant of the group is plot. 
#Potential problem not handle yet : some pathway might have name including "/". But we also use "/" to seprate name of pathway within group of pathway with same geneID! 
#Output : for the moment over_represented_GO is the term that have been found in one of the GO database at least, and over_represented_pathway are the term found in other database than GO
#coocurrence matrix and plot are by default no performed, be carefull, if gene list is big, it is computationaly intensive
#coocurrence plot only display the top 100 most present genes 
over_representation_analysis <- function(geneList = c(),
                                         go_of_interest = c("BP", "MF", "CC"),
                                         max_GO_similarity_accepted = 0.5,
                                         cutoff = 0.05,
                                         max_gene_per_item = 20,
                                         max_item_plot = 20,
                                         coocurrence = TRUE,
                                         max_coocurrence = 100,
                                         path_to_store_results = ""){
  start <- Sys.time()
  #Libraries
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
  library(writexl)
  library(tidyverse)
  library(ReactomePA)
  library(DOSE)
  library(ComplexHeatmap)
  require(circlize)
  
  
  Nbr_go_of_interest <- length(go_of_interest)
  dir.create(path_to_store_results)
  options(bitmapType="cairo")
  
  
  ##############################################################################################
  ############################################### GLOBAL GENE SET ANALYSIS   ###################
  ##############################################################################################
  ids_geneList <- bitr(geneList, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  
  if(length(unique(ids_geneList[,1])) < length(unique(geneList))){
    gene_not_handled <- setdiff(geneList, unique(ids_geneList[,1]))
    writexl::write_xlsx(data.frame(gene_not_handled), path = paste(path_to_store_results, "list_gene_not_handled_for_analysis.xlsx", sep = "/")) 
    geneList <- geneList[-which(geneList %in% gene_not_handled)]
  }
  
  
  #GO OVER-REPRESENTATION TEST (on unique gene list : see below to use original dataframe)
  path_GO_over_representation <- paste(path_to_store_results, "GO_over_representation", sep = "/")
  dir.create(path_GO_over_representation)
  short_ego_track <- list()
  for (j in 1:length(go_of_interest)) {
    current_folder <- paste(path_GO_over_representation, go_of_interest[j], sep = "/")
    dir.create(current_folder)
    current_folder <- paste(current_folder, "/", sep = "")
    ego <- clusterProfiler::enrichGO(gene = base::unique(ids_geneList[,2]),
                                     OrgDb = org.Hs.eg.db,
                                     ont = go_of_interest[j],
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = cutoff,
                                     readable = TRUE)
    ego_result <- ego@result
    
    
    #Keep only significant results
    short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)
    
    #Remove GO term with more than max_gene_per_item included (big GO term filter)
    id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
    if(length(id_too_big_rm) > 0){
      remove_too_big <-  short_ego[id_too_big_rm,]
      writexl::write_xlsx(remove_too_big, path = paste(current_folder, "/", go_of_interest[j], "_GO_removed_too_big.xlsx", sep = ""))
      short_ego <- short_ego[-id_too_big_rm,] 
    }
    
    #Remove terms with only one gene from the gene list
    short_ego <- short_ego[which(short_ego$Count > 1),]
    
    
    #Simplify results
    if(dim(short_ego)[1] > 0){
      #Simplify the GO terms, by removing redundant terms based on their semantic similarities
      hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=go_of_interest[j])
      
      similarity_matrix  <- mgoSim(short_ego$ID, 
                                   short_ego$ID, 
                                   semData = hsGO, 
                                   combine = NULL)
      for (l in 1:dim(similarity_matrix)[1]) {
        for (k in 1:dim(similarity_matrix)[2]) {
          if(l == k){
            similarity_matrix[l,k] <- -1
          }
        }
      }
      over_similarity <- base::colSums(similarity_matrix > max_GO_similarity_accepted)
      tmp_over_similarity <- over_similarity
      tmp_similarity_matrix <- similarity_matrix
      current_max_sim <- max(tmp_similarity_matrix)
      idx_to_remove <- c()
      nbr_loop <- 0
      while (current_max_sim > max_GO_similarity_accepted) {
        nbr_loop <- nbr_loop + 1
        current_idx <- which(tmp_over_similarity == max(tmp_over_similarity))
        if(length(current_idx) > 1){
          current_idx <- current_idx[length(current_idx)]
        }
        tmp_idx <- which(short_ego$ID == colnames(tmp_similarity_matrix)[current_idx])
        tmp_similarity_matrix <- tmp_similarity_matrix[-current_idx,-current_idx]
        if(is.null(dim(tmp_similarity_matrix))){
          tmp_over_similarity <- 0
        }else{
          tmp_over_similarity <- base::colSums(tmp_similarity_matrix > max_GO_similarity_accepted)
        }
        idx_to_remove <- c(idx_to_remove, tmp_idx)
        current_max_sim <- max(tmp_similarity_matrix)
      }
      if(length(idx_to_remove) > 0){
        removed_from_short_ego <- short_ego[idx_to_remove,]
        short_ego <- short_ego[-idx_to_remove,]
        writexl::write_xlsx(removed_from_short_ego, path = paste(current_folder, "/", go_of_interest[j], "_GO_removed_similar_representation.xlsx", sep = ""))
      }
      
      
      #Remove redundant geneID list from results, keep only most significant one
      duplicate_geneID <- duplicated(short_ego$geneID)
      to_rm <- which(duplicate_geneID)
      if(length(to_rm) > 0){
        removed_duplicate <- short_ego[to_rm,]
        short_ego <- short_ego[which(!duplicate_geneID),]
        writexl::write_xlsx(removed_duplicate, path = paste(current_folder, "/", go_of_interest[j], "_GO_removed_duplicated_gene_list.xlsx", sep = ""))
        
        #Create group for GO terms with same geneID list 
        dup_gene_list <- unique(removed_duplicate$geneID)
        nbr_dup_gene_list <- length(dup_gene_list)
        group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
        colnames(group_description) <- c("group_ID", "Nbr_GO_term", "geneID", "Description")
        for (n in 1:nbr_dup_gene_list) {
          new_name <- paste(short_ego$Description[which(short_ego$geneID == dup_gene_list[n])], 
                            paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                            sep = "/")
          tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
          short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- paste(go_of_interest[j], "group", n, sep = "_")
          group_description$group_ID[n] <- paste(go_of_interest[j], "group", n, sep = "_")
          group_description$Nbr_GO_term[n] <- tmp_nbr
          group_description$geneID[n] <- dup_gene_list[n]
          group_description$Description[n] <- new_name
        }
        writexl::write_xlsx(group_description, path = paste(current_folder, "/", go_of_interest[j], "_group_description.xlsx", sep = ""))
      }
      
      
      #Plotting the results 
      ego_simplified <- ego
      ego_simplified@result <- short_ego
      current_item_plot <- dim(ego_simplified@result)[1]
      if(current_item_plot > max_item_plot){
        current_item_plot <- max_item_plot
      }
      png(paste(current_folder, "Dot_plot_", go_of_interest[j], "_GO_overrepresentation.png", sep = ""),
          width = 1000,
          height = 1000)
      print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
      dev.off()
      
      #TODO why bug when only one term to plot
      if(dim(short_ego)[1] > 1){
        png(paste(current_folder, "GO_Plot_", go_of_interest[j], "_GO_overrepresentation.png", sep = ""),
            width = 1000,
            height = 1000)
        print(enrichplot::goplot(ego_simplified, showCategory = current_item_plot))
        dev.off() 
      }
      
      tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
      tmp_height <- 600 + current_item_plot * 10
      png(paste(current_folder, "Plot_Gene_occurence_GO_", go_of_interest[j], ".png", sep = ""),
          width = tmp_width,
          height = tmp_height)
      print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
      dev.off()
      
      
      #transform results if group are formed
      if(length(to_rm) > 0){
        short_ego$ID[grep(short_ego$Description, pattern = paste(go_of_interest[j], "group", sep = "_"))] <- short_ego$Description[grep(short_ego$Description, pattern = paste(go_of_interest[j], "group", sep = "_"))]
        short_ego$Database <- paste("GO", go_of_interest[j], sep = "_")
      }
      
      
      #Save results
      writexl::write_xlsx(short_ego, path = paste(current_folder, "/", go_of_interest[j], "_GO_over_representation.xlsx", sep = ""))
      
    }#IF dim(short_ego)[1] > 0
    
    short_ego_track[[j]] <- short_ego
  }#For loop for each BP, CC, MF GO terms
  
  
  ###Others Databas for pathway analysis
  pathway_analysis <- paste(path_to_store_results, "pathway_over_representation", sep = "/")
  dir.create(pathway_analysis)
  short_pathway_track <- list()
  
  
  ###KEGG analysis
  KEGG_path <- paste(pathway_analysis, "KEGG", sep = "/")
  dir.create(KEGG_path)
  ego <- clusterProfiler::enrichKEGG(gene = base::unique(ids_geneList[,2]),
                                     organism = "hsa",
                                     keyType = "kegg",
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = cutoff,
                                     use_internal_data = FALSE)
  ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
  ego_result <- ego@result
  
  
  #Keep only significant results
  short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)
  
  
  #Remove pathways with more than max_gene_per_item included (big pathway filter)
  id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
  if(length(id_too_big_rm) > 0){
    remove_too_big <-  short_ego[id_too_big_rm,]
    writexl::write_xlsx(remove_too_big, path = paste(KEGG_path, "KEGG_removed_too_big.xlsx", sep = "/"))
    short_ego <- short_ego[-id_too_big_rm,] 
  }
  
  #Remove terms with only one gene from the gene list
  short_ego <- short_ego[which(short_ego$Count > 1),]
  
  
  #Simplify results
  if(dim(short_ego)[1] > 0){
    #Remove redundant geneID list from results, keep only most significant one
    duplicate_geneID <- duplicated(short_ego$geneID)
    to_rm <- which(duplicate_geneID)
    if(length(to_rm) > 0){
      removed_duplicate <- short_ego[to_rm,]
      short_ego <- short_ego[which(!duplicate_geneID),]
      writexl::write_xlsx(removed_duplicate, path = paste(KEGG_path, "KEGG_removed_duplicated_gene_list.xlsx", sep = "/"))
      
      #Create group for pathway terms with same geneID list
      dup_gene_list <- unique(removed_duplicate$geneID)
      nbr_dup_gene_list <- length(dup_gene_list)
      group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
      colnames(group_description) <- c("group_ID", "Nbr_pathway", "geneID", "Description")
      for (n in 1:nbr_dup_gene_list) {
        new_name <- paste(short_ego$Description[which(short_ego$geneID == dup_gene_list[n])], 
                          paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                          sep = "/")
        tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
        short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- paste("KEGG_group", n, sep = "_")
        group_description$group_ID[n] <- paste("KEGG_group", n, sep = "_")
        group_description$Nbr_pathway[n] <- tmp_nbr
        group_description$geneID[n] <- dup_gene_list[n]
        group_description$Description[n] <- new_name
      }
      writexl::write_xlsx(group_description, path = paste(KEGG_path, "KEGG_group_description.xlsx", sep = "/"))
    }
    
    
    #Plotting the results 
    ego_simplified <- ego
    ego_simplified@result <- short_ego
    current_item_plot <- dim(ego_simplified@result)[1]
    if(current_item_plot > max_item_plot){
      current_item_plot <- max_item_plot
    }
    png(paste(KEGG_path, "Dot_plot_KEGG_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    png(paste(KEGG_path, "Network_KEGG_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    png(paste(KEGG_path, "Network_Circle_KEGG_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  circular = TRUE,
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
    tmp_height <- 600 + current_item_plot * 10
    png(paste(KEGG_path, "Plot_Gene_occurence_KEGG_pathway.png", sep = "/"),
        width = tmp_width,
        height = tmp_height)
    print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    
    #transform results if group are formed
    if(length(to_rm) > 0){
      short_ego$ID[grep(short_ego$Description, pattern = "KEGG_group_")] <- short_ego$Description[grep(short_ego$Description, pattern = "KEGG_group_")]
      short_ego$Database <- "Pathway_KEGG"
    }
    
    
    #Save results
    writexl::write_xlsx(short_ego, path = paste(KEGG_path, "KEGG_over_representation.xlsx", sep = "/"))
    
  }#IF dim(short_ego)[1] > 0
  
  short_pathway_track[[1]] <- short_ego
  
  
  ###WikiPathways analysis
  WikiPathways_path <- paste(pathway_analysis, "WikiPathways", sep = "/")
  dir.create(WikiPathways_path)
  ego <- clusterProfiler::enrichWP(gene = base::unique(ids_geneList[,2]),
                                   organism = "Homo sapiens",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = cutoff)
  ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
  ego_result <- ego@result
  
  
  #Keep only significant results
  short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)
  
  
  #Remove pathways with more than max_gene_per_item included (big pathway filter)
  id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
  if(length(id_too_big_rm) > 0){
    remove_too_big <-  short_ego[id_too_big_rm,]
    writexl::write_xlsx(remove_too_big, path = paste(WikiPathways_path, "WikiPathways_removed_too_big.xlsx", sep = "/"))
    short_ego <- short_ego[-id_too_big_rm,] 
  }
  
  #Remove terms with only one gene from the gene list
  short_ego <- short_ego[which(short_ego$Count > 1),]
  
  #Simplify results
  if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
    #Remove redundant geneID list from results, keep only most significant one
    duplicate_geneID <- duplicated(short_ego$geneID)
    to_rm <- which(duplicate_geneID)
    if(length(to_rm) > 0){
      removed_duplicate <- short_ego[to_rm,]
      short_ego <- short_ego[which(!duplicate_geneID),]
      writexl::write_xlsx(removed_duplicate, path = paste(WikiPathways_path, "WikiPathways_removed_duplicated_gene_list.xlsx", sep = "/"))
      
      #Create group for pathway terms with same geneID list
      dup_gene_list <- unique(removed_duplicate$geneID)
      nbr_dup_gene_list <- length(dup_gene_list)
      group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
      colnames(group_description) <- c("group_ID", "Nbr_pathway", "geneID", "Description")
      for (n in 1:nbr_dup_gene_list) {
        new_name <- paste(short_ego$Description[which(short_ego$geneID == dup_gene_list[n])], 
                          paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                          sep = "/")
        tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
        short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- paste("WikiPathways_group", n, sep = "_")
        group_description$group_ID[n] <- paste("WikiPathways_group", n, sep = "_")
        group_description$Nbr_pathway[n] <- tmp_nbr
        group_description$geneID[n] <- dup_gene_list[n]
        group_description$Description[n] <- new_name
      }
      writexl::write_xlsx(group_description, path = paste(WikiPathways_path, "WikiPathways_group_description.xlsx", sep = "/"))
    }
    
    
    #Plotting the results 
    ego_simplified <- ego
    ego_simplified@result <- short_ego
    current_item_plot <- dim(ego_simplified@result)[1]
    if(current_item_plot > max_item_plot){
      current_item_plot <- max_item_plot
    }
    png(paste(WikiPathways_path, "Dot_plot_WikiPathways_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    png(paste(WikiPathways_path, "Network_WikiPathways_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    png(paste(WikiPathways_path, "Network_Circle_WikiPathways_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  circular = TRUE,
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
    tmp_height <- 600 + current_item_plot * 10
    png(paste(WikiPathways_path, "Plot_Gene_occurence_WikiPathways_pathway.png", sep = "/"),
        width = tmp_width,
        height = tmp_height)
    print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    
    #transform results if group are formed
    if(length(to_rm) > 0){
      short_ego$ID[grep(short_ego$Description, pattern = "WikiPathways_group_")] <- short_ego$Description[grep(short_ego$Description, pattern = "WikiPathways_group_")]
      short_ego$Database <- "Pathway_WikiPathways"
    }
    
    
    #Save results
    writexl::write_xlsx(short_ego, path = paste(WikiPathways_path, "WikiPathways_over_representation.xlsx", sep = "/"))
    
  }#IF dim(short_ego)[1] > 0
  
  short_pathway_track[[2]] <- short_ego
  
  
  ###Reactome analysis
  Reactome_path <- paste(pathway_analysis, "Reactome", sep = "/")
  dir.create(Reactome_path)
  ego <- ReactomePA::enrichPathway(gene = base::unique(ids_geneList[,2]),
                                   organism = "human",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = cutoff,
                                   readable = TRUE)
  ego_result <- ego@result
  
  
  #Keep only significant results
  short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)
  
  
  #Remove pathways with more than max_gene_per_item included (big pathway filter)
  id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
  if(length(id_too_big_rm) > 0){
    remove_too_big <-  short_ego[id_too_big_rm,]
    writexl::write_xlsx(remove_too_big, path = paste(Reactome_path, "Reactome_removed_too_big.xlsx", sep = "/"))
    short_ego <- short_ego[-id_too_big_rm,] 
  }
  
  #Remove terms with only one gene from the gene list
  short_ego <- short_ego[which(short_ego$Count > 1),]
  
  
  #Simplify results
  if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
    #Remove redundant geneID list from results, keep only most significant one
    duplicate_geneID <- duplicated(short_ego$geneID)
    to_rm <- which(duplicate_geneID)
    if(length(to_rm) > 0){
      removed_duplicate <- short_ego[to_rm,]
      short_ego <- short_ego[which(!duplicate_geneID),]
      writexl::write_xlsx(removed_duplicate, path = paste(Reactome_path, "Reactome_removed_duplicated_gene_list.xlsx", sep = "/"))
      
      #Create group for pathway terms with same geneID list
      dup_gene_list <- unique(removed_duplicate$geneID)
      nbr_dup_gene_list <- length(dup_gene_list)
      group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
      colnames(group_description) <- c("group_ID", "Nbr_pathway", "geneID", "Description")
      for (n in 1:nbr_dup_gene_list) {
        new_name <- paste(short_ego$Description[which(short_ego$geneID == dup_gene_list[n])], 
                          paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                          sep = "/")
        tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
        short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- paste("Reactome_group", n, sep = "_")
        group_description$group_ID[n] <- paste("Reactome_group", n, sep = "_")
        group_description$Nbr_pathway[n] <- tmp_nbr
        group_description$geneID[n] <- dup_gene_list[n]
        group_description$Description[n] <- new_name
      }
      writexl::write_xlsx(group_description, path = paste(Reactome_path, "Reactome_group_description.xlsx", sep = "/"))
    }
    
    
    #Plotting the results 
    ego_simplified <- ego
    ego_simplified@result <- short_ego
    current_item_plot <- dim(ego_simplified@result)[1]
    if(current_item_plot > max_item_plot){
      current_item_plot <- max_item_plot
    }
    png(paste(Reactome_path, "Dot_plot_Reactome_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    png(paste(Reactome_path, "Network_Reactome_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    png(paste(Reactome_path, "Network_Circle_Reactome_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  circular = TRUE,
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
    tmp_height <- 600 + current_item_plot * 10
    png(paste(Reactome_path, "Plot_Gene_occurence_Reactome_pathway.png", sep = "/"),
        width = tmp_width,
        height = tmp_height)
    print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    
    #transform results if group are formed
    if(length(to_rm) > 0){
      short_ego$ID[grep(short_ego$Description, pattern = "Reactome_group_")] <- short_ego$Description[grep(short_ego$Description, pattern = "Reactome_group_")]
      short_ego$Database <- "Pathway_Reactome"
    }
    
    
    #Save results
    writexl::write_xlsx(short_ego, path = paste(Reactome_path, "Reactome_over_representation.xlsx", sep = "/"))
    
  }#IF dim(short_ego)[1] > 0
  
  short_pathway_track[[3]] <- short_ego
  
  
  ###DiseaseOntology analysis
  DiseaseOntology_path <- paste(pathway_analysis, "DiseaseOntology", sep = "/")
  dir.create(DiseaseOntology_path)
  ego <- DOSE::enrichDO(gene = base::unique(ids_geneList[,2]),
                        ont = "DO",
                        pAdjustMethod = "BH",
                        qvalueCutoff = cutoff,
                        readable = TRUE)
  ego_result <- ego@result
  
  
  #Keep only significant results
  short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)
  
  
  #Remove pathways with more than max_gene_per_item included (big pathway filter)
  id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
  if(length(id_too_big_rm) > 0){
    remove_too_big <-  short_ego[id_too_big_rm,]
    writexl::write_xlsx(remove_too_big, path = paste(DiseaseOntology_path, "DiseaseOntology_removed_too_big.xlsx", sep = "/"))
    short_ego <- short_ego[-id_too_big_rm,] 
  }
  
  #Remove terms with only one gene from the gene list
  short_ego <- short_ego[which(short_ego$Count > 1),]
  
  
  #Simplify results
  if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
    #Remove redundant geneID list from results, keep only most significant one
    duplicate_geneID <- duplicated(short_ego$geneID)
    to_rm <- which(duplicate_geneID)
    if(length(to_rm) > 0){
      removed_duplicate <- short_ego[to_rm,]
      short_ego <- short_ego[which(!duplicate_geneID),]
      writexl::write_xlsx(removed_duplicate, path = paste(DiseaseOntology_path, "DiseaseOntology_removed_duplicated_gene_list.xlsx", sep = "/"))
      
      #Create group for pathway terms with same geneID list
      dup_gene_list <- unique(removed_duplicate$geneID)
      nbr_dup_gene_list <- length(dup_gene_list)
      group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
      colnames(group_description) <- c("group_ID", "Nbr_pathway", "geneID", "Description")
      for (n in 1:nbr_dup_gene_list) {
        new_name <- paste(short_ego$Description[which(short_ego$geneID == dup_gene_list[n])], 
                          paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                          sep = "/")
        tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
        short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- paste("DiseaseOntology_group", n, sep = "_")
        group_description$group_ID[n] <- paste("DiseaseOntology_group", n, sep = "_")
        group_description$Nbr_pathway[n] <- tmp_nbr
        group_description$geneID[n] <- dup_gene_list[n]
        group_description$Description[n] <- new_name
      }
      writexl::write_xlsx(group_description, path = paste(DiseaseOntology_path, "DiseaseOntology_group_description.xlsx", sep = "/"))
    }
    
    
    #Plotting the results 
    ego_simplified <- ego
    ego_simplified@result <- short_ego
    current_item_plot <- dim(ego_simplified@result)[1]
    if(current_item_plot > max_item_plot){
      current_item_plot <- max_item_plot
    }
    png(paste(DiseaseOntology_path, "Dot_plot_DiseaseOntology_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    png(paste(DiseaseOntology_path, "Network_DiseaseOntology_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    png(paste(DiseaseOntology_path, "Network_Circle_DiseaseOntology_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  circular = TRUE,
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
    tmp_height <- 600 + current_item_plot * 10
    png(paste(DiseaseOntology_path, "Plot_Gene_occurence_DiseaseOntology_pathway.png", sep = "/"),
        width = tmp_width,
        height = tmp_height)
    print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    
    #transform results if group are formed
    if(length(to_rm) > 0){
      short_ego$ID[grep(short_ego$Description, pattern = "DiseaseOntology_group_")] <- short_ego$Description[grep(short_ego$Description, pattern = "DiseaseOntology_group_")]
      short_ego$Database <- "Pathway_DiseaseOntology"
    }
    
    
    #Save results
    writexl::write_xlsx(short_ego, path = paste(DiseaseOntology_path, "DiseaseOntology_over_representation.xlsx", sep = "/"))
    
  }#IF dim(short_ego)[1] > 0
  
  short_pathway_track[[4]] <- short_ego
  
  
  ###NetCanGen analysis : Network of cancer gene database
  NetCanGen_path <- paste(pathway_analysis, "NetCanGen", sep = "/")
  dir.create(NetCanGen_path)
  ego <- DOSE::enrichNCG(gene = base::unique(ids_geneList[,2]),
                         pAdjustMethod = "BH",
                         qvalueCutoff = cutoff,
                         readable = TRUE)
  ego_result <- ego@result
  
  
  #Keep only significant results
  short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)
  
  
  #Remove pathways with more than max_gene_per_item included (big pathway filter)
  id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
  if(length(id_too_big_rm) > 0){
    remove_too_big <-  short_ego[id_too_big_rm,]
    writexl::write_xlsx(remove_too_big, path = paste(NetCanGen_path, "NetCanGen_removed_too_big.xlsx", sep = "/"))
    short_ego <- short_ego[-id_too_big_rm,] 
  }
  
  #Remove terms with only one gene from the gene list
  short_ego <- short_ego[which(short_ego$Count > 1),]
  
  
  #Simplify results
  if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
    #Remove redundant geneID list from results, keep only most significant one
    duplicate_geneID <- duplicated(short_ego$geneID)
    to_rm <- which(duplicate_geneID)
    if(length(to_rm) > 0){
      removed_duplicate <- short_ego[to_rm,]
      short_ego <- short_ego[which(!duplicate_geneID),]
      writexl::write_xlsx(removed_duplicate, path = paste(NetCanGen_path, "NetCanGen_removed_duplicated_gene_list.xlsx", sep = "/"))
      
      #Create group for pathway terms with same geneID list
      dup_gene_list <- unique(removed_duplicate$geneID)
      nbr_dup_gene_list <- length(dup_gene_list)
      group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
      colnames(group_description) <- c("group_ID", "Nbr_pathway", "geneID", "Description")
      for (n in 1:nbr_dup_gene_list) {
        new_name <- paste(short_ego$Description[which(short_ego$geneID == dup_gene_list[n])], 
                          paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                          sep = "/")
        tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
        short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- paste("NetCanGen_group", n, sep = "_")
        group_description$group_ID[n] <- paste("NetCanGen_group", n, sep = "_")
        group_description$Nbr_pathway[n] <- tmp_nbr
        group_description$geneID[n] <- dup_gene_list[n]
        group_description$Description[n] <- new_name
      }
      writexl::write_xlsx(group_description, path = paste(NetCanGen_path, "NetCanGen_group_description.xlsx", sep = "/"))
    }
    
    
    #Plotting the results 
    ego_simplified <- ego
    ego_simplified@result <- short_ego
    current_item_plot <- dim(ego_simplified@result)[1]
    if(current_item_plot > max_item_plot){
      current_item_plot <- max_item_plot
    }
    png(paste(NetCanGen_path, "Dot_plot_NetCanGen_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    png(paste(NetCanGen_path, "Network_NetCanGen_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    png(paste(NetCanGen_path, "Network_Circle_NetCanGen_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  circular = TRUE,
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
    tmp_height <- 600 + current_item_plot * 10
    png(paste(NetCanGen_path, "Plot_Gene_occurence_NetCanGen_pathway.png", sep = "/"),
        width = tmp_width,
        height = tmp_height)
    print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    
    #transform results if group are formed
    if(length(to_rm) > 0){
      short_ego$ID[grep(short_ego$Description, pattern = "NetCanGen_group_")] <- short_ego$Description[grep(short_ego$Description, pattern = "NetCanGen_group_")]
      short_ego$Database <- "Pathway_NetCanGen"
    }
    
    
    #Save results
    writexl::write_xlsx(short_ego, path = paste(NetCanGen_path, "NetCanGen_over_representation.xlsx", sep = "/"))
    
  }#IF dim(short_ego)[1] > 0
  
  short_pathway_track[[5]] <- short_ego
  
  
  ###DisGeNET analysis : Disease Gene Network database
  DisGeNET_path <- paste(pathway_analysis, "DisGeNET", sep = "/")
  dir.create(DisGeNET_path)
  ego <- DOSE::enrichDGN(gene = base::unique(ids_geneList[,2]),
                         pAdjustMethod = "BH",
                         qvalueCutoff = cutoff,
                         readable = TRUE)
  ego_result <- ego@result
  
  
  #Keep only significant results
  short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)
  
  
  #Remove pathways with more than max_gene_per_item included (big pathway filter)
  id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
  if(length(id_too_big_rm) > 0){
    remove_too_big <-  short_ego[id_too_big_rm,]
    writexl::write_xlsx(remove_too_big, path = paste(DisGeNET_path, "DisGeNET_removed_too_big.xlsx", sep = "/"))
    short_ego <- short_ego[-id_too_big_rm,] 
  }
  
  #Remove terms with only one gene from the gene list
  short_ego <- short_ego[which(short_ego$Count > 1),]
  
  
  #Simplify results
  if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
    #Remove redundant geneID list from results, keep only most significant one
    duplicate_geneID <- duplicated(short_ego$geneID)
    to_rm <- which(duplicate_geneID)
    if(length(to_rm) > 0){
      removed_duplicate <- short_ego[to_rm,]
      short_ego <- short_ego[which(!duplicate_geneID),]
      writexl::write_xlsx(removed_duplicate, path = paste(DisGeNET_path, "DisGeNET_removed_duplicated_gene_list.xlsx", sep = "/"))
      
      #Create group for pathway terms with same geneID list
      dup_gene_list <- unique(removed_duplicate$geneID)
      nbr_dup_gene_list <- length(dup_gene_list)
      group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
      colnames(group_description) <- c("group_ID", "Nbr_pathway", "geneID", "Description")
      for (n in 1:nbr_dup_gene_list) {
        new_name <- paste(short_ego$Description[which(short_ego$geneID == dup_gene_list[n])], 
                          paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                          sep = "/")
        tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
        short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- paste("DisGeNET_group", n, sep = "_")
        group_description$group_ID[n] <- paste("DisGeNET_group", n, sep = "_")
        group_description$Nbr_pathway[n] <- tmp_nbr
        group_description$geneID[n] <- dup_gene_list[n]
        group_description$Description[n] <- new_name
      }
      writexl::write_xlsx(group_description, path = paste(DisGeNET_path, "DisGeNET_group_description.xlsx", sep = "/"))
    }
    
    
    #Plotting the results 
    ego_simplified <- ego
    ego_simplified@result <- short_ego
    current_item_plot <- dim(ego_simplified@result)[1]
    if(current_item_plot > max_item_plot){
      current_item_plot <- max_item_plot
    }
    png(paste(DisGeNET_path, "Dot_plot_DisGeNET_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    png(paste(DisGeNET_path, "Network_DisGeNET_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    png(paste(DisGeNET_path, "Network_Circle_DisGeNET_overrepresentation.png", sep = "/"),
        width = 1000,
        height = 1000)
    p <- cnetplot(ego_simplified, 
                  showCategory = current_item_plot,
                  categorySize="pvalue", 
                  circular = TRUE,
                  colorEdge = TRUE, 
                  color_category='firebrick', 
                  color_gene='steelblue')
    print(p)
    dev.off()
    
    tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
    tmp_height <- 600 + current_item_plot * 10
    png(paste(DisGeNET_path, "Plot_Gene_occurence_DisGeNET_pathway.png", sep = "/"),
        width = tmp_width,
        height = tmp_height)
    print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
    dev.off()
    
    
    #transform results if group are formed
    if(length(to_rm) > 0){
      short_ego$ID[grep(short_ego$Description, pattern = "DisGeNET_group_")] <- short_ego$Description[grep(short_ego$Description, pattern = "DisGeNET_group_")]
      short_ego$Database <- "Pathway_DisGeNET"
    }
    
    
    #Save results
    writexl::write_xlsx(short_ego, path = paste(DisGeNET_path, "DisGeNET_over_representation.xlsx", sep = "/"))
    
  }#IF dim(short_ego)[1] > 0
  
  short_pathway_track[[6]] <- short_ego
  
  
  ###Results combination
  nbr_result <- length(short_ego_track) + length(short_pathway_track)
  name_column <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Database")
  pathway_database <- c("Pathway_KEGG", "Pathway_WikiPathways", "Pathway_Reactome", "Pathway_DiseaseOntology", "Pathway_NetCanGen", "Pathway_DisGeNET")
  combined_results <- data.frame(matrix(NA, ncol = length(name_column), nrow = 0))
  colnames(combined_results) <- name_column
  count <- 1
  for (i in 1:nbr_result) {
    if(i <= Nbr_go_of_interest){
      #Retrieving GO term results
      tmp <- short_ego_track[[i]]
      if(dim(tmp)[2] < 10 && dim(tmp)[1] > 0){
        tmp$Database <- paste("GO", go_of_interest[i], sep = "_")
      }
      if(dim(tmp)[1] > 0){
        combined_results <- rbind(combined_results, tmp) 
      }
    }else{
      #Retrieving Pathway results
      tmp <- short_pathway_track[[count]]
      if(dim(tmp)[2] < 10 && dim(tmp)[1] > 0){
        tmp$Database <- pathway_database[count]
      }
      if(dim(tmp)[1] > 0){
        combined_results <- rbind(combined_results, tmp) 
      }
      count <- count + 1
    }
  }
  
  
  #Remove redundant geneID list from combined results, keep only the first one
  duplicate_geneID <- duplicated(combined_results$geneID)
  to_rm <- which(duplicate_geneID)
  if(length(to_rm) > 0){
    removed_duplicate <- combined_results[to_rm,]
    combined_results <- combined_results[which(!duplicate_geneID),]
    writexl::write_xlsx(removed_duplicate, path = paste(path_to_store_results, "Combined_removed_duplicated_gene_list.xlsx", sep = "/"))
    
    #Create group for pathway terms with same geneID list
    dup_gene_list <- unique(removed_duplicate$geneID)
    nbr_dup_gene_list <- length(dup_gene_list)
    group_description <- data.frame(matrix(NA, ncol = 4, nrow = nbr_dup_gene_list))
    colnames(group_description) <- c("group_ID", "Nbr_database", "geneID", "Description")
    for (n in 1:nbr_dup_gene_list) {
      new_name <- paste(combined_results$Description[which(combined_results$geneID == dup_gene_list[n])], 
                        paste(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                        sep = "/")
      new_name_db <- paste(combined_results$Database[which(combined_results$geneID == dup_gene_list[n])], 
                           paste(removed_duplicate$Database[which(removed_duplicate$geneID == dup_gene_list[n])], collapse = "/"),
                           sep = "/")
      tmp_nbr <- length(removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]) + 1
      combined_results$Description[which(combined_results$geneID == dup_gene_list[n])] <- paste("Combined_group", n, sep = "_")
      combined_results$ID[which(combined_results$geneID == dup_gene_list[n])] <- paste("Combined_group", n, sep = "_")
      combined_results$Database[which(combined_results$geneID == dup_gene_list[n])] <- new_name_db
      group_description$group_ID[n] <- paste("Combined_group", n, sep = "_")
      group_description$Nbr_database[n] <- tmp_nbr
      group_description$geneID[n] <- dup_gene_list[n]
      group_description$Description[n] <- new_name
    }
    writexl::write_xlsx(group_description, path = paste(path_to_store_results, "Combined_group_description.xlsx", sep = "/"))
  }
  
  
  #Plotting the results 
  ego_simplified <- ego
  ego_simplified@result <- combined_results
  ego_simplified@result$ID <- rownames(ego_simplified@result)
  current_item_plot <- dim(ego_simplified@result)[1]
  if(current_item_plot > max_item_plot){
    current_item_plot <- max_item_plot
  }
  png(paste(path_to_store_results, "Dot_plot_combined_overrepresentation.png", sep = "/"),
      width = 1000,
      height = 1000)
  print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
  dev.off()
  
  png(paste(path_to_store_results, "Network_combined_overrepresentation.png", sep = "/"),
      width = 1000,
      height = 1000)
  p <- cnetplot(ego_simplified, 
                showCategory = current_item_plot,
                categorySize="pvalue", 
                colorEdge = TRUE, 
                color_category='firebrick', 
                color_gene='steelblue')
  print(p)
  dev.off()
  
  png(paste(path_to_store_results, "Network_Circle_combined_overrepresentation.png", sep = "/"),
      width = 1000,
      height = 1000)
  p <- cnetplot(ego_simplified, 
                showCategory = current_item_plot,
                categorySize="pvalue", 
                circular = TRUE,
                colorEdge = TRUE, 
                color_category='firebrick', 
                color_gene='steelblue')
  print(p)
  dev.off()
  
  tmp_width <- 600 + length(unique(strsplit(paste(combined_results$geneID, collapse = "/"), split = "/")[[1]])) * 10
  tmp_height <- 600 + current_item_plot * 10
  png(paste(path_to_store_results, "Plot_Gene_occurence_combined_pathway.png", sep = "/"),
      width = tmp_width,
      height = tmp_height)
  print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
  dev.off()
  
  
  ###Summary
  ORA_summary_col <- c("geneList", "GeneRatio", "BestQvalue", "ID", "Description", "Database")
  ORA_summary <- data.frame(matrix(NA, ncol = length(ORA_summary_col), nrow = dim(combined_results)[1]))
  colnames(ORA_summary) <- ORA_summary_col
  ORA_summary$geneList <- combined_results$geneID
  ORA_summary$GeneRatio <- combined_results$GeneRatio
  ORA_summary$BestQvalue <- combined_results$qvalue
  ORA_summary$ID <- combined_results$ID
  ORA_summary$Description <- combined_results$Description
  ORA_summary$Database <- combined_results$Database
  writexl::write_xlsx(ORA_summary, path = paste(path_to_store_results, "Summary_results_over_representation.xlsx", sep = "/"))
  
  
  if(coocurrence){
    start_matrix <- Sys.time()
    #restrict analysis on top most significant GO/pathway
    ORA_summary_best <- ORA_summary %>%
      arrange(desc(BestQvalue))
    if(dim(ORA_summary_best)[1] > max_coocurrence){
      ORA_summary_best <- ORA_summary_best[c(1:max_coocurrence),]
    }
    
    
    #By gene based on input geneList
    nbr_gene <- length(geneList)
    geneBygene_matrix <- matrix(0, ncol = nbr_gene, nrow = nbr_gene) 
    colnames(geneBygene_matrix) <- geneList
    rownames(geneBygene_matrix) <- geneList
    for (i in 1:dim(ORA_summary_best)[1]) {
      tmp_list_gene <- strsplit(ORA_summary_best$geneList[i], split = "/")[[1]]
      for (j in 1:length(tmp_list_gene)) {
        tmp_gene <- tmp_list_gene[j]
        geneBygene_matrix[which(rownames(geneBygene_matrix) == tmp_gene),which(colnames(geneBygene_matrix) %in% tmp_list_gene)] <- geneBygene_matrix[which(rownames(geneBygene_matrix) == tmp_gene),which(colnames(geneBygene_matrix) %in% tmp_list_gene)] + 1  
      }
    }
    id_rm <- which(colSums(geneBygene_matrix) == 0)
    if(length(id_rm) > 0){
      geneBygene_matrix <- geneBygene_matrix[-id_rm, -id_rm] 
    }
    geneBygene_data <- data.frame(Gene = rownames(geneBygene_matrix))
    geneBygene_data <- cbind(geneBygene_data, as.data.frame(geneBygene_matrix))
    writexl::write_xlsx(geneBygene_data, path = paste(path_to_store_results, "Cooccurence_genes.xlsx", sep = "/"))
    end_matrix <- Sys.time()
    
    
    start_plot <- Sys.time()
    if(dim(geneBygene_matrix)[2] > 100){
      max_id <- colSums(geneBygene_matrix)
      condition_plot <- min(max_id[order(max_id, decreasing = TRUE)][1:100])
      to_keep_plot <- which(max_id >= condition_plot)[1:100]
      geneBygene_matrix_reduce <- geneBygene_matrix[to_keep_plot,to_keep_plot]
    }else{
      geneBygene_matrix_reduce <- geneBygene_matrix
    }
    
    #ploting coocurrence of genes
    name_of_plot <- paste(path_to_store_results, "Coocurrence_genes.png", sep = "/")
    png(name_of_plot,
        width = 1500, height = 1200)
    
    h1 <- Heatmap(geneBygene_matrix_reduce,
                  heatmap_legend_param = list(
                    title = "CoOcurrence"),
                  column_km = 1,
                  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                    if(geneBygene_matrix_reduce[i, j] >= 1)
                      grid.text(sprintf("%d", geneBygene_matrix_reduce[i, j]), x, y, gp = gpar(fontsize = 10))
                  }
    )
    draw(h1)
    dev.off()
    
    end_plot <- Sys.time()
  }#enf of If coocurrence
  
  
  #Gene never over represented in original input
  id_rm <- which(!geneList %in% unique(strsplit(paste(combined_results$geneID, collapse = "/"), split = "/")[[1]]))
  if(length(id_rm) > 0){
    name_no_OR <- geneList[id_rm]
    writexl::write_xlsx(data.frame(Gene = name_no_OR), path = paste(path_to_store_results, "Gene_not_over_represented.xlsx", sep = "/")) 
  }
  
  
  #Returning results
  #GO terms
  over_represented_GO_id <- grep(pattern = "GO_", combined_results$Database)
  Description <- c(combined_results$Description[over_represented_GO_id])
  geneID <- c(combined_results$geneID[over_represented_GO_id])
  Database <- c(combined_results$Database[over_represented_GO_id])
  over_represented_GO <- list(Description, geneID, Database)
  #Pathways
  over_represented_pathway_id <- which(!seq(1:dim(combined_results)[1]) %in% over_represented_GO_id)
  Description <- c(combined_results$Description[over_represented_pathway_id])
  geneID <- c(combined_results$geneID[over_represented_pathway_id])
  Database <- c(combined_results$Database[over_represented_pathway_id])
  over_represented_pathway <- list(Description, geneID, Database)
  
  stop <- Sys.time()
  print(stop - start)
  return(list(over_represented_pathway,
              over_represented_GO))
  
}#end of function
