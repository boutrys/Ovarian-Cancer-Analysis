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

#' data_pretreatment
#'
#' given clinical data, transform the values into categories
#' computes final CA125 and PCI and remove irrelevent column and save the new clinical data into a excell file called new_clinical_data
#' 
#' @param path where to find the clinical_data.xlsx file on your computer  
#' @param to_remove name of the column you want to get rid of
#' @param output where to store the results
#'
#' @return clinical data with categorical variable, CA 125 ans PCI pre-post computed and put into categories and irrelevent columns removed
#'

data_pretreatment <- function(full_path = getwd(),
                              to_remove = c(),
                              output = getwd()){
  ###Load needed libraries
  library(readxl)
  library(writexl)
  library(gtsummary)
  library(webshot2)
  
  
  ###Load data
  data <- read_xlsx(paste(full_path, "clinical_data.xlsx", sep = ""))
  data_noFactor <- data
  
  
  ###treat data based on clinical legend
  #BRCA
  if(!is.na(to_remove)){
    data$BRCA[which(data$BRCA == 1)] <- "Yes"
    data$BRCA[which(data$BRCA == 2)] <- "No"
    data$BRCA[which(data$BRCA == 3)] <- "Unknown"
  }
  data_noFactor$BRCA <- data$BRCA
  data$BRCA <- as.factor(data$BRCA)
  
  #FIGO stage
  if(!is.na(to_remove)){
  data$FIGO_stage[which(data$FIGO_stage == 11)] <- "IIIB"
  data$FIGO_stage[which(data$FIGO_stage == 12)] <- "IIIc"
  data$FIGO_stage[which(data$FIGO_stage == 13)] <- "IVA"
  data$FIGO_stage[which(data$FIGO_stage == 14)] <- "IVB"
  }
  data_noFactor$FIGO_stage <- data$FIGO_stage
  data$FIGO_stage <- as.factor(data$FIGO_stage)
  
  #Anapath
  if(!is.na(to_remove)){
  data$anapath[which(data$anapath == 1)] <- "serous"
  data$anapath[which(data$anapath == 2)] <- "carcinosacome"
  }
  data_noFactor$anapath <- data$anapath
  data$anapath <- as.factor(data$anapath)
  
  #Chemo type
  if(!is.na(to_remove)){
  data$Chemo_type[which(data$Chemo_type == 1)] <- "Taxol"
  data$Chemo_type[which(data$Chemo_type == 2)] <- "Gemci"
  }
  data_noFactor$Chemo_type <- data$Chemo_type
  data$Chemo_type <- as.factor(data$Chemo_type)
  
  #Reponse biologique
  if(!is.na(to_remove)){
  data$Reponse_biologique[which(data$Reponse_biologique == 1)] <- "Good"
  data$Reponse_biologique[which(data$Reponse_biologique == 2)] <- "Intermediate"
  data$Reponse_biologique[which(data$Reponse_biologique == 3)] <- "Bad"
  data_noFactor$Reponse_biologique <- data$Reponse_biologique
  data$Reponse_biologique <- as.factor(data$Reponse_biologique)
  }
  
  #CA 125 PRE-POST
  if(!is.na(to_remove)){
  new_CA125_post <- data$CA125_post
  opt_CA <- which(data$CA125_post < 35)
  new_CA125_post[opt_CA] <- 1
  opt_CA_2 <- data$CA125_initial * 0.1 > data$CA125_post
  for (i in 1:length(new_CA125_post)) {
    if(is.na(new_CA125_post[i])){
      
    }else if(new_CA125_post[i] > 1){
      if(opt_CA_2[i]){
        new_CA125_post[i] <- 2
      }else{
        new_CA125_post[i] <- 3 
      }
    }
  }
  data$CA125_post <- new_CA125_post
  data$CA125_post[which(data$CA125_post == 1)] <- "Good"
  data$CA125_post[which(data$CA125_post == 2)] <- "Intermediate"
  data$CA125_post[which(data$CA125_post == 3)] <- "Bad"
  }
  data_noFactor$CA125_post <- data$CA125_post
  data$CA125_post <- as.factor(data$CA125_post)
  
  #Reponse per op
  if(!is.na(to_remove)){
  data$reponse_per_op[which(data$reponse_per_op == 1)] <- "Good"
  data$reponse_per_op[which(data$reponse_per_op == 2)] <- "Intermediate"
  data$reponse_per_op[which(data$reponse_per_op == 3)] <- "Bad"
  }
  
  #AINS per op
  data_noFactor$AINS_per_op <- data$AINS_per_op
  data$AINS_per_op <- as.factor(data$AINS_per_op)
  
  #Difference PCI PRE-POST
  if(!is.na(to_remove)){
  new_PCI_diff <- data$Diff_PCI
  new_PCI_diff[which(new_PCI_diff %in% c(2:10))] <- 2
  new_PCI_diff[which(new_PCI_diff < 2)] <- 3
  new_PCI_diff[which(new_PCI_diff > 10)] <- 1
  data$Diff_PCI <- new_PCI_diff
  data$Diff_PCI[which(data$Diff_PCI == 1)] <- "Good"
  data$Diff_PCI[which(data$Diff_PCI == 2)] <- "Intermediate"
  data$Diff_PCI[which(data$Diff_PCI == 3)] <- "Bad"
  
  #If missing values in Diff PCI, assign value from response_per_op if available
  missing_values <- which(is.na(data$Diff_PCI))
  data$Diff_PCI[missing_values] <- data$reponse_per_op[missing_values]
  data_noFactor$reponse_per_op <- data$reponse_per_op
  data$reponse_per_op <- as.factor(data$reponse_per_op)
  }
  data_noFactor$Diff_PCI <- data$Diff_PCI
  data$Diff_PCI <- as.factor(data$Diff_PCI)
  
  #Amount RD TODO when mathieu explained me the data
  if(!is.na(to_remove)){
  data$amount_RD[which(data$amount_RD == 0)] <- "NO residual"
  data$amount_RD[which(data$amount_RD == 1)] <- "Milimeter"
  data$amount_RD[which(data$amount_RD == 2)] <- "Supra milimeter"
  data$amount_RD[which(data$amount_RD == 3)] <- "+5cm"
  }
  data_noFactor$amount_RD <- data$amount_RD
  data$amount_RD <- as.factor(data$amount_RD)
  
  #Anapath response 
  if(!is.na(to_remove)){
  data$Anapath_response[which(data$Anapath_response == 1)] <- "Good"
  data$Anapath_response[which(data$Anapath_response == 2)] <- "Intermediate"
  data$Anapath_response[which(data$Anapath_response == 3)] <- "Bad"
  }
  data_noFactor$Anapath_response <- data$Anapath_response
  data$Anapath_response <- as.factor(data$Anapath_response) 
  
  #PFI
  if(!is.na(to_remove)){
  new_PFI <- data$PFI
  new_PFI[which(data$PFI %in% c(0:6))] <- "platinum resistant"
  new_PFI[which(data$PFI %in% c(7:12))] <- "semi sensitive"
  new_PFI[which(data$PFI > 12)] <- "platinum sensitive"
  data$PFI <- new_PFI
  }
  data_noFactor$PFI <- data$PFI
  data$PFI <- as.factor(data$PFI)
  
  #remove columns not use for the downstream analysis
  if(!is.na(to_remove)){
  data <- data[,-which(colnames(data) %in% to_remove)]
  data_noFactor <- data_noFactor[,-which(colnames(data_noFactor) %in% to_remove)]
  }
  
  #save data pretreated
  write_xlsx(data, paste(output, "clinical_data_used_for_analysis.xlsx", sep = ""))
  
  table_summary <- tbl_summary(data[,-1]) %>%
    as_gt() %>%
    gt::gtsave(filename = paste(output, "summary_clinical_data.pdf", sep = ""))

  
  #return the new clinical data
  return(list(data, data_noFactor))
}