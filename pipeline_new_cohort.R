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

#PLEASE READ TODO in this script before running it (you might find them by ctrl+F and search for TODO)
#' Pipeline 
#'
#' Within the pipeline folder, go to the input folder and create a new folder for your analysis, and assign that name to name_analysis in parameters here bellow
#' Make sure to put all necessary input into name_analysis folder 
#' data containing the clinical data must be named clinical.data (and has to be an .xlsx file) Make sure the columns of your file have the right column name (see input folder example)
#' 
#' 
#' @param old_path If you wanna reuse an existing output folder to store results please give path here bellow and set old_reuse = TRUE
#' @param old_reuse by default old_reuse = FALSE, and old_path is not use, a new ouput folder will be created automatically
#' 
#' @param do_pretreatment If you don't want to pre-treat your data, set to FALSE
#' @param name_analysis folder name that must be placed in the input folder within the Pipeline and that must contains all the files required (with the right names)
#' @param path_to_pipeline where did you put the pipeline folder on your computer (that folder must contain the following folders : data, functions, input, results)
#' folder functions must contain all necessary functions listed here bellow
#' data_pretreatment.R, 
#' folder input must contain all necessary input (see input folder example to know how to structure your input data)
#' @param clinical_annotation_to_remove name of the column in your clinical data that are irrelevent for your analysis
#' 
#' @param do_clustering If you don't want to perform a clustering analysis, set to FALSE
#' @param input_data Where you have stored the output file "intersection" from Raphael Helaers genomicRegionIntersection java program
#' @param clinical_data where to find the clinical data
#' @param nbr_group_cluster How many clustered do you want to split your data
#' 
#' @param do_survival If you don't want to perform a survival analysis, set to FALSE
#' @param variable_for_survival variable of the form c(event, time) (Must be in that order!!) in order to perform survival analysis based on these criteria
#' @param to_investigate_survival clinical annotation to investigate for the survival analysis
#' 
#' @param do_regression If you don't want to perform a regression analysis, set to FALSE
#' @param patient_to_exclude name of patient to be removed for the regression computation (because of missing data)
#' @param dependent_variable name of variables to be used as the dependent variable in the regression model
#' @param independent_variable name of variables to be used as independent variables in the regression model
#' 
#' @param do_WisecondorX If you don't want to perform a WisecondorX analysis, set to FALSE
#' @param path_input_wisecondorx where to find results from wisecondorX (see example folder to see how to organize the data)
#' @param path_annotsv where to find results from annotsv (see example folder to see how to organize the data)
#' 
#' @param do_QDNASeq If you don't want to perform a QDNASeq analysis, set to FALSE
#' @param path_input_qdnaseq where to find results from QDNASeq (see example folder to see how to organize the data)
#' 
#' @param do_cni_score_computation If you already have computed your CNI score and data are store in input folder, set to FALSE
#' @param software type of software used, can be "QDNASeq", "WisecondorX" or c("QDNASeq", "WisecondorX")
#' @param control_name name of the group to be used as control for the cni computation
#' @param group_name name of the group on which to compute the cni score
#' 
#' @param do_cni_global If you don't want to perform a CNI score global analysis, set to FALSE
#' @param group_name_global_cni which group of cni score to be compared (TuPo TuPr) 
#' @param path_input_global_cni where to retrieve clinical data and the cni score from different software
#' @param path_analysis_global_cni where to store the new annotated clinical data
#' 
#' @param do_cni_analysis If you don't want to perform a CNI score analysis, set to FALSE
#' @param group_name_analysis which group of cni score to be compared (TuPo TuPr)
#' @param software_CNI which software has been used to compute such cni score
#' @param path_input_cni_analysis where to retrieve the cni score
#' @param clinical_data_cni_analysis where to find the clinical data
#' @param var_of_interest_cni_analysis which columns names of clinical data to be used to analyse the cni score
#' 
#' @param do_facets_analysis If you don't want to perform a facets analysis, set to FALSE
#' @param path_input_facets_analysis where to retrieve facets data
#' @param patient_to_change list of patient with bad name in mutation_highlander file
#' @param patient_new_val list of new patient name to replace in mutation_highlander file
#' 
#' @param do_mutation_analysis If you don't want to perform a mutation analysis, set to FALSE
#' @param mutation_data full path + file name.tsv of the mutation file retrieved from Highlander
#' @param patient_to_remove list of patient to remove from mutation_highlander file
#' 
#' @param do_combined_analysis_mutation_and_facets If you don't want to perform combined analysis mutation and facets, set to FALSE
#' @param path_clean_mutation path to mutation file as outputed by Mutation_analysis() function
#' @param path_indi_gene path to patient_indiv genes file as outputed by Facets_analysis() function
#' @param path_cancer_gene path to file with cancer gene list with mandatory columns Hugo_Symbol and for analysis "oncoGene_and_tumorSupr" also Is_Oncogene Is_Tumor_Suppressor_Gene
#' @param analysis_to_do name of analysis to perform, or Default, or oncoGene_and_tumorSupr or both
#' @param group_to_plot which group of should be analyzed (ASC, TuPo, TuPr) 
#' @param min_occurence_gene minimum number of time a gene is present in combined data (mutation + Facets) to be taken into account in the analysis
#' 
#' 
#' Requiered packages and libraries (make sure they are installed prior running this pipeline)
#' @readxl
#' @writexl
#' @ComplexHeatmap
#' @circlize
#' @survminer
#' @survival
#' @broom
#' @rcompanion
#' @nnet
#' @regclass
#' @tidyverse
#' @IRanges
#' @biomaRt
#' @RColorBrewer
#' @zoo
#' @plotly
#' @htmlwidgets
#' @stats
#' @ggpubr
#' @openxlsx
#' @cowplot
#' @clusterProfiler
#' @org.Hs.eg.db
#' @GOSemSim
#' @ReactomePA
#' @DOSE
#' 
#' 
#' The pipeline will create a folder within the results folder called YEAR_MONTH_DAY_HOUR_MINUTE_SECOND_NAME_ANALYSIS 
#' where the YEAR_MONTH_DAY_HOUR_MINUTE_SECOND will be taken from the moment you launch the analysis. 
#' All results will be found in this folder and with one sub folder for each analysis
#' Survival Analysis folder contains the results for the survival analysis 
#'
#'
#' By Simon Boutry 09/12/2022

############################################################################################################################
############################ Change parameters in this section to fit your analysis             ############################ 
############################################################################################################################
#path to pipeline 
path_to_pipeline <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Pipeline/"


#If you wanna reuse an existing output folder to store results please give path here bellow and set old_reuse = TRUE
#by default old_reuse = FALSE, and old_path is not use, a new ouput folder will be created automatically
old_path <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Results_final_06_05_2024/"
name_analysis <- "New_cohort" 
old_reuse <- TRUE


#Pre-treatment of clinical data
do_pretreatment <- TRUE #If you don't want to pre-treat your data, set to FALSE
clinical_annotation_to_remove <- NA #when no clinical annotation to remove from clinical data and clinical data are already in the desired shape
#clinical_annotation_to_remove <- c("CA125_initial", "PCI_initial", "Reponse_biologique", "PCI_post", "reponse_per_op")


#clustering based on output data from Raphael Helaers genomicRegionIntersection java program
do_clustering <- FALSE
input_data <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Clustering_facets/intersection_TuPr.tsv"
clinical_data <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Pipeline/input/New_cohort/clinical_data.xlsx"
nbr_group_cluster <- c(1,2,3)


#Survival Analysis
do_survival <- FALSE #If you don't want to perform a survival analysis, set to FALSE
variable_for_survival <- c("censored_DFS", "DFS") #other valid example c("Censored_OS", "OS")
to_investigate_survival <- c("BRCA", "FIGO_stage", "Chemo_type", "Nombre_de_cures_NACT", "CA125_post", 
                             "NLR", "AINS_per_op", "Diff_PCI", "Anapath_response", "cluster_2_fac", "cluster_3_fac")

#Regression
do_regression <- FALSE #If you don't want to perform a regression analysis, set to FALSE
patient_to_exclude <- c("IRIS", "WIMI", "TRUF", "ZONA", "WILE")
dependent_variable <-  c("DFS", "PFI", "cluster_2_fac", "cluster_3_fac")
restricted_model <- c("CA125_post", "Diff_PCI", "Anapath_response")
full_model <- c("BRCA", "FIGO_stage", "Chemo_type", "Nombre_de_cures_NACT", "CA125_post", 
                "NLR", "AINS_per_op", "Diff_PCI", "amount_RD", "Anapath_response")


#WisecondorX data analysis 
do_WisecondorX <- FALSE #If you don't want to perform a WisecondorX analysis, set to FALSE
path_input_wisecondorx <- paste(path_to_pipeline, "input/", name_analysis, "/wisecondor/", sep = "")
#path_annotsv <- paste(path_to_pipeline, "input/", name_analysis, "/Annotsv/", sep = "")


#QDNASeq data analysis 
do_QDNASeq <- FALSE #If you don't want to perform a QDNASeq analysis, set to FALSE
path_input_qdnaseq <- paste(path_to_pipeline, "input/", name_analysis, "/QDNASeq/", sep = "")


#CNI score computation (if CNI score have not been computed before)! This might take time!!!
do_cni_score_computation <- FALSE #If you already have computed your CNI score, set to FALSE
software <- c("WisecondorX") #c("QDNASeq","WisecondorX")
control_name = c("San") #c("controls", "San")
group_name = c("Asc", "TuPo", "TuPr")


#CNI score global analysis
#You need to have computed cni score based on at least two software to use that function!!!
do_cni_global <- FALSE #If you don't want to perform a CNI score global analysis, set to FALSE
group_name_global_cni <- c("TuPr", "TuPo")
path_input_global_cni <- "" #where CNI score are stored, if left empty will check in your results folder
#This function also need the control_name input parameter as used in CNI score computation


#CNI analysis
do_cni_analysis <- FALSE #If you don't want to perform a CNI score analysis, set to FALSE
group_name_analysis <- c("TuPr", "TuPo")
software_CNI <- c("QDNASeq","WisecondorX")#"QDNASeq"
#TODO put the path where to retrieve the cni score (only 1 per software_CNI)
path_input_cni_analysis <- c("C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Results_final/QDNASeq_CNI_score/San/",
                             "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Results_final/WisecondorX_CNI_score/San/")
#TODO put the right path for the new clinical data in your results folder, these clinical data have been created by data_pretreatment() function and will be used by all functions bellow
clinical_data_cni_analysis <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Results_final_06_05_2024/clinical_data_used_for_analysis.xlsx"
var_of_interest_cni_analysis <- c("Anapath_response", "DFS", "PFI")  #c("BRCA", "FIGO_stage", "CA125_post", "Diff_PCI", "Anapath_response", 
# "censored_DFS", "Censored_OS", "DFS", "PFI")


#Facets analysis
do_facets_analysis <- FALSE #If you don't want to perform a facets analysis, set to FALSE
path_input_facets_analysis <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Pipeline/input/New_cohort/facets_data/"
patient_to_change = c("CYLO-TuPrP", "DADA-TuAPr", "DZSC-TuPrP", "NACH-TUPRE", "POWR-TUPO", "KDAU-TUPRE")
patient_new_val = c("CYLO-TuPr", "DADA-TuPr", "DZSC-TuPr", "NACH-TuPr", "POWR-TuPo", "KDAU-TuPr")
#also need clinical data as input argument, we use the clinical_data_cni_analysis above


#Mutations analysis
do_mutation_analysis <- FALSE
mutation_data <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Pipeline/input/New_cohort/mutation/updated_mutations.tsv"
patient_to_remove = c("MYKO-TuPo", "MYKO-TuPr", "IRIS-TuPo", "IRIS-TuPr", "IRIS-ASC", "WIMI-TuPo", "WIMI-TuPr",
                      "TURF-TuPr", "XARA-FIPRE","WILE-SGPRE")

#TODO remove KDAU from analysis
no_mutation <- c("ADES-TuPo", "ATAL-TuPo", "BOBA-TuPo", "FETT-TuPo", "LOUF-TuPo", "RENT-TuPo", "XARA-TuPo")  #c("ADES-TuPo", "FETT-TuPo")
#also need clinical data as input argument, we use the clinical_data_cni_analysis above
#also need patient_to_change and patient_new_val


#Combined analysis mutation and facets
do_combined_analysis_mutation_and_facets <- TRUE
path_clean_mutation = "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Results_final_06_05_2024/Mutations_analysis/mutation_data_cleaned.tsv"
path_sig_mutation <- "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Results_final_06_05_2024/Mutations_analysis/Only_sig_enrich_mutation_data.tsv"
path_indi_gene = "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Results_final_06_05_2024/Facets_analysis/patients_indiv genes.tsv"
path_cancer_gene = "C:/Users/boutrys/OneDrive - UCL/1_POSTDOC/Pipeline/input/Last_input_modified_on_13_01_2023/mutation/cancerGeneList.tsv"
analysis_to_do = c("Default")
group_to_plot = c("ASC", "TuPr", "TuPo")
min_occurence_gene = 2
#also need clinical data as input argument, we use the clinical_data_cni_analysis above


############################################################################################################################
############################              Pipeline variables initialization                       ########################## 
############################################################################################################################
#Load all necessary functions necessary that must be in the folder functions 
path_functions <- paste(path_to_pipeline, "functions/", sep = "") 
function_to_load <- list.files(path_functions)
for (i in 1:length(function_to_load)) {
  source(paste(path_functions, function_to_load[i], sep = ""))
}

if(old_reuse){
  path_output <- old_path
}else{
  #Create a new folder into results folder to store results of analysis. Name of results folder is YEAR_MONTH_DAY_HOUR_MINUTE_SECOND_NAME_ANALYSIS
  tmp_date <- strsplit(as.character(Sys.time()), split = " ")[[1]]
  path_output <- paste(path_to_pipeline, "results/",
                       paste(paste(strsplit(tmp_date[[1]], split = "-")[[1]], collapse = "_"), 
                             paste(strsplit(tmp_date[[2]], split = ":")[[1]], collapse = "_"), 
                             name_analysis, 
                             sep = "_"),
                       "/",
                       sep = "")
  dir.create(path_output)
}


############################################################################################################################
############################              Pre-treatment of clinical data                        ############################ 
############################################################################################################################
if(do_pretreatment){
  path_input <- paste(path_to_pipeline, "input/", name_analysis, "/", sep = "")
  res <- data_pretreatment(full_path = path_input,
                           to_remove = clinical_annotation_to_remove,
                           output = path_output)
  data <- res[[1]]
  data_noFactor <- res[[2]]
}


############################################################################################################################
############################              Clustering Analysis                                   ############################ 
############################################################################################################################
if(do_clustering){
  path_output_clustering <- paste(path_output, "clustering_analysis/", sep = "")
  dir.create(path_output_clustering)
  clustering_analysis(input_data_cluster = input_data,
                      clinical_data = clinical_data,
                      nbr_group_cluster = nbr_group_cluster,
                      output = path_output_clustering)  
}


############################################################################################################################
############################              Survival Analysis                                     ############################ 
############################################################################################################################
if(do_survival){
  path_output_survival <- paste(path_output, "survival_analysis/", sep = "")
  dir.create(path_output_survival)
  result_survival <- survival_analysis(data_survival = data,
                                       variable_for_survival = variable_for_survival,
                                       to_investigate_survival = to_investigate_survival,
                                       output = path_output_survival)  
}


############################################################################################################################
############################              Regression Analysis                                   ############################ 
############################################################################################################################
if(do_regression){
  path_output_regression <- paste(path_output, "regression_analysis/", sep = "")
  dir.create(path_output_regression)
  result_regression <- regression_analysis(data_regression = data,
                                           patient_to_exclude = patient_to_exclude,
                                           dependent_variable = dependent_variable,
                                           restricted_model = restricted_model,
                                           full_model = full_model,
                                           output = path_output_regression) 
}


############################################################################################################################
############################              WisecondorX data analysis                             ############################ 
############################################################################################################################
if(do_WisecondorX){
  path_output_wisecondorX <- paste(path_output, "wisecondorx_analysis/", sep = "")
  dir.create(path_output_wisecondorX)
  result_wisecondorx <- wisecondorx_analysis(path = path_input_wisecondorx,
                                             path_output = path_output_wisecondorX) 
}


############################################################################################################################
############################              QDNASeq data analysis                                 ############################ 
############################################################################################################################
if(do_QDNASeq){
  path_output_qdnaseq <- paste(path_output, "qdnaseq_analysis/", sep = "")
  dir.create(path_output_qdnaseq)
  result_qdnaseq <- qdnaseq_analysis(path_input = path_input_qdnaseq,
                                     path_output = path_output_qdnaseq,
                                     do_segmentation = FALSE) 
}


############################################################################################################################
############################              CNI score computation                                 ############################ 
############################################################################################################################
if(do_cni_score_computation){
  #lauch the computation for each software 
  for (i in 1:length(software)) {
    if(software[i] == "QDNASeq"){
      path_input_cni <- path_input_qdnaseq
    }else if(software[i] == "WisecondorX"){
      path_input_cni <- path_input_wisecondorx
    }else{
      print("softwar must be QDNASeq or WisecondorX")
      stop()
    }
    path_output_CNI_score <- paste(path_output, software[i], "_CNI_score/", sep = "")
    dir.create(path_output_CNI_score)
    result_CNI_score <- CNI_score_computation(path_input = path_input_cni,
                                              software = software[i],
                                              control_name = control_name,
                                              group_name = group_name,
                                              path_output = path_output_CNI_score)
  } 
}


############################################################################################################################
############################              CNI score global analysis                             ############################ 
############################################################################################################################
if(do_cni_global){
  path_output_global_cni <- paste(path_output, "CNI_global_analysis/", sep = "")
  if(path_input_global_cni == ""){
    path_input_global_cni <- path_output
  }
  CNI_global_analysis(group_name = group_name_global_cni,
                      control_group_name = control_name,
                      path_input = path_input_global_cni,
                      path_output = path_output_global_cni,
                      path_analysis = path_output) 
}


############################################################################################################################
############################              CNI score analysis                                    ############################ 
############################################################################################################################
if(do_cni_analysis){
  for (i in 1:length(software_CNI)) {
    path_output_cni_analysis <- paste(path_output, software_CNI[i], "_CNI_score_analysis/", sep = "")
    CNI_analysis(group_name = group_name_analysis,
                 software = software_CNI[i],
                 path_input = path_input_cni_analysis[i],
                 path_output = path_output_cni_analysis,
                 clinical_data = clinical_data_cni_analysis,
                 var_of_interest = var_of_interest_cni_analysis)
  }
}


############################################################################################################################
############################              FACETS  analysis                                      ############################ 
############################################################################################################################
if(do_facets_analysis){
  path_output_facets_analysis <- paste(path_output, "Facets_analysis/", sep = "")
  Facets_analysis(path = path_input_facets_analysis,
                  to_change = patient_to_change,
                  new_val = patient_new_val,
                  clinical_data_path = clinical_data_cni_analysis,
                  path_output = path_output_facets_analysis)  
}


############################################################################################################################
############################              Mutations  analysis                                   ############################ 
############################################################################################################################
if(do_mutation_analysis){
  path_output_mutation_analysis <- paste(path_output, "Mutations_analysis/", sep = "")
  Mutation_analysis(path_mutation_data = mutation_data,
                    no_mutation = no_mutation,
                    clinical_data_path = clinical_data_cni_analysis,
                    to_change = patient_to_change,
                    new_val = patient_new_val,
                    to_remove = patient_to_remove,
                    path_output = path_output_mutation_analysis)
}



############################################################################################################################
############################              Combined analysis Mutation and Facets                 ############################ 
############################################################################################################################
if(do_combined_analysis_mutation_and_facets){
  path_output_combined_analysis_mutation_and_facets <- paste(path_output, "combined_analysis/", sep = "")
  combined_analysis_mutation_and_facets(path_clean_mutation = path_clean_mutation,
                                        path_sig_mutation = path_sig_mutation,
                                        path_indi_gene = path_indi_gene,
                                        clinical_data_path = clinical_data_cni_analysis,
                                        path_cancer_gene = path_cancer_gene,
                                        analysis_to_do = analysis_to_do,
                                        group_to_plot = group_to_plot,
                                        min_occurence_gene = min_occurence_gene,
                                        path_output = path_output_combined_analysis_mutation_and_facets)
}


