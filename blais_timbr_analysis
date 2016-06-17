# This script is involved in generating biomarker predictions 
# based on gene expression changes using the new algorithm,
# TIMBR (Transcriptionally-Inferred Biomarker Response).
# 
# This script anaylzes raw TIMBR predictions generated using
# rat and human metabolic networks and gene expression changes 
# from rat and human hepatocytes treated with various 
# pharmaceutical compounds and environmental toxicants.

options(java.parameters = "-Xmx4g")
options(stringsAsFactors = FALSE)
source("ncomm_helper.R")

library(reshape2)
library(dplyr)
library(ggplot2)
library(xlsx)

# Specify path to downlaoded Supplementary Tables/Datasets
path.manuscript.directory = "C:/Users/edikb/Dropbox/2016_01_ncomm_submission/"
# Specify path to saved normalized gene expression datasets
path.eset.directory = "C:/Users/edikb/Desktop/tox/"
# Specify path to saved gene expression changes
path.efit.directory = "C:/Users/edikb/Google Drive/results/tox/"
# Specify path to saved raw TIMBR predictions
path.weights.directory = ""
path.timbr.directory = "C:/Users/edikb/Google Drive/results/ncomm_submission/timbr_predictions/"
# path.timbr.directory = "C:/Users/edikb/Google Drive/results/ncomm_submission/timbr_directory/"

rxn.info.load = read.xlsx2(paste0(path.manuscript.directory,"ncomm_blais_supplementary_table3.xlsx"),
                           sheetIndex = 1, startRow = 2) %>% as.tbl


rxn.info = rxn.info.load %>% 
  mutate(enabled = ifelse(!is.na(enabled), enabled == "true", F)) %>% 
  mutate(enabled_hsa = ifelse(!is.na(enabled_hsa), enabled_hsa == "true", F)) %>% 
  mutate(enabled_rno = ifelse(!is.na(enabled_rno), enabled_rno == "true", F))%>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))  %>%
  mutate(met = gsub(" exchange| demand","",rxn_name)) %>%
  mutate(met = ifelse(rxn_id == "RCR30122","vitamin C",met))

limma.info = paste0(path.efit.directory, "ncomm_blais_limma_info.txt.gz") %>% 
  read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F) %>% ef_df %>% 
  filter(eset_part == "hepatocyte") %>% 
  mutate(eset_folder = path.eset.directory) %>% 
  ef_dataset_file_check %>% 
  mutate(efit_folder = path.efit.directory) %>% 
  ef_limma_file_check %>% 
  mutate(drug_id = eset_drug_id, time_id = eset_time_id, dose_id = efit_contrast) %>% 
  mutate(limma_id = paste0(organism_id, "_", drug_id, "_", time_id, "_", dose_id),
         limma_ok = efit_ok)

rno.timbr.weights = paste0(path.weights.directory,"ncomm_blais_timbr_weights_rno.txt") %>% 
  read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F) %>% ef_df
hsa.timbr.weights = paste0(path.weights.directory,"ncomm_blais_timbr_weights_hsa.txt") %>% 
  read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F) %>% ef_df


limma.info %>% count(limma_ok)

timbr.info = bind_rows(list(
  limma.info %>% mutate(condition_id = "ctl", timbr_id = paste0(limma_id, "_ctl")),
  limma.info %>% mutate(condition_id = "trt", timbr_id = paste0(limma_id, "_trt")))) %>% 
  mutate(timbr_folder = path.timbr.directory,
         timbr_file = paste0(timbr_folder, "/", "ncomm_blais_timbr_", timbr_id, ".txt"))


##### IN PROGRESS #######
# I recently changed the .txt file formats for TIMBR reaction weights and TIMBR predictions 
# Before TIMBR predictions for trt and ctl conditions were stored in one file
# representing one set of gene expression changes (dataset/contrast). This code was bulky.
# Now, TIMBR predictions for trt and ctl conditions are stored as separate files
# representing one set of TIMBR reaction weights (dataset/contrast/condition). This code is easier to read.

# I still need to update the code that reads in this data used to generate manuscript figures
