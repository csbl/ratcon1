# This script is involved in generating biomarker predictions 
# based on gene expression changes using the new algorithm,
# TIMBR (Transcriptionally-Inferred Biomarker Response).

# This script calculates gene expression changes using
# gene expression microarrays of rat and human hepatocytes treated
# with various pharmaceutical compounds and environmental toxicants.
# After preprocessing raw gene expression microarray data,
# we can perform differential expression analysis 
# to determine significant/insignificant gene expression changes. 
# Gene expression changes can be summarized into 
# TIMBR reaction weights based with the following outputs:
# gene_id: Entrez gene identifier (integer)
# logfc: log2 fold changes represent average differences 
# in log2-transformed gene expression values between treatment and control samples
# fdr: false discovery-rate adjusted q-values (FDR) represent 
# statistical significance of gene expression changes. 
# For gene expression changes, fdr < 0.1 was considered
# significantly differentially expressed. 

# Most data/code needed to determine gene expression changes are available 
# as Supplementary Tables/Datasets or on the ratcon GitHub website:
# www.github.com/edikblais/ratcon

# Raw gene expression microarray data can be obtained from Open TG-GATEs:
# Open Toxicogenomics Project-Genomics Assisted Toxicity Evaluation System
# http://www.ncbi.nlm.nih.gov/pubmed/25313160.
# http://toxico.nibio.go.jp/english/index.html
# Raw gene expression microarray data are also available from ArrayExpress:
# E-MTAB-797 - Transcription profiling by array of rat hepatocytes 
# treated with approximately 130 chemicals in vitro
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-797/
# E-MTAB-798 - Transcription profiling by array of human hepatocytes 
# treated with approximately 130 chemicals in vitro
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-798/

# TIMBR predictions can be reproduced by running the following scripts:
# 1. ncomm_blais_timbr_expression.R
#   Inputs raw expression data
#   Outputs gene expression changes
#   Helper functions: ncomm_helper.R
# 2. ncomm_blais_timbr_weights.R
#   Inputs gene expression changes, rat and human GPR rules
#   Outputs TIMBR reaction weights
#   Helper functions: ncomm_helper.R
# 3. ncomm_blais_timbr_predictions.m
#   Inputs TIMBR reaction weights, rat and human metabolic networks
#   Outputs raw TIMBR predictions
#   Helper functions: timbr.m
# 4. ncomm_blais_timbr_analysis.R
#   Inputs raw TIMBR predictions
#   Outputs normalized TIMBR production scores, manuscript figures
#   Helper functions: ncomm_helper.R

options(java.parameters = "-Xmx4g")
options(stringsAsFactors = FALSE)

source("ncomm_helper.R")

library(reshape2)
library(plyr)
library(dplyr)
library(oligo)
library(genefilter)
library(rat2302.db)
library(pd.rat230.2)
library(hgu133plus2.db)
library(pd.hg.u133.plus.2)
library(limma)

# Make sure dplyr functions are used
select <- dplyr::select
rename <- dplyr::rename
desc <- dplyr::desc
collapse <- dplyr::collapse
combine <- dplyr::combine
intersect <- dplyr::intersect
setdiff <- dplyr::setdiff
union <- dplyr::union
summarize <- dplyr::summarize
slice <- dplyr::slice


# Specify path to R/Biocondcutor working directory:
path.rcode.directory = "C:/Users/edikb/Dropbox/2016_01_ncomm_submission/"
# Specify path to downlaoded Supplementary Tables/Datasets
path.manuscript.directory = "C:/Users/edikb/Dropbox/2016_01_ncomm_submission/"
# Specify path to raw Affymetrix .CEL files
path.affymetrix.directory = "D:/tggates/"
# Specify directory for saving normalized gene expression datasets
path.eset.directory = "C:/Users/edikb/Desktop/tox/"
# Specify directory for saving gene expression changes
path.efit.directory = "C:/Users/edikb/Google Drive/results/tox/"

# Load phenotype information for Affymetrix samples.
sample.info = read.table(paste0(path.manuscript.directory, "ncomm_blais_affymetrix_sample_info.txt.gz"), 
                         sep = "\t", quote = "", header = T,
                         stringsAsFactors = F, check.names = F, fill = T) %>%
  mutate(sample_file_full = paste0(path.affymetrix.directory,sample_id),
         sample_file_exists = sample_file_full %>% file.exists) %>%
  mutate(eset_dose_strategy = c("one" = "s.d.","rep" = "r.d.")[eset_strategy]) %>%
  mutate(eset_species = c("hsa" = "human","rno" = "rat")[eset_organism]) %>%
  mutate(eset_name = paste0(eset_drug_name, " (",eset_drug_abbrev," ",
                            eset_dose_range, " ", eset_dose_strategy,")",
                            " @ ", eset_time_name," in ", eset_species, " ", eset_part)) %>%
  mutate(dataset_id = eset_id, dataset_name = eset_name) %>% ef_df

# sample.info includes sample information for microarray samples obtained from ArrayExpress:

# E-MTAB-797 - Transcription profiling by array of rat hepatocytes 
# treated with approximately 130 chemicals in vitro
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-797/
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-798/E-MTAB-798.sdrf.txt

# E-MTAB-798 - Transcription profiling by array of human hepatocytes 
# treated with approximately 130 chemicals in vitro
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-798/
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-798/E-MTAB-798.sdrf.txt

# Sample information from E-MTAB-799 and E-MTAB-800 
# were also included but not analyzed in this study

# sample.info also includes download links to 
# compressed Affymetrix .CEL files from ArrayExpress.

# The abbrevation, 'eset', stands for Expression Set and typically refers 
# to an individual dataset preprocessed using the oligo package
# The abbrevation, 'efit', stands for Expression Fit and typically refers 
# to an individual comparison within a specific dataset analyzed using limma package

# Verify that .CEL files exist in the right location before preprocessing datasets
sample.info %>% count(sample_id,eset_id) %>% ungroup %>% arrange(desc(n))
sample.info %>% count(eset_id,sample_file_exists)
sample.info %>% count(sample_file_exists)
sample.info %>% filter(!sample_file_exists) %>% count(eset_id)

# Due to a large number of samples, we preprocessed 
# gene expression data separately as individual datasets 
# for each organism/drug/timepoint. As a result, each preprocessed 
# expression set (eset) was considered an independent dataset 
# with unique control samples. 
dataset.info = sample.info %>% select(starts_with("eset_")) %>% distinct %>%
  mutate(eset_version = "gene", eset_folder = path.eset.directory, 
         eset_type = "core", eset_file_ext = ".txt.gz") %>%
  ef_dataset_file_check()

dataset.design = sample.info %>% select(eset_id) %>% 
  distinct %>% with(setNames(eset_id, eset_id)) %>% 
  lapply(function(x) (~ dose_id))
dataset.contrasts = sample.info %>% 
  select(eset_id, dose_id) %>% distinct %>% 
  filter(dose_id != "d0") %>% 
  ef_df_slice("eset_id") %>% 
  lapply(function(x) setNames(x$dose_id, x$dose_id))

# We preprocessed raw gene expression data using the Robust Multichip Average (RMA) 
# method implemented in the oligo package. RMA involves multiple normalized steps:
# background subtraction, quantile normalization, and summarization via median-polish.

# The function, ef_dataset_preprocess, saves each individual dataset as 3 compressed text files:
# normalized values stored as an expression matrix (*_exprs.txt.gz)
# gene information stored as feature data (*_fdata.txt.gz) 
# sample information stored as phenotype data (*_pdata.txt.gz)
# The filename for each dataset is based on the following information separated by "_": 
#   organism (hsa = human; rno = rat)
#   cell type (hep = hepatocytes)
#   treatment duration (t1 = 2 hours; t2 = 8 hours; t3 = 24 hours)
#   treatment dosing strategy (one = single dose; rep = repeated dose)
#   treatment compound (e.g. caffeine; acetaminophen)
#   feature type (gene)
#   genefilter method (eset = featureFilter; nset = nsFilter; core = no filtering)

dataset.status = dataset.info %>% #filter(!eset_ok) %>% 
  with(eset_id) %>% unique %>% 
  lapply(ef_dataset_preprocess,sample.info,x.version = "gene",x.folder = path.eset.directory) %>%#
  bind_rows

## Calculate gene expression changes
# We performed differential expression analysis 
# using Linear Models for Microarrays (limma package) 
# after preprocessing gene expression datasets.
# Independently for each dataset, we determined gene expression changes 
# by comparing control samples (d0) to samples treated 
# with either a low dose (d1), medium dose (d2) or high dose (d3). 
# Concentrations for each dose varied by compound and by organism.

# The function, ef_limma_preprocess, performs differential expression analysis
# and saves each set of gene expression changes as a separate text file.
# The filename is based on the following information separated by "_": 
#   organism (hsa = human; rno = rat)
#   cell type (hep = hepatocytes)
#   treatment duration (t1 = 2 hours; t2 = 8 hours; t3 = 24 hours)
#   treatment dosing strategy (one = single dose; rep = repeated dose)
#   treatment compound (e.g. caffeine; acetaminophen)
#   treatment dose (d1 = low; d2 = medium; d3 = high)
#   feature type (gene)
#   genefilter method (eset = featureFilter; nset = nsFilter; core = no filtering)
#   limma method (robust = robust linear fit; ls = least squares fit)

library(limma)

limma.setup = ef_limma_setup(dataset.contrasts, dataset.design, sample.info, 
                             path.eset.directory, path.efit.directory, "eset", "gene")
limma.info = limma.setup %>% 
  select(organism_id, dataset_id, contrast_id, contrast_string, contrast_abbrev,
         starts_with("eset_"), starts_with("efit_")) %>% distinct

limma.info %>% select(eset_id, eset_ok, eset_type, eset_version) %>% distinct %>% count(eset_ok)
limma.info %>% select(efit_id, efit_ok, efit_type, efit_version, efit_lm) %>% distinct %>% count(efit_ok)

limma.status = limma.info %>% #filter(!efit_ok) %>% 
  filter(eset_ok) %>% arrange(eset_root, efit_root) %>%
  ef_df_slice("eset_root") %>% 
  lapply(ef_limma_preprocess) %>% 
  bind_rows

# Save limma info so gene expression changes can be 
# used to generate TIMBR reaction weights in the next script:
# ncomm_blais_timbr_weights.R
limma.info %>% select(-ends_with("_ok"), -ends_with("_exists"), -ends_with("_folder"), -ends_with("_file"),
                      -ends_with("_exprs"), -ends_with("_fdata"), -ends_with("_pdata"),
                      -one_of(c("eset_file_efits"))) %>% 
  write.table(paste0(path.efit.directory, "ncomm_blais_limma_info.txt.gz"), sep = "\t", quote = F, row.names = F)

# End
