# Ratcon1 database
Code availability for the manuscript titled "Reconciled rat and human metabolic networks for comparative toxicogenomics analyses and biomarker predictions" by Blais *et al*.

The article is now published and available at http://www.nature.com/articles/ncomms14250

## Simulating metabolic tasks with *iRno* and *iHsa*

### blais_check_tasks.m

This MATLAB script imports *iRno* and *iHsa* into the RAVEN toolbox and performs the checkTasks function for 327 metabolic tasks described in Supplementary Table 4. Most tasks capture metabolic functions that are common to both rats and humans. Some tasks are meant to fail in both organisms while others are meant to fail in only one organism, consistent with known species-specific differences between rat and human metabolism.

## Generating a draft of *iRno* based on *iHsa*

### blais_gpr_conversion.R

This script was used to generate Fig. 2a and Supplementary Fig. 3. For Supplementary Fig. 3, this R script inputs orthology annotations from multiple databases by reading in data from Supplementary Table 2 to generate rat gene-protein-reaction (GPR) relationship rules based on human GPR rules from HMR2 (after replacing Ensembl Gene Identifiers with Entrez gene identifiers). For Fig. 2a, this R script compares the complexity of gene-protein-reaction (GPR) directly between rat and human metabolic models by reading in data from Supplementary Table 3. This script inputs gene-protein-reaction (GPR) relationship information from a superset of reactions included in *iRno* and *iHsa* to compare the global distribution of GPR sizes between the rat and human models. The goal of this script was to demonstrate how GPR sizes varied between species when filtering orthology annotations for the automated conversion of *iHsa* to *iRno* as well as after manual curation. 

## Predicting biomarkers based on gene expression changes with __TIMBR__

### blais_timbr_expression.R

  Inputs raw gene expression microarray data
  
  Outputs gene expression changes
  
  Helper functions: ncomm_helper.R

### blais_timbr_weights.R

  Inputs gene expression changes, rat and human GPR rules
  
  Outputs TIMBR reaction weights
  
  Helper functions: ncomm_helper.R

### blais_timbr_predictions.m

  Inputs TIMBR reaction weights, rat and human metabolic networks
  
  Outputs raw TIMBR predictions
  
  Helper functions: timbr.m and ncomm_blais_model2irrev.m

### blais_timbr_analysis.R

  Inputs raw TIMBR predictions
  
  Outputs normalized TIMBR production scores, manuscript figures
  
  Helper functions: ncomm_helper.R

## Helper function files:

### ncomm_helper.R
  R source file that defines several helper functions used in the R scripts above.

### timbr.m
  MATLAB implementation of the TIMBR (transcriptionally-inferred metabolic biomarker response) algorithm

### ncomm_blais_xls2model.m
  Modified version of xls2model from the COBRA toolbox (www.github.com/opencobra/cobratoolbox) 
  that fixes a minor bug associated with importing Excel-formatted versions 
  of rat and human metabolic networks in the ratcon database

### ncomm_blais_model2irrev.m
  Modified version of convertToIrreversible from the COBRA toolbox (www.github.com/opencobra/cobratoolbox) 
  that facilitates the mapping of TIMBR reaction weights to irreversible reactions

