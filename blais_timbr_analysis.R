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
path.manuscript.directory = ""
# path.manuscript.directory = "C:/Users/edikb/Dropbox/2016_01_ncomm_submission/"
# path.manuscript.directory = "E:/sync/Dropbox/2016_01_ncomm_submission/"

# Specify path to saved raw TIMBR predictions
path.timbr.directory = ""
# path.timbr.directory = "C:/Users/edikb/Google Drive/results/ncomm_submission/timbr_predictions/"
# path.timbr.directory = "E:/sync/GoogleDrive/results/ncomm_submission/timbr_predictions/"


rxn.info.load = read.xlsx2(paste0(path.manuscript.directory,"ncomm_blais_supplementary_table3.xlsx"),
                           sheetIndex = 1, startRow = 2) %>% as.tbl

timbr.predictions.load = read.xlsx2(paste0(path.manuscript.directory,"ncomm_blais_supplementary_table7.xlsx"),
                                    sheetIndex = 1, startRow = 2) %>% as.tbl

expression.summary.load = read.xlsx2(paste0(path.manuscript.directory,"ncomm_blais_supplementary_table6.xlsx"),
                                     sheetIndex = 1, startRow = 2) %>% as.tbl

rxn.info = rxn.info.load %>% 
  mutate(enabled = ifelse(!is.na(enabled), enabled == "true", F)) %>% 
  mutate(enabled_hsa = ifelse(!is.na(enabled_hsa), enabled_hsa == "true", F)) %>% 
  mutate(enabled_rno = ifelse(!is.na(enabled_rno), enabled_rno == "true", F))%>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))  %>%
  mutate(met = gsub(" exchange| demand","",trimws(rxn_name))) %>%
  mutate(met = ifelse(rxn_id == "RCR30122","vitamin C",met))

timbr.predictions = timbr.predictions.load %>%
  mutate(production_score = as.numeric(as.character(production_score))) %>% 
  left_join(rxn.info %>% select(rxn_id, rxn_name, met))
expression.summary = expression.summary.load %>% 
  mutate(drug_selected = ifelse(!is.na(drug_selected), grepl("true", tolower(drug_selected)), NA)) %>%
  mutate(dose_selected = ifelse(!is.na(dose_selected), grepl("true", tolower(dose_selected)), NA))

timbr.compare = timbr.predictions %>% 
  semi_join(expression.summary %>% filter(drug_selected, dose_selected) %>% 
              select(organism_id, drug_id, time_id, dose_id) %>% distinct) %>%
  # left_join(rxn.info %>% select(rxn_id, rxn_name, met)) %>% 
  dcast(drug_id + rxn_id + rxn_name + met ~ organism_id, 
        value.var = "production_score", fill = NA) %>% as.tbl %>% 
  arrange(rno) %>% group_by(drug_id) %>% mutate(rno_rank = 1:n(), rno_pctile = 100 * rno_rank / n()) %>% ungroup %>%
  arrange(hsa) %>% group_by(drug_id) %>% mutate(hsa_rank = 1:n(), hsa_pctile = 100 * hsa_rank / n()) %>% ungroup %>%
  mutate(hsa_qrtile = ifelse(hsa_pctile > 75, 1, ifelse(hsa_pctile < 25, -1, 0))) %>% 
  mutate(rno_qrtile = ifelse(rno_pctile > 75, 1, ifelse(rno_pctile < 25, -1, 0))) %>% 
  mutate(hsa_border = ifelse(hsa_pctile > 75 | hsa_pctile < 25, "#000000", ifelse(hsa_pctile < 25, "#000000", NA))) %>% 
  mutate(rno_border = ifelse(rno_pctile > 75 | rno_pctile < 25, "#000000", ifelse(rno_pctile < 25, "#000000", NA))) %>% 
  left_join(expression.summary %>% select(drug_id, drug_name) %>% mutate(drug_abbrev = drug_id) %>% distinct)

timbr.drug.summary = timbr.compare %>%
  group_by(drug_id, drug_name, drug_abbrev) %>% 
  summarize(drug_cor = cor.test(hsa, rno, eps = 1e-200)$estimate,
            drug_pval = cor.test(hsa, rno, eps = 1e-200)$p.value,
            drug_sd_rno = sd(rno),
            drug_sd_hsa = sd(hsa),
            drug_ave_rno = mean(rno),
            drug_ave_hsa = mean(hsa)) %>% ungroup %>% 
  mutate(drug_ave = (drug_ave_rno + drug_ave_hsa) / 2) %>% 
  mutate(drug_cor = ifelse(!is.na(drug_cor), drug_cor, 0)) %>% 
  mutate(drug_pval = ifelse(!is.na(drug_pval), drug_pval, 1)) %>% 
  mutate(drug_fdr = p.adjust(drug_pval)) %>% 
  arrange(desc(drug_cor), drug_fdr, drug_pval) %>% 
  mutate(drug_index = 1:n())

timbr.rxn.summary = timbr.compare %>%
  group_by(rxn_id, rxn_name, met) %>% 
  summarize(rxn_cor = cor.test(hsa, rno, eps = 1e-200)$estimate,
            rxn_pval = cor.test(hsa, rno, eps = 1e-200)$p.value,
            rxn_ave_rno = mean(rno),
            rxn_ave_hsa = mean(hsa)) %>% ungroup %>% 
  mutate(rxn_ave = (rxn_ave_rno + rxn_ave_hsa) / 2) %>% 
  mutate(rxn_cor = ifelse(!is.na(rxn_cor), rxn_cor, 0)) %>% 
  mutate(rxn_pval = ifelse(!is.na(rxn_pval), rxn_pval, 1)) %>% 
  mutate(rxn_fdr = p.adjust(rxn_pval)) %>% 
  arrange(desc(rxn_cor), rxn_fdr, rxn_pval) %>% 
  mutate(rxn_index = 1:n())

antipyretic.validation.mets = c(
  "prostaglandin E2")
caffeine.validation.mets = c(
  "ornithine", "arginine", # FP
  "urea", "citrulline", "aspartate", "glutamate",# TP
  "glutamine", "L-lactate", "pyruvate", "acetoacetate", "(R)-3-hydroxybutanoate") # TN
theophylline.validation.mets = c(
  "urea", "citrulline", "glutamate", "glucose", "urate")
extra.mets = c(
  "albumin", "ribose", "uracil", 
  "nicotinate", "inosine", "orotate", 
  "palmitate","oleate", "stearate", "4-aminobutyrate", 
  "cholesterol", "acetylcholine", 
  "cholate", "taurocholate", "glycocholate", 
  "chenodeoxycholate")# "taurochenodeoxycholate", "glycochenodeoxycholate"
select.common.mets = c(
  "GSH", "acetate",
  "creatinine", "creatine",
  "glycine", "alanine")
conditionally.essential.amino.acids = c(
  "proline", "glycine", "cysteine", "glutamine", "tyrosine", "arginine")
non.essential.amino.acids = c(
  "aspartate", "glutamate", "serine", "asparagine", "alanine")

original.mets = c(
  "urea","ornithine","citrulline","creatinine","creatine", "urate", 
  "glycine","aspartate","arginine","glutamine","glutamate","alanine",
  "glucose","pyruvate","L-lactate","acetate", "(R)-3-hydroxybutanoate",
  "GSH", "prostaglandin E2") %>% unique


original.drugs = c(
  "acetaminophen", "ibuprofen","aspirin",
  "diclofenac", "fluphenazine", "indomethacin", 
  "mefenamicacid", "naproxen", "nimesulide", "sulindac",
  "phenylbutazone", "benziodarone", "benzbromarone", 
  "colchicine", "caffeine", "theophylline") %>% unique

simplified.mets = c(
  #select.common.mets,
  caffeine.validation.mets,
  theophylline.validation.mets,
  antipyretic.validation.mets) %>% unique

simplified.drugs = c(
  "acetaminophen", "ibuprofen", "aspirin",
  "phenylbutazone", "colchicine", "benziodarone", "benzbromarone", 
  "caffeine", "theophylline") %>% unique


expanded.mets = c(
  extra.mets, 
  conditionally.essential.amino.acids,
  non.essential.amino.acids,
  select.common.mets,
  caffeine.validation.mets,
  theophylline.validation.mets,
  antipyretic.validation.mets) %>% unique

expanded.drugs = timbr.drug.summary %>% filter(drug_sd_rno > 0.7, drug_sd_hsa > 0.7) %>%
  with(drug_id) %>% c(original.drugs) %>% unique

heatmap.original = ef_timbr_heatmap(original.drugs, original.mets)# + theme(aspect.ratio = 1)
ggsave("ncomm_blais_bonus_original_fig6a.pdf", heatmap.original, width = 9.5, height = 11)

heatmap.drugs = original.drugs %>% unique
heatmap.mets = original.mets %>% c("acetoacetate") %>% setdiff(c("acetate")) %>% unique

heatmap.aspect.ratio = length(heatmap.mets) / length(heatmap.drugs)

heatmap.final = ef_timbr_heatmap(heatmap.drugs, heatmap.mets) + 
  theme(aspect.ratio = heatmap.aspect.ratio)
heatmap.final

ggsave("ncomm_blais_final_fig6a.pdf", heatmap.final, width = 8, height = 10)

heatmap.simplified = ef_timbr_heatmap(simplified.drugs, simplified.mets) + 
  theme(aspect.ratio = length(simplified.mets) / length(simplified.drugs))
heatmap.simplified

ggsave("ncomm_blais_bonus_simplified_fig6a.pdf", heatmap.simplified, width = 6, height = 9)

heatmap.expanded = ef_timbr_heatmap(expanded.drugs, expanded.mets, base.size = 24) + 
  theme(aspect.ratio = length(expanded.mets) / length(expanded.drugs))
heatmap.expanded

ggsave("ncomm_blais_bonus_expanded_fig6a.pdf", heatmap.expanded, width = 24, height = 18)

axis.alpha = 0.3
axis.size = 1

timbr.pge2.plot = timbr.compare %>% filter(met %in% c("prostaglandin E2")) %>% 
  mutate(label = ifelse(drug_id %in% c("acetaminophen", "ibuprofen"), drug_name, "")) %>% 
  mutate(color = ifelse(drug_id %in% c("acetaminophen", "ibuprofen"), ggcolor["opposite"], "#BBBBBB")) %>% 
  mutate(shape = 16) %>% # filled circle no outline
  mutate(size = ifelse(label != "", 9, 5)) %>% 
  mutate(alpha = ifelse(label != "", .7, .3)) %>% 
  mutate(hjust = ifelse(hsa > 0, 1, -0.05)) %>% 
  mutate(vjust = 1) %>% arrange(size) %>%
  ggplot(aes(x = hsa, y = rno, label = label, color = color)) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_vline(xintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_point(aes(shape = shape, size = size, alpha = alpha)) + 
  scale_shape_identity() + scale_alpha_identity() + 
  scale_size_identity() + scale_color_identity() + 
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  theme_minimal(base_size = 24) + 
  theme(legend.position = "none") + xlab("") + ylab("")
timbr.pge2.plot
ggsave("ncomm_blais_final_fig6b.pdf", timbr.pge2.plot, width = 6, height = 6)
ggsave("ncomm_blais_final_fig6b_labeled.pdf", timbr.pge2.plot + 
         geom_text(vjust = 1, aes(hjust = ifelse(hsa > 1, 1, -0.05)), size = 6) , width = 6, height = 6)

# Correlation coefficient for Figure 6b
timbr.pge2.cor = timbr.compare %>% 
  filter(met %in% c("prostaglandin E2")) %>% group_by(met) %>% 
  summarize(biomarker_cor = cor.test(hsa, rno, eps = 1e-200)$estimate,
            biomarker_pval = cor.test(hsa, rno, eps = 1e-200)$p.value) %>% ungroup
timbr.pge2.cor

timbr.urate.plot = timbr.compare %>% filter(met %in% c("urate")) %>% 
  mutate(label = ifelse(drug_id %in% c(
    "theophylline", "caffeine", "benzbromarone", 
    "benziodarone", "colchicine", "phenylbutazone"), drug_name, "")) %>% 
  mutate(color = ifelse(drug_id %in% c(
      "benzbromarone","benziodarone","colchicine", "phenylbutazone", "caffeine"), ggcolor["opposite"], 
      ifelse(drug_id %in% c("theophylline"), ggcolor["similar"],"#BBBBBB"))) %>%
  mutate(shape = 16) %>% # filled circle no outline
  mutate(size = ifelse(label != "", 9, 5)) %>% 
  mutate(alpha = ifelse(label != "", .7, .3)) %>% 
  mutate(hjust = ifelse(hsa > 0, 1, -0.05)) %>% 
  mutate(vjust = 1) %>% arrange(size) %>%
  ggplot(aes(x = hsa, y = rno, label = label, color = color)) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_vline(xintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_point(aes(shape = shape, size = size, alpha = alpha)) + 
  scale_shape_identity() + scale_alpha_identity() + 
  scale_size_identity() + scale_color_identity() + 
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  theme_minimal(base_size = 24) + 
  theme(legend.position = "none") + xlab("") + ylab("")
timbr.urate.plot
ggsave("ncomm_blais_final_fig6c.pdf", timbr.urate.plot, width = 6, height = 6)
ggsave("ncomm_blais_final_fig6c_labeled.pdf", timbr.urate.plot + 
         geom_text(vjust = 1, aes(hjust = ifelse(hsa > 1, 1, -0.05)), size = 6), width = 6, height = 6)

# Correlation coefficient for Figure 6c
timbr.urate.cor = timbr.compare %>% 
  filter(met %in% c("urate")) %>% group_by(met) %>% 
  summarize(biomarker_cor = cor.test(hsa, rno, eps = 1e-200)$estimate,
            biomarker_pval = cor.test(hsa, rno, eps = 1e-200)$p.value) %>% ungroup
timbr.urate.cor

# Join caffeine and theophylline TIMBR production scores together side-by-side
timbr.xanthine.derivative.data = timbr.predictions %>% 
  filter(drug_id %in% c("caffeine","theophylline")) %>% 
  mutate(drug = ifelse(drug_id == "caffeine","caff","theo")) %>% 
  dcast(organism_id + rxn_id + rxn_name + met ~ drug, 
        value.var = "production_score", fill = NA) %>% ef_df %>%
  mutate(label = ifelse(met %in% c("urea", "urate", "citrulline","glutamate", "glucose"),met,"")) %>%
  mutate(color = ifelse(label != "", ggcolor[organism_id], "#BBBBBB")) %>% 
  mutate(shape = 16) %>% # filled circle no outline
  mutate(size = ifelse(label != "", 9, 5)) %>% 
  mutate(alpha = ifelse(label != "", .7, .3)) %>% 
  mutate(hjust = ifelse(caff > 0, 1, -0.05)) %>% 
  mutate(vjust = 1) %>% arrange(size)

timbr.xanthine.derivative.plot.rno = timbr.xanthine.derivative.data %>% 
  filter(organism_id == "rno")   %>%
  mutate(theo = pmin(pmax(theo, -2.5), 2.5)) %>%
  mutate(caff = pmin(pmax(caff, -2.5), 2.5)) %>%
  ggplot(aes(x = caff, y = theo, label = label, size = size, color = color)) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_vline(xintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_point(aes(shape = shape, alpha = alpha)) + 
  scale_shape_identity() + scale_alpha_identity() + 
  scale_size_identity() + scale_color_identity() + 
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  theme_minimal(base_size = 24) + 
  theme(legend.position = "none") + xlab("") + ylab("")

timbr.xanthine.derivative.plot.rno
ggsave(paste0("ncomm_blais_final_fig6d.pdf"), timbr.xanthine.derivative.plot.rno, width = 6, height = 6)
ggsave(paste0("ncomm_blais_final_fig6d_labeled.pdf"), timbr.xanthine.derivative.plot.rno + 
         geom_text(aes(label = label, hjust = hjust, vjust = vjust), size = 6, color = "black"), width = 6, height = 6)

timbr.xanthine.derivative.plot.hsa = timbr.xanthine.derivative.data %>% 
  filter(organism_id == "hsa")  %>% 
  mutate(theo = pmin(pmax(theo, -2.5), 2.5)) %>%
  mutate(caff = pmin(pmax(caff, -2.5), 2.5)) %>%
  ggplot(aes(x = caff, y = theo, label = label, size = size, color = color)) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_vline(xintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_point(aes(shape = shape, alpha = alpha)) + 
  scale_shape_identity() + scale_alpha_identity() + 
  scale_size_identity() + scale_color_identity() + 
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  theme_minimal(base_size = 24) + 
  theme(legend.position = "none") + xlab("") + ylab("")

timbr.xanthine.derivative.plot.hsa
ggsave(paste0("ncomm_blais_final_fig6e.pdf"), timbr.xanthine.derivative.plot.hsa, width = 6, height = 6)
ggsave(paste0("ncomm_blais_final_fig6e_labeled.pdf"), timbr.xanthine.derivative.plot.hsa + 
         geom_text(aes(label = label, hjust = hjust, vjust = vjust), size = 6, color = "black"), width = 6, height = 6)

# Correlation coefficients for Figure 6d-e
timbr.xanthine.derivative.cor = timbr.xanthine.derivative.data %>% group_by(organism_id) %>% 
  summarize(organism_cor = cor.test(theo, caff, method = "pearson", eps = 1e-200)$estimate,
            organism_pval = cor.test(theo, caff, method = "pearson", eps = 1e-200)$p.value) %>% ungroup
timbr.xanthine.derivative.cor


# Reproduce scatterplot comparing rat and human TIMBR production scores in response to caffeine 
timbr.caffeine.plot = timbr.compare %>% filter(drug_id %in% c("caffeine")) %>%
  mutate(label = ifelse(met %in% c("urea", "urate", "citrulline","glutamate", "glucose"),met,"")) %>%
  mutate(color = ifelse(label != "", ifelse(
    rno > 0, ggcolor["Elevated"], ggcolor["Reduced"]), "#BBBBBB")) %>% 
  mutate(shape = 16) %>% # filled circle no outline
  mutate(size = ifelse(label != "", 9, 5)) %>% 
  mutate(alpha = ifelse(label != "", .7, .3)) %>% 
  mutate(hjust = ifelse(hsa > 0, 1, -0.05)) %>% 
  mutate(vjust = 1) %>% arrange(size) %>%
  mutate(hsa = pmin(pmax(hsa,-2.5),2.5)) %>% 
  mutate(rno = pmin(pmax(rno,-2.5),2.5)) %>% 
  ggplot(aes(x = hsa, y = rno, label = label, size = size, color = color)) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_vline(xintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_point(aes(shape = shape, alpha = alpha)) + 
  scale_shape_identity() + scale_alpha_identity() + 
  scale_size_identity() + scale_color_identity() + 
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  theme_minimal(base_size = 24) + 
  theme(legend.position = "none") + xlab("") + ylab("")
timbr.caffeine.plot
# Reproduce scatterplot comparing rat and human TIMBR production scores in response to theophylline 
timbr.theophylline.plot = timbr.compare %>% filter(drug_id %in% c("theophylline")) %>%
  mutate(label = ifelse(met %in% c("urea", "urate", "citrulline","glutamate", "glucose"),met,"")) %>%
  mutate(color = ifelse(label != "", ifelse(
    rno > 0, ggcolor["Elevated"], ggcolor["Reduced"]), "#BBBBBB")) %>% 
  mutate(shape = 16) %>% # filled circle no outline
  mutate(size = ifelse(label != "", 9, 5)) %>% 
  mutate(alpha = ifelse(label != "", .7, .3)) %>% 
  mutate(hjust = ifelse(hsa > 0, 1, -0.05)) %>% 
  mutate(vjust = 1) %>% arrange(size) %>%
  mutate(hsa = pmin(pmax(hsa,-2.5),2.5)) %>% 
  mutate(rno = pmin(pmax(rno,-2.5),2.5)) %>% 
  ggplot(aes(x = hsa, y = rno, label = label, size = size, color = color)) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_vline(xintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_point(aes(shape = shape, alpha = alpha)) + 
  scale_shape_identity() + scale_alpha_identity() + 
  scale_size_identity() + scale_color_identity() + 
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2 )) + 
  theme_minimal(base_size = 24) + 
  theme(legend.position = "none") + xlab("") + ylab("")
ggsave(paste0("ncomm_blais_bonus_fig6g.pdf"), timbr.caffeine.plot, width = 6, height = 6)
ggsave(paste0("ncomm_blais_bonus_fig6h.pdf"), timbr.theophylline.plot, width = 6, height = 6)
ggsave(paste0("ncomm_blais_bonus_fig6g_labeled.pdf"), timbr.caffeine.plot + 
         geom_text(aes(label = label, hjust = hjust, vjust = vjust), size = 6, color = "black"), width = 6, height = 6)
ggsave(paste0("ncomm_blais_bonus_fig6h_labeled.pdf"), timbr.theophylline.plot + 
         geom_text(aes(label = label, hjust = hjust, vjust = vjust), size = 6, color = "black"), width = 6, height = 6)

# Correlation coefficients
timbr.drug.summary %>% 
  filter(drug_id %in% c("caffeine", "theophylline")) %>% 
  select(starts_with("drug_")) %>% distinct


caffeine.validation.data.entry = bind_rows(list( 
  data_frame(met = "ornithine", ctl = 224, trt = 219, sig = 0, pval = 1, ctl_sd = 28, trt_sd = 30, pergramliver = 1e-9, unit = "nmol/g liver", pubmed_src = "table1"),
  data_frame(met = "citrulline", ctl = 80, trt = 91, sig = 1, pval = 0.05, ctl_sd = 8, trt_sd = 12, pergramliver = 1e-9, unit = "nmol/g liver", pubmed_src = "table1"),
  data_frame(met = "aspartate", ctl = 614, trt = 704, sig = 1, pval = 0.002, ctl_sd = 51, trt_sd = 55, pergramliver = 1e-9, unit = "nmol/g liver", pubmed_src = "table1"),
  data_frame(met = "arginine", ctl = 70, trt = 72, sig = 0, pval = 1, ctl_sd = 7, trt_sd = 5, pergramliver = 1e-9, unit = "nmol/g liver", pubmed_src = "table1"),
  data_frame(met = "urea", ctl = 2.9, trt = 4.9, sig = 1, pval = 0.001, ctl_sd = .4, trt_sd = .3, pergramliver = 1e-6, unit = "umol/ min / g wet wt", pubmed_src = "table1"),
  data_frame(met = "glutamine", ctl = 5.0, trt = 5.1, sig = 0, pval = 1, ctl_sd = 0.7, trt_sd = 0.5, pergramliver = 1e-6, unit = "umol/g liver", pubmed_src = "table2"),
  data_frame(met = "glutamate", ctl = 8.7, trt = 6.4, sig = -1, pval = 0.001, ctl_sd = 0.4, trt_sd = 0.4, pergramliver = 1e-6, unit = "umol/g liver", pubmed_src = "table2"),
  data_frame(met = "pyruvate", ctl = .28, trt = .25, sig = 0, pval = 1, ctl_sd = 0.03, trt_sd = 0.03, pergramliver = 1e-6, unit = "umol/g liver", pubmed_src = "table3"),
  data_frame(met = "L-lactate", ctl = 3.3, trt = 3.5, sig = 0, pval = 1, ctl_sd = 0.4, trt_sd = 0.4, pergramliver = 1e-6, unit = "umol/g liver", pubmed_src = "table3"),
  data_frame(met = "acetoacetate", ctl = .07, trt = .06, sig = 0, pval = 1, ctl_sd = 0.01, trt_sd = 0.01, pergramliver = 1e-6, unit = "umol/g liver", pubmed_src = "table3"),
  data_frame(met = "(R)-3-hydroxybutanoate", ctl = .16, trt = .17, sig = 0, pval = 1, ctl_sd = 0.02, trt_sd = 0.03, pergramliver = 1e-6, unit = "umol/g liver", pubmed_src = "table3"),
  data_frame(met = "alanine", ctl = 3.5, trt = 2.0, sig = -1,  ctl_sd = 0.5, trt_sd = 0.5, pergramliver = 1e-6, unit = "umol/g liver", pubmed_src = "table4", excluded = T,
             pubmed_note = "excluded from correlation because article mentions alanine decreased in vivo but was the only metabolite to not change similarly in vitro"))) %>% 
  mutate(drug_id = "caffeine", organism_id = "rno", pubmed_id = "2764993") %>%
  mutate(included = !(!is.na(excluded) & excluded)) %>% filter(included) %>%
  mutate(experimental_logfc = log2(trt / ctl))

caffeine.validation.prediction = timbr.predictions %>% 
  select(organism_id, drug_id, rxn_id, rxn_name, met, production_score) %>% 
  arrange(production_score) %>% 
  inner_join(caffeine.validation.data.entry) %>% 
  mutate(prediction = production_score, validation = experimental_logfc)

caffeine.validation.prediction %>% group_by(drug_id, organism_id) %>% 
  summarize(cor_coef = cor.test(prediction, validation)$estimate,
            cor_pvalue = cor.test(prediction, validation)$p.value) %>% ungroup

caffeine.validation.plot = caffeine.validation.prediction %>% 
  mutate(label = ifelse(met == "(R)-3-hydroxybutanoate", "3-hydroxybutyrate", met)) %>% 
  mutate(label = ifelse(met == "L-lactate", "lactate", label)) %>% 
  mutate(color = ifelse(sig > 0, ggcolor["similar"], ifelse(sig < 0, ggcolor["opposite"], "#BBBBBB"))) %>%
  mutate(shape = 16) %>%
  mutate(size = ifelse(label != "", 9, 5)) %>% 
  mutate(alpha = ifelse(label != "", 0.7, 0.7)) %>% 
  mutate(vjust = 1, hjust = ifelse(validation > 0, 1, -0.05)) %>% 
  ggplot(aes(x = validation, y = prediction, label = label, size = size, color = color)) + 
  geom_hline(yintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_vline(xintercept = 0, linetype = "solid", alpha = axis.alpha, size = axis.size) + 
  geom_point(aes(shape = shape, alpha = alpha)) + 
  scale_shape_identity() + scale_alpha_identity() + 
  scale_size_identity() + scale_color_identity() + 
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(limits = c(-2.6, 2.6), breaks = c(-2, -1, 0, 1, 2)) +
  theme_minimal(base_size = 24) + 
  theme(legend.position = "none", aspect.ratio = 1.5) + xlab("") + ylab("")
caffeine.validation.plot
ggsave(paste0("ncomm_blais_final_fig5.pdf"), caffeine.validation.plot, 
       width = 6, height = 8)
ggsave(paste0("ncomm_blais_final_fig5_labeled.pdf"), caffeine.validation.plot + 
         geom_text(aes(vjust = vjust, hjust = hjust), size = 6, alpha = 1), 
       width = 6, height = 8)
