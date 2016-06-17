# This script reproduces parts of Figures 2 and 3
# of the manuscript. 

# All helper functions defined in this source file
# have the prefix: ef (edik function)

options(java.parameters = "-Xmx4g")
options(stringsAsFactors = FALSE)
source("ncomm_helper.R")

library(reshape2)
library(dplyr)
library(ggplot2)

path.manuscript.directory = ""

gpr.color = c(`non-enzymatic` = "#7F7F7F", shared = "#674EA7", 
              `rat-specific` = "#C0504D", `human-specific` = "#4F81BD",
              rno = "#C0504D", hsa = "#4F81BD", 
              similar = "#674EA7", opposite = "#F79646")

# paste0(path.rdata.directory, "ncomm_blais_supplementary_table_s1.txt", sep = "\t", quote = "", header = T, check.names = F)
library(xlsx)

# Orthology information mapped to reactions from the original HMR2 reconstruction in Supplementary Table 2
hmr2.rxn.gene = read.xlsx2(paste0(path.manuscript.directory,"ncomm_blais_supplementary_table2.xlsx"),
                           sheetIndex = 1, startRow = 2) %>% as.tbl

# Curated reaction information in Supplementary Table 3
rxn.info.load = read.xlsx2(paste0(path.manuscript.directory,"ncomm_blais_supplementary_table3.xlsx"),
                           sheetIndex = 1, startRow = 2) %>% as.tbl

rxn.info = rxn.info.load %>% 
  mutate(enabled = ifelse(!is.na(enabled), enabled == "true", F)) %>% 
  mutate(enabled_hsa = ifelse(!is.na(enabled_hsa), enabled_hsa == "true", F)) %>% 
  mutate(enabled_rno = ifelse(!is.na(enabled_rno), enabled_rno == "true", F))%>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))



gpr.size.curated = rxn.info %>% filter(enabled) %>% select(rxn_id,n_hsa, n_rno)
gpr.limit = 9

gpr.size.curated.count = gpr.size.curated %>% 
  mutate(rno = paste0(pmin(n_rno, gpr.limit), ifelse(n_rno >= gpr.limit, "+", "")),
         hsa = paste0(pmin(n_hsa, gpr.limit), ifelse(n_hsa >= gpr.limit, "+", ""))) %>%
  count(rno,hsa) %>% ungroup %>% rename(n_gpr = n) %>%
  mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
    rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))

fig2a.plot = gpr.size.curated.count %>% 
  ggplot(aes(x = factor(hsa),y = factor(rno))) + 
  geom_tile(aes(fill = organism, alpha = log10(n_gpr), width = 0.9, height = 0.9)) + #width = 0.9, height = 0.9) + #, height = 0.9
  scale_fill_manual(values = gpr.color) + scale_color_manual(values = gpr.color) +
  scale_alpha_continuous(range = c(.2,.9)) + 
  scale_x_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  scale_y_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  xlab("Human GPR size") + ylab("Rat GPR size") + 
  theme_bw(base_size = 12) +
  theme(title = element_text(size = 12), legend.position = "none",
        axis.text.y = element_text(hjust = .7),
        panel.grid.major = element_line(size = NA),
        panel.grid.minor = element_line(size = NA))
fig2a.plot
ggsave(filename = paste0("ncomm_blais_fig2a.pdf"), plot = fig2a.plot, width = 2.2, height = 2.2)



hmr2.rxns = rxn.info %>% filter(enabled) %>% 
  select(rxn_id,hmr2_id) %>% 
  filter(!is.na(hmr2_id)) %>% filter(nchar(hmr2_id) > 0)

orthology.data = hmr2.rxn.gene %>% 
  mutate(hsa = as.character(hsa_gene), rno = as.character(rno_gene),
         orthology_src = orthology_database_list,
         orthology_count = as.numeric(orthology_database_count),
         hsa_score = as.numeric(hsa_gene_score), 
         rno_score = as.numeric(rno_gene_score)) %>%
  mutate(orthology_count = ifelse(!is.na(orthology_count), orthology_count, -1)) %>%
  mutate(hsa_rank = ifelse(!is.na(hsa_rank), hsa_rank, Inf)) %>%
  mutate(hsa_score = ifelse(!is.na(hsa_score), hsa_score, -1)) %>%
  mutate(rno_rank = ifelse(!is.na(rno_rank), rno_rank, Inf)) %>%
  mutate(rno_score = ifelse(!is.na(rno_score), rno_score, -1))

# GPR conversion algorithm filters orthology annotations
# based on two parameters:
########################## database.minimum
# database.minimum = the mininum number of database 
# in which a unique ortholog pair must occur.
# orthologs were obtained from 5 databases,
# so we assigned a higher confidence to 
# human genes mapped to rat genes in multiple databases
# this value ranged from 0 - 5

########################## rank.maximum
# rank.maximum = maximum number of allowed 
# rat orthologs to be replaced by a human gene
# the score that was used to rank rat genes 
# was determined by the following equation:
# gene score = 5*(status+annotation+evidence+type)
# this value ranged from 0 - 100

# The prioritization of orthologs depended on bulk downloads 
# from multiple databases. The final scores are included, 
# but we did not receive authorization to distribute 
# these databases in whole so the methods that we used are
# described below:
###### type = NCBI gene type
# protein-coding evidence = 5
# no protein-coding evidence = 0
###### status = Ensembl gene summary status
# KNOWN = 5
# KNOWN_BY_PROJECTION = 4
# NOVEL = 3
# PUTATIVE = 1
# no data = 0
###### annotation = UniProt annotation score 
# value assigned was same as UniProt annotation score, otherwise 0
###### evidence = UniProt gene evidence
# Evidence at protein level = 5
# Evidence at transcript level = 4
# Inferred from homology = 3
# Predicted = 2
# Uncertain = 1
# no data = 0

# Apply a reasonable range of cutoffs for our two parameters:
gpr.conversion.result = lapply(c(1:9),function(rank.maximum) {
  lapply(c(1:5),function(database.minimum) {
    hmr2.rxns %>% select(rxn_id) %>% distinct %>% 
      mutate(cutoff_rank = rank.maximum, 
             cutoff_database = database.minimum) %>% 
      left_join(orthology.data %>% 
                  mutate(rno = ifelse(hsa_rank <= rank.maximum, rno, "0"),
                         rno = ifelse(orthology_count >= database.minimum, rno, "0")) %>% 
                  select(rxn_id, hsa, rno, hsa_score,rno_score), by = "rxn_id")}) %>% 
    rbind_all}) %>% rbind_all %>%
  group_by(rxn_id) %>% mutate(n_hsa = length(unique(setdiff(hsa,c(NA, "0", ""))))) %>% ungroup

gpr.orthology.unfiltered = hmr2.rxns %>% select(rxn_id) %>% distinct %>% 
  mutate(cutoff_rank = Inf, cutoff_database = 0) %>% 
  left_join(orthology.data %>% select(rxn_id, hsa, rno, hsa_score,rno_score), by = "rxn_id") %>% 
  group_by(rxn_id) %>% 
  mutate(n_hsa = length(unique(setdiff(hsa,c(NA, "0", "")))),
         n_rno = length(unique(setdiff(rno,c(NA, "0", ""))))) %>% ungroup

gpr.size.limit = 5
gpr.conversion.count = gpr.conversion.result %>% 
  select(rxn_id,n_hsa,cutoff_database,cutoff_rank) %>% distinct %>% 
  left_join(gpr.conversion.result %>% 
              mutate(rno = ifelse(!is.na(rno) & nchar(rno)>0,rno,"0")) %>%
              filter(rno != "0") %>% 
              group_by(rxn_id,n_hsa,cutoff_database,cutoff_rank) %>% 
              summarize(n_rno = n_distinct(rno)) %>% ungroup) %>%
  mutate(n_rno = ifelse(is.na(n_rno),0,n_rno),n_hsa = ifelse(is.na(n_hsa),0,n_hsa),
         n_diff = n_hsa - n_rno,
         facet_database = paste0(cutoff_database, "+ databases"), 
         facet_rank = paste0("Top ", cutoff_rank, " orthologs"),
         rno = paste0(pmin(n_rno, gpr.size.limit), ifelse(n_rno >= gpr.size.limit, "+", "")),
         hsa = paste0(pmin(n_hsa, gpr.size.limit), ifelse(n_hsa >= gpr.size.limit, "+", "")),
         rxn_organism = ifelse(n_hsa > 0 & n_rno > 0, "shared",ifelse(
           n_rno > 0, "rat-specific",ifelse(n_hsa > 0, "human-specific", "non-enzymatic"))),
         n_bias = ifelse(n_rno > 0 & n_hsa > 0, sign(n_diff), 0),
         dir_bias = ifelse(n_rno > 0 & n_hsa > 0,ifelse(n_diff == 0,"shared",ifelse(n_diff > 0,"hsa","rno")),"none"))

# Heatmap version of Supplementary Figure 3
gpr.conversion.count %>% filter(n_rno > 0, n_hsa > 0) %>% 
  count(cutoff_database, cutoff_rank, dir_bias) %>% ungroup %>% 
  dcast(cutoff_database + cutoff_rank ~ dir_bias,fill = 0,value.var = "n") %>% as.tbl %>%
  ggplot(aes(x = factor(cutoff_rank), y = factor(cutoff_database))) + 
  geom_tile(aes(fill = hsa - rno)) + 
  scale_fill_gradient2(limits=c(-3000, 3000)) + 
  geom_text(aes(label = rno - hsa))

gpr.conversion.summary = gpr.conversion.count %>% 
  filter(n_rno > 0,n_hsa > 0) %>% 
  select(facet_database,facet_rank,rxn_id,n_hsa,n_rno) %>% distinct %>%
  group_by(facet_database,facet_rank) %>% 
  summarize(mean_gpr_ratio = mean(n_rno / n_hsa),
            mean_gpr_difference = mean(n_rno - n_hsa),
            n_gpr_rno_bigger = sum((n_rno > n_hsa)),
            n_gpr_hsa_bigger = sum((n_rno < n_hsa))) %>%
  ungroup

# The cutoffs were selected based on this summary statistic:
gpr.conversion.summary %>% arrange(abs(n_gpr_rno_bigger - n_gpr_hsa_bigger))

# An ideal balance was found when:
# orthology was annotated in 2+ databases and,
# up to rat 2 orthologs was selected.
gpr.conversion.threshold = gpr.conversion.count %>% 
  mutate(shared = (n_rno > 0) & (n_hsa > 0)) %>% 
  group_by(facet_database,facet_rank) %>% 
  mutate(gpr_shared = sum(shared),
         gpr_ratio = mean(n_rno[shared] / n_hsa[shared]),
         gpr_diff = mean(n_rno[shared] - n_hsa[shared]),
         n_gpr_rno_bigger = sum((n_rno[shared] > n_hsa[shared])),
         n_gpr_hsa_bigger = sum((n_rno[shared] < n_hsa[shared])),
         # sum_hsa = sum(n_hsa[n_rno > 0 & n_hsa > 0]),
         # sum_rno = sum(n_rno[n_rno > 0 & n_hsa > 0]),
         gpr_human_specific = sum(n_rno == 0 & n_hsa > 0),
         # threshold_label = paste0("",round(gpr_ratio,2)),
         threshold_label = paste0("",round(n_gpr_rno_bigger - n_gpr_hsa_bigger,2)),
         threshold_label_x = ifelse(gpr_ratio < 1,"1","4"),
         threshold_label_y = ifelse(gpr_ratio < 1,"5+","0")) %>% ungroup %>% 
  arrange(rxn_id,cutoff_database, cutoff_rank) %>% 
  mutate(cutoff_selected = cutoff_database == 2 & cutoff_rank == 2)

# Create plot labels for different thresholds
gpr.threshold.label = gpr.conversion.threshold %>% 
  filter(cutoff_rank <= 5) %>% 
  mutate(sum_diff = n_gpr_hsa_bigger - n_gpr_rno_bigger,
         txt_value = round(abs(sum_diff), 2), txt_dir = sign(sum_diff), 
         txt_color = ifelse(txt_dir < 0, "rat-specific", "human-specific")) %>%
  select(facet_database, facet_rank, hsa = threshold_label_x, rno = threshold_label_y,
         label = threshold_label, txt_dir, txt_value, txt_color) %>% distinct

plot.gpr.sensitivity = gpr.conversion.threshold %>% 
  filter(cutoff_rank < 6) %>% 
  ggplot(aes(x = factor(hsa),y = factor(rno))) + 
  geom_point(aes(color = rxn_organism), alpha = 0.3, position = "jitter") +
  geom_abline(linetype = "dashed", alpha= 0.5) +
  geom_text(data = gpr.threshold.label, aes(label = label, color = txt_color),size = 10) +
  scale_color_manual(values = gpr.color) +
  theme_minimal(base_size = 32) + 
  theme(legend.position = "none", strip.background = element_rect(fill = 'white')) + 
  xlab("Human GPR size") + ylab("Rat GPR size") +
  facet_grid(facet_database~facet_rank, as.table = F)
plot.gpr.sensitivity
ggsave(filename = "ncomm_blais_supplementary_fig3.pdf",plot = plot.gpr.sensitivity, width = 15,height = 15)

##### Bonus analyses
## Without filtering orthology data, GPR sizes are disproportionately large in rats
# Calculate gpr sizes obtained by replacing all human genes with rat orthologs

gpr.size.filtered = hmr2.rxns %>% 
  left_join(gpr.conversion.threshold %>% filter(cutoff_database == 2, cutoff_rank == 2))

gpr.size.unfiltered = hmr2.rxns %>% 
  left_join(orthology.data) %>% 
  group_by(rxn_id) %>% 
  summarize(n_hsa = length(unique(hsa[!hsa %in% c("NA","0","",NA)])),
            n_rno = length(unique(rno[!rno %in% c("NA","0","",NA)]))) %>% 
  ungroup

gpr.size.unfiltered.count = gpr.size.unfiltered %>% 
  mutate(rno = paste0(pmin(n_rno, gpr.limit), ifelse(n_rno >= gpr.limit, "+", "")),
         hsa = paste0(pmin(n_hsa, gpr.limit), ifelse(n_hsa >= gpr.limit, "+", ""))) %>%
  count(rno,hsa) %>% ungroup %>% rename(n_gpr = n) %>%
  mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
    rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))

gpr.size.unfiltered.plot = gpr.size.unfiltered.count %>% 
  ggplot(aes(x = factor(hsa),y = factor(rno))) + 
  geom_tile(aes(fill = organism, alpha = log10(n_gpr), width = 0.9, height = 0.9)) + #width = 0.9, height = 0.9) + #, height = 0.9
  scale_fill_manual(values = gpr.color) + scale_color_manual(values = gpr.color) +
  scale_alpha_continuous(range = c(.2,.9)) + 
  scale_x_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  scale_y_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  xlab("Human GPR size") + ylab("Rat GPR size") + 
  theme_bw(base_size = 12) +
  theme(title = element_text(size = 12), legend.position = "none",
        axis.text.y = element_text(hjust = .7),
        panel.grid.major = element_line(size = NA),
        panel.grid.minor = element_line(size = NA))
gpr.size.unfiltered.plot

ggsave(filename = "ncomm_blais_bonus_fig2a_unfiltered.pdf",plot = gpr.size.unfiltered.plot, width = 2.2,height = 2.2)

gpr.size.filtered.count = gpr.size.filtered %>% 
  mutate(rno = paste0(pmin(n_rno, gpr.limit), ifelse(n_rno >= gpr.limit, "+", "")),
         hsa = paste0(pmin(n_hsa, gpr.limit), ifelse(n_hsa >= gpr.limit, "+", ""))) %>%
  count(rno,hsa) %>% ungroup %>% rename(n_gpr = n) %>%
  mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
    rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))

gpr.size.filtered.plot = gpr.size.filtered.count %>% 
  ggplot(aes(x = factor(hsa),y = factor(rno))) + 
  geom_tile(aes(fill = organism, alpha = log10(n_gpr), width = 0.9, height = 0.9)) + #width = 0.9, height = 0.9) + #, height = 0.9
  scale_fill_manual(values = gpr.color) + scale_color_manual(values = gpr.color) +
  scale_alpha_continuous(range = c(.2,.9)) + 
  scale_x_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  scale_y_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  xlab("Human GPR size") + ylab("Rat GPR size") + 
  theme_bw(base_size = 12) +
  theme(title = element_text(size = 12), legend.position = "none",
        axis.text.y = element_text(hjust = .7),
        panel.grid.major = element_line(size = NA),
        panel.grid.minor = element_line(size = NA))
gpr.size.filtered.plot

ggsave(filename = "ncomm_blais_bonus_fig2a_filtered.pdf",plot = gpr.size.filtered.plot, width = 2.2,height = 2.2)


## An alternative choice would have been to just pick one database.
# If we had to pick one, it would have been Homologene.
database.ids = c("kegg", "uniprot", "homologene", "ensembl", "rgd")
orthology.databases = database.ids %>% setNames(database.ids) %>% 
  lapply(function(x.db) orthology.data %>% filter(grepl(x.db,orthology_src))) %>%
  bind_rows(.id = "orthology_db")



# GPR sizes and rules from HMR2 after converting Ensembl genes into Entrez genes.
gpr.size.database.hsa = hmr2.rxns %>% select(rxn_id) %>% distinct %>% mutate(organism_id = "hsa") %>% 
  left_join(orthology.databases %>% select(rxn_id, hsa) %>% distinct %>% 
              filter(!is.na(hsa)) %>% filter(hsa != "") %>% filter(hsa != "0") %>%
              select(rxn_id, hsa) %>% distinct %>% 
              arrange(rxn_id, hsa) %>% 
              group_by(rxn_id) %>% 
              summarize(n_gene = length(unique(hsa)), gpr_rule = paste0(unique(hsa), collapse=  ";")) %>% ungroup %>% 
              mutate(organism_id = "hsa")) %>%
  mutate(n_gene = ifelse(!is.na(n_gene), n_gene, 0), gpr_rule = ifelse(!is.na(gpr_rule), gpr_rule, ""))

gpr.size.database.rno = database.ids %>% setNames(database.ids) %>% 
  lapply(function(x.db) { 
    hmr2.rxns %>% select(rxn_id) %>% distinct %>% mutate(organism_id = "rno", orthology_db = x.db) %>% 
      left_join(orthology.databases %>% select(orthology_db, rxn_id, rno) %>% distinct %>% 
                  filter(!is.na(rno)) %>% filter(rno != "") %>% filter(rno != "0") %>%
                  select(orthology_db, rxn_id, rno) %>% distinct %>% 
                  arrange(orthology_db, rxn_id, rno) %>% 
                  group_by(orthology_db, rxn_id) %>% 
                  summarize(n_gene = length(unique(rno)), gpr_rule = paste0(unique(rno), collapse=  ";")) %>% ungroup %>% 
                  mutate(organism_id = "rno")) %>%
      mutate(n_gene = ifelse(!is.na(n_gene), n_gene, 0), gpr_rule = ifelse(!is.na(gpr_rule), gpr_rule, ""))
  }) %>% bind_rows()

gpr.size.database.hsa = database.ids %>% setNames(database.ids) %>% 
  lapply(function(x.db) gpr.size.unfiltered %>% select(rxn_id, n_gene = n_hsa) %>% mutate(organism_id = "hsa", orthology_db = x.db)) %>% bind_rows()

gpr.size.database = bind_rows(list(
  bind_rows(list(gpr.size.database.hsa, gpr.size.database.rno)) %>% 
    mutate(variable = paste0("n_",organism_id)) %>% 
    dcast(rxn_id + orthology_db ~ variable, value.var = "n_gene", fill = 0, drop = F) %>% ef_df,
  gpr.size.unfiltered %>% mutate(orthology_db = "all databases")))


# Plot GPR size comparisons for individual databases
# gpr.size.limit = 5
# gpr.size.database = gpr.size.database %>%
#   mutate(shared = (n_hsa > 0) & (n_rno > 0), fail = (n_hsa > 0) & (n_rno <= 0)) %>% 
#   group_by(orthology_db) %>% 
#   mutate(rno = paste0(pmin(n_rno, gpr.size.limit), ifelse(n_rno >= gpr.size.limit, "+", "")),
#          hsa = paste0(pmin(n_hsa, gpr.size.limit), ifelse(n_hsa >= gpr.size.limit, "+", "")),
#          rxn_organism = ifelse((n_hsa > 0) & (n_rno > 0), "shared",ifelse(
#            n_rno > 0, "rat-specific",ifelse(n_hsa > 0, "human-specific", "non-enzymatic"))))
# 
# plot.gpr.db.conversion = gpr.db.orthology  %>%
#   ggplot(aes(x = factor(hsa),y = factor(rno))) + 
#   geom_point(aes(color = rxn_organism), alpha = 0.3, size = 2, position = "jitter") +
#   geom_abline(linetype = "dashed",alpha= 0.5) +
#   scale_color_manual(values = gpr.color) +
#   theme_minimal(base_size = 32) + 
#   theme(legend.position = "none") + #, strip.background = element_rect(fill = 'white')
#   xlab("Human GPR size") + ylab("Rat GPR size") +
#   facet_wrap(~orthology_db)
# plot.gpr.db.conversion
# # It looks like Homologene would have been the best single database choice;
# # however, any individual database would produce a lot of human-specific reactions


gpr.size.database.count = gpr.size.database %>% 
  mutate(rno = paste0(pmin(n_rno, gpr.limit), ifelse(n_rno >= gpr.limit, "+", "")),
         hsa = paste0(pmin(n_hsa, gpr.limit), ifelse(n_hsa >= gpr.limit, "+", ""))) %>%
  count(orthology_db, rno,hsa) %>% ungroup %>% rename(n_gpr = n) %>%
  mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
    rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))


gpr.size.database.plot = gpr.size.database.count %>% 
  ggplot(aes(x = factor(hsa),y = factor(rno))) + 
  geom_tile(aes(fill = organism, alpha = log10(n_gpr), width = 0.9, height = 0.9)) + 
  scale_fill_manual(values = gpr.color) + scale_color_manual(values = gpr.color) +
  scale_alpha_continuous(range = c(.2,.9)) + 
  scale_x_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  scale_y_discrete(breaks = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit), 
                   labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
  xlab("Human GPR size") + ylab("Rat GPR size") + 
  theme_minimal(base_size = 12) +
  theme(title = element_text(size = 12), legend.position = "none",
        axis.text.y = element_text(hjust = .7),
        panel.grid.major = element_line(size = NA),
        panel.grid.minor = element_line(size = NA)) + facet_wrap(~orthology_db)
gpr.size.database.plot
ggsave(filename = "ncomm_blais_bonus_fig2a_database.pdf",plot = gpr.size.database.plot, width = 6.6, height = 4.4)


gpr.size.database.scatterplot = gpr.size.database  %>%
  mutate(rno = paste0(pmin(n_rno, gpr.limit), ifelse(n_rno >= gpr.limit, "+", "")),
         hsa = paste0(pmin(n_hsa, gpr.limit), ifelse(n_hsa >= gpr.limit, "+", ""))) %>%
  mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
    rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic")))) %>%
  ggplot(aes(x = factor(hsa),y = factor(rno))) + 
  geom_point(aes(color = organism), alpha = 0.3, size = 2, position = "jitter") +
  geom_abline(linetype = "dashed",alpha= 0.5) +
  scale_color_manual(values = gpr.color) +
  theme_minimal(base_size = 32) + 
  theme(legend.position = "none") + 
  xlab("Human GPR size") + ylab("Rat GPR size") +
  facet_wrap(~orthology_db)
ggsave(filename = "ncomm_blais_bonus_fig2a_database_scatterplot.pdf",plot = gpr.size.database.scatterplot, width = 6.6, height = 4.4)

# The analyses below were used to count the numbers of orthologs/genes for reference in the text:

orthology.filtered = orthology.data %>% 
  semi_join(rxn.info %>% filter(enabled) %>% select(rxn_id))  %>% 
  anti_join(rxn.info %>% filter(n_gene_hsa == 0, n_gene_rno == 0) %>% select(rxn_id))  %>% 
  filter(orthology_count >= 2, hsa_rank <= 2)

orthology.unfiltered = orthology.data %>% 
  semi_join(rxn.info %>% filter(enabled) %>% select(rxn_id))  %>% 
  anti_join(rxn.info %>% filter(n_gene_hsa == 0, n_gene_rno == 0) %>% select(rxn_id))

# count the numbers of ortholog pairs, rat genes, and human genes before filtering
data_frame(hsa_genes = orthology.unfiltered %>% filter(hsa != "0") %>% count(hsa) %>% nrow,
           rno_orthologs = orthology.unfiltered %>% filter(rno != "0") %>% count(rno) %>% nrow,
           orthology_pairs = orthology.unfiltered %>% filter(rno != "0", hsa != "0") %>% count(hsa, rno) %>% ungroup %>% nrow)
# count the numbers of ortholog pairs, rat genes, and human genes after filtering
data_frame(hsa_genes = orthology.filtered %>% filter(hsa != "0") %>% count(hsa) %>% nrow,
           rno_orthologs = orthology.filtered %>% filter(rno != "0") %>% count(rno) %>% nrow,
           orthology_pairs = orthology.filtered %>% filter(rno != "0", hsa != "0") %>% count(hsa, rno) %>% ungroup %>% nrow)

# human-specific reactions after filtering
gpr.size.filtered %>% select(rxn_id, hsa, rno) %>% 
  inner_join(rxn.info %>% filter(enabled) %>% 
               select(rxn_id, rxn_formula_name, subsystem_id, subsystem_name, 
                      n_gene_hsa, n_gene_rno, enabled_hsa, enabled_rno)) %>% 
  filter(hsa != "0", rno == "0") %>% data.frame

# human-specific reactions after filtering
gpr.size.filtered %>% select(rxn_id, hsa, rno) %>% 
  inner_join(rxn.info %>% filter(enabled) %>% 
               select(rxn_id, rxn_formula_name, subsystem_id, subsystem_name, 
                      n_gene_hsa, n_gene_rno, enabled_hsa, enabled_rno)) %>% 
  count("human-nonenzymatic" = hsa == "0", "rat-nonenzymatic" = rno == "0")

# human-specific reactions after manual curation
gpr.size.filtered %>% select(rxn_id, hsa, rno) %>% 
  inner_join(rxn.info %>% filter(enabled, subsystem_id != "Biomass") %>% select(
    rxn_id, rxn_formula_name, subsystem_id, subsystem_name, n_gene_hsa, n_gene_rno, enabled_hsa, enabled_rno)) %>% 
  filter(enabled_hsa, enabled_rno, hsa != "0", rno == "0") %>% 
  arrange(subsystem_id, subsystem_name)

gpr.size.filtered %>% select(rxn_id, hsa, rno) %>% 
  inner_join(rxn.info %>% filter(enabled, subsystem_id != "Biomass") %>% 
               select(rxn_id, enabled_hsa, enabled_rno)) %>% 
  count(filtered_status = paste0(ifelse(hsa != "0","hsa",""), ifelse(rno != "0","rno","")), 
        curated_status = paste0(ifelse(enabled_hsa,"hsa",""), ifelse(enabled_rno,"rno","")))
# 5 that were non-enzymatic are now rat-specific; 7 that were human-specific are still human-specific
# 12 that were human specific are now shared

# subsystems with human-specific reactions after filtering orthology
gpr.size.filtered %>% 
  inner_join(rxn.info %>% filter(enabled) %>% select(rxn_id, subsystem_id, subsystem_name)) %>% 
  filter(hsa != "0", rno == "0") %>% 
  select(rxn_id,subsystem_id, subsystem_name) %>% distinct %>% count(subsystem_id, subsystem_name)

# total number of unique ortholog pairs before filtering
orthology.unfiltered %>% filter(rno != "0", hsa != "0") %>% count(hsa, rno) %>% nrow
# total number of unique ortholog pairs after filtering
orthology.filtered %>% filter(rno != "0", hsa != "0") %>% count(hsa, rno) %>% nrow

# total number of human genes with rat orthologs before filtering
orthology.unfiltered %>% filter(rno != "0", hsa != "0") %>% count(hsa) %>% nrow
# total number of human genes with rat orthologs after filtering
orthology.filtered %>% filter(rno != "0", hsa != "0") %>% count(hsa) %>% nrow

# total number of rat orthologs mapped to human genes before filtering
orthology.unfiltered %>% filter(rno != "0", hsa != "0") %>% count(rno) %>% nrow
# total number of rat orthologs mapped to human genes after filtering
orthology.filtered %>% filter(rno != "0", hsa != "0") %>% count(rno) %>% nrow

# GPR size bias before filtering
gpr.size.unfiltered %>% filter(n_hsa > 0 & n_rno > 0) %>% 
  count(n_hsa > n_rno, n_rno > n_hsa)
# GPR size bias after filtering
gpr.size.filtered %>% filter(n_hsa > 0 & n_rno > 0) %>% 
  count(n_hsa > n_rno, n_rno > n_hsa)
# GPR size bias after filtering and manual curation
gpr.size.curated %>% filter(n_hsa > 0 & n_rno > 0) %>% 
  count(n_hsa > n_rno, n_rno > n_hsa)


