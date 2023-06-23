
# Full analysis completed within Summers et al. 
# Driver genes in endometrial hyperplasia and carcinoma exhibit stage-specific mutation rates, changing selection intensities, and antagonistic epistasis

# Code collaboratively built by Mary Summers, Sem Asmelash, J. Nick Fisk, Jeffrey D. Mandell, and Vincent L. Cannataro 





# load cancer effect size and necessary packages ----

library(cancereffectsizeR) # v2.7.0 https://townsend-lab-yale.github.io/cancereffectsizeR/
library(ces.refset.hg19)
library(data.table)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(scales)


# loading all VCF and TCGA data data ----

# load VCF1 and VCF2 data sets from Li et al.
# data provided by Li et al.

# Li L, Yue P, Song Q, et al. Genome-wide mutation analysis in precancerous lesions of endometrial carcinoma. J Pathol 2021; 253: 119-128

VCF1 <- read_csv("input_data/VCF1_r.csv")
VCF2 <- read_csv("input_data/VCF2_r.csv")

# loading clinical data and selecting necessary columns
clinical <- fread("input_data/clinical.tsv")
clinical <- select(clinical, 
                   case_id, 
                   figo_stage)

# loading TCGA UCEC data
if(!file.exists("input_data/TCGA_ucec_data.maf")){
  get_TCGA_project_MAF(
    project = "TCGA-UCEC",
    filename = "input_data/TCGA_ucec_data.maf",
    test_run = FALSE,
    exclude_TCGA_nonprimary = TRUE)
}

TCGA_ucec_data <- read_tsv("input_data/TCGA_ucec_data.maf", comment = "#")

# processing VCF data ----

# removing extra rows
VCF1 <- VCF1[-c(1:6),]

# rename sample ID column of VCF2
VCF2 <- dplyr::rename(VCF2, "#Sample_ID" = "#SAHmple_ID")

# combining VCF data sets
VCF_all <- rbind(VCF1,VCF2)

# removing duplicate samples
VCF_all <- VCF_all |> 
  filter(!`#Sample_ID` %in% c("AH13(1)", "AH1(1)")) 

#removing samples where column Gene_Name has value "Gene_Name"
VCF_all <- subset(VCF_all, Gene_Name != "Gene_Name")


# renaming columns of VCF_all to be compatible with cancer effect size r
VCF_all <- dplyr::rename(VCF_all, "Start_Position" = "Start_position", 
                         "Reference_Allele" = "Ref", 
                         "Tumor_Sample_Barcode" = "#Sample_ID", 
                         "Tumor_Allele" = "Alt")


# assinging stage column to VCF_all, assigned as "AH" or "EC"
VCF_all <- VCF_all %>% 
  mutate(Stage = stringr::str_sub(Tumor_Sample_Barcode,1,2))

# pre load maf VCF and TCGA ----

# pre load VCF data
VCF_all_maf_data <- preload_maf(maf = VCF_all, 
                                refset = "ces.refset.hg19", 
                                keep_extra_columns = c("Stage", "Gene_Name"))
# removing samples where column Problem is equal to NA
VCF_all_maf_data <- VCF_all_maf_data[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
VCF_all_maf_data <- VCF_all_maf_data[germline_variant_site == F]

# keeping only samples that do not occur in repetitive regions 
VCF_all_maf_data <- VCF_all_maf_data[(repetitive_region == F | cosmic_site_tier %in% 1:3)]

# removing EC samples from VCF data because the figo stage is not specified
VCF_all_maf_data <- subset(VCF_all_maf_data, Stage == "AH")

# pre laod TCGA data
maf_ucec <- preload_maf(maf = TCGA_ucec_data, 
                        refset = "ces.refset.hg19", 
                        chain_file = "input_data/hg38ToHg19.over.chain", 
                        keep_extra_columns = c("case_id"))

# removing samples where column Problem is equal to NA
maf_ucec <- maf_ucec[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
maf_ucec <- maf_ucec[germline_variant_site == F]

# keeping only samples that do not occur in repetitive regions 
maf_ucec <- maf_ucec[(repetitive_region == F | cosmic_site_tier %in% 1:3)]


# combinging ucec clinical and maf data
maf_ucec <- maf_ucec %>% 
  left_join(clinical, by = "case_id")

# removing the word "Stage" from all values in column figo_stage
maf_ucec <- maf_ucec |> 
  mutate(figo_stage = gsub(pattern = 'Stage ',replacement= "", x = figo_stage))

# specifying Stage 1 and Rest of Stages
maf_ucec <- maf_ucec |> 
  filter(!is.na(figo_stage)) |> 
  mutate(Stage = case_when(
    figo_stage %in% c("I", "IA", "IB", "IC") ~ "Stage1",
    TRUE ~ "RestOfStages"))

# loading and processing CPTAC data ----

# loading sample and clinical data from download, instructions to be written in read.me file
# loading all maf data from CPTAC 3 project


if(!file.exists("input_data/cptac.maf")){
  get_TCGA_project_MAF(
    project = "CPTAC-3",
    filename = "input_data/cptac.maf",
    test_run = FALSE,
    exclude_TCGA_nonprimary = TRUE)
}


# loading maf data
cptac.maf <- read_tsv("cptac.maf", comment = "#")
manifest <- read_tsv("input_data/gdc_manifest.2023-02-28.txt")
cptac_clinical <- read_tsv("input_data/clinical_cart/clinical.tsv")


# selecting data from CPTAC maf file with tumor sample barcodes matching the 102 samples from manifest
maf_endo_cptac <- cptac.maf |> filter(source_file_id %in% manifest$id)


# joining clinical and sample data sets by case id
maf_endo_cptac <- left_join(x = maf_endo_cptac, y = cptac_clinical, by = "case_id")


# pre load maf of clinical and sample data 
maf_cptac <- preload_maf(maf = maf_endo_cptac, 
                         refset = "ces.refset.hg19", 
                         chain_file = "input_data/hg38ToHg19.over.chain", 
                         keep_extra_columns = c("case_id", "ajcc_pathologic_stage"))


# removing samples where column Problem is equal to NA
maf_cptac <- maf_cptac[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
maf_cptac <- maf_cptac[germline_variant_site == F]

# keeping only samples that do not occur in repetitive regions 
maf_cptac <- maf_cptac[(repetitive_region == F | cosmic_site_tier %in% 1:3)]

# adding stage column and specifying Stage 1 and Rest of Stages
maf_cptac <- maf_cptac |> 
  filter(!is.na(ajcc_pathologic_stage)) |> 
  mutate(Stage = case_when(
    ajcc_pathologic_stage %in% c("Stage I") ~ "Stage1",
    TRUE ~ "RestOfStages"))





# stage specific analysis of selection intensity ----

# initiate analysis
cesa <- CESAnalysis(refset = "ces.refset.hg19")

# load in data from CPTAC, TCGA and Li projects
cesa <- load_maf(cesa = cesa, maf = maf_cptac,sample_data_cols = "Stage")
cesa <- load_maf(cesa = cesa, maf = VCF_all_maf_data,sample_data_cols = "Stage")
cesa <- load_maf(cesa = cesa, maf = maf_ucec,sample_data_cols = "Stage")


# estimating mutation rates
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "UCEC", samples = cesa$samples[Stage == "AH"],save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "UCEC", samples = cesa$samples[Stage == "Stage1"], save_all_dndscv_output = T)
# cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "UCEC", samples = cesa$samples[Stage == "RestOfStages"], save_all_dndscv_output = T)


# selecting genes of interest
selected_genes <- c("KRAS", "PTEN", "ARID1A", "CTCF", "CTNNB1", "PIK3CA", "CHD4", "FGFR2", "PIK3R1")


library(ces.refset.hg19)
RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- cesa_samples_by_groups$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# selecting mutation rate data for samples in hyperplasia 
samples_in_hyperplasia <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

# selecting mutation rate data for samples in stage 1
samples_in_cancer1 <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_2$annotmuts$sampleID ))

# creating a data frame with mutation rate data for hyperplasia and stage 1
mut_rate_df <- tibble(gene = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_hyp_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_cancer1_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(hyp_mu = (exp_hyp_mu / n_syn_sites) / samples_in_hyperplasia) %>%
  mutate(cancer1_mu = (exp_cancer1_mu / n_syn_sites) / samples_in_cancer1) %>%
  mutate(cancer_greater = cancer1_mu > hyp_mu) -> 
  mut_rate_df

# defining rate 1 and rate 2 as mutation rates for hyperplasia and stage 1
rate_1 <- mut_rate_df|>
  select(gene, hyp_mu)
rate_2 <- mut_rate_df|>
  select(gene, cancer1_mu)

# change in mutation rate across stages
mut_rate_df <- mut_rate_df %>% 
  select(gene, hyp_mu, cancer1_mu) %>% 
  mutate(p_1 = hyp_mu / cancer1_mu) %>% 
  mutate(p_2 = 1 - p_1)

# saving "last" gene mutation rates into separate data frame, "last" rates meaning from last stage cancer1_mu
set_cancer_rates <- mut_rate_df %>%
  select(gene, cancer1_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates to highest rates from cancer1_mu
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = set_cancer_rates, missing_genes_take_nearest = T) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "UCEC")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)


# selecting PIK3CA and PIK3R1 variants with maf prevalence > 1, including stop mutations for PIK3R1 (TSG https://www.oncokb.org/cancerGenes)
pik3ca_variants <- cesa_samples_by_groups$variants[gene == "PIK3CA" & maf_prevalence > 1 & intergenic == F]
pik3r1_variants <- cesa_samples_by_groups$variants[gene == "PIK3R1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
# selecting KRAS and FGFR2 variants with maf prevalence > 1
kras_variants <- cesa_samples_by_groups$variants[gene == "KRAS" & maf_prevalence > 1 & intergenic == F]
fgfr2_variants <- cesa_samples_by_groups$variants[gene == "FGFR2" & maf_prevalence > 1 & intergenic == F]
pten_variants <- cesa_samples_by_groups$variants[gene == "PTEN" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
arid1a_variants <- cesa_samples_by_groups$variants[gene == "ARID1A" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
ctcf_variants <- cesa_samples_by_groups$variants[gene == "CTCF" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
ctnnb1_variants <- cesa_samples_by_groups$variants[gene == "CTNNB1" & maf_prevalence > 1 & intergenic == F]
chd4_variants <- cesa_samples_by_groups$variants[gene == "CHD4" & maf_prevalence > 1 & intergenic == F]

stage_specific_variants <- rbind(pik3ca_variants,
                                 pik3r1_variants,
                                 kras_variants,
                                 fgfr2_variants,
                                 pten_variants,
                                 arid1a_variants,
                                 ctcf_variants,
                                 ctnnb1_variants,
                                 chd4_variants)



compound <- define_compound_variants(cesa = cesa_samples_by_groups, 
                                     variant_table = stage_specific_variants,
                                     by = "gene", merge_distance = Inf)



# feed compound variants into loop which calculates stage specific selection for each gene using "p" values from above

source("R/new_sequential_lik.R")

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  cesa_samples_by_groups <- ces_variant(cesa = cesa_samples_by_groups, variants = this_comp, model = sequential_lik_dev, 
                                        ordering_col = 'Stage', ordering = c('AH', 'Stage1'), 
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}


# creating selection plots with stage specific model ----

scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}

# selecting necessary data

selection_data <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_data <- selection_data |> 
  select(variant_name, starts_with("si"), starts_with("ci")) |>
  pivot_longer(cols = -variant_name, names_to = "data_type") |>
  mutate(stage = stringr::word(string = data_type, sep = "_",start = -1)) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(si_or_ci = stringr::word(string = data_type, sep = "_",start = 1, end=3)) |>
  mutate( si_or_ci = case_when(is.na(si_or_ci) ~ "si", TRUE ~ si_or_ci)) |>
  mutate (value = case_when (is.na(value)~0, TRUE~value))

# pivoting data set to create columns for gene, stage, si, and CIs
selection_data <- selection_data|> 
  select(-data_type) |>
  pivot_wider(values_from = value, names_from = si_or_ci)

# defining stages to be plotted
selection_data$stage <- factor(selection_data$stage, levels = c("AH","Stage1"))

# plotting results
selection_data |>
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), labels = c("Normal through AH", "AH through Stage 1")) + 
  theme_bw() + xlab("") + ylab("Cancer effect size") +
  theme(legend.position = "bottom", axis.text.x = element_blank(), legend.text=element_text(size=10)) +
  scale_y_continuous(labels = scientific) +
  theme(text = element_text(size = 20)) +
  expand_limits (y = 0)

# ggsave(filename = "figures/figure2.png", width = 7, height = 10)


# selection_data

# ,plot.margin = margin(0, 0, 0, 0, "cm")

limits1 <- c(0,1.5e3)
text_size <- 13

ARID1A_plot <- selection_data %>% 
  filter(variant_name == "ARID1A") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) + 
  theme_bw() + xlab("") + ylab("Cancer effect size") + labs(title = "ARID1A") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits1,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank())

CHD4_plot <- selection_data %>% 
  filter(variant_name == "CHD4") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) +  
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "CHD4") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits1,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank())



CTCF_plot <- selection_data %>% 
  filter(variant_name == "CTCF") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) +  
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "CTCF") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits1,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank())



FGFR2_plot <- selection_data %>% 
  filter(variant_name == "FGFR2") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) +  
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "FGFR2") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits1,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank())


PIK3CA_plot <- selection_data %>% 
  filter(variant_name == "PIK3CA") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) +  
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "PIK3CA") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits1,expand = c(0.01, 0)) +
  theme(axis.title.y = element_blank(),text = element_text(size = text_size),axis.ticks.x = element_blank()) 

PIK3R1_plot <- selection_data %>% 
  filter(variant_name == "PIK3R1") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) + 
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "PIK3R1") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits1,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank())


limits2 <- c(0,5e3)

KRAS_plot <- selection_data %>% 
  filter(variant_name == "KRAS") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) + 
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "KRAS") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits2,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) 

CTNNB1_plot <- selection_data %>% 
  filter(variant_name == "CTNNB1") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Normal through atypical hyperplasia", "Atypical hyperplasia Stage-1 endometrial carcinoma")) + 
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "CTNNB1") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits2,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank())

PTEN_plot <- selection_data %>% 
  filter(variant_name == "PTEN") %>% 
  ggplot(aes(x = stage, y = si, color = stage)) + 
  geom_point(size = 2.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) + 
  # facet_wrap(~variant_name, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("dark green", "purple"), 
                     labels = c("Progression from normal to atypical hyperplasia", 
                                "Progression from atypical hyperplasia to Stage-1 endometrial carcinoma")) +  
  theme_bw() + xlab("") + ylab("Cancer effect size") + 
  labs(title = "PTEN", color = NULL) + 
  theme(legend.position = "right", axis.text.x = element_blank(),
        legend.text=element_text(size=text_size)) +
  scale_y_continuous(labels = scientific,limits = limits2,expand = c(0.01, 0)) +
  theme(text = element_text(size = text_size),axis.ticks.x = element_blank()) + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank())



all_plots <- PIK3CA_plot + 
  PIK3R1_plot + 
  FGFR2_plot +  
  ARID1A_plot + 
  CHD4_plot  +  
  CTCF_plot + 
  KRAS_plot + 
  CTNNB1_plot + 
  PTEN_plot + 
  guide_area()  + 
  plot_layout(nrow = 2, ncol = 6,guides = "collect") 

ylab <- KRAS_plot$labels$y
labp <- ggplot(data.frame(l = ylab, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90,size=text_size * (5/14)) + 
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
coord_cartesian(clip = "off") 

KRAS_plot$labels$y <- " "




all_plots <-   KRAS_plot + 
  CTNNB1_plot + 
  PTEN_plot + 
  plot_spacer() + 
  guide_area() + 
  plot_spacer() + 
  PIK3CA_plot + 
  PIK3R1_plot + 
  FGFR2_plot +  
  ARID1A_plot + 
  CTCF_plot + 
  CHD4_plot  +
  plot_layout(nrow = 2, ncol = 6,guides = "collect") 

plot_all <- labp + all_plots + plot_layout(widths = c(1, 75))


ggsave(filename = "figures/figure_2_stages_selection.png", width = 14, height = 5, plot = plot_all)

# + 
# plot_annotation(tag_levels = "A")


# prevalences in stages ----

source("R/compound_prevalence_calc.R")

prevalence_df <- variant_prevalence_function(compound_obj = compound, cesa = cesa_samples_by_groups)


prevalence_df <- prevalence_df %>% 
  mutate(stages = case_when(
    stages == "AH" ~ "Atypical hyperplasia", 
    stages == "Stage1" ~ "Stage-1 EC"))

prevalence_df <- prevalence_df |>
  filter(stages %in% c("Atypical hyperplasia","Stage-1 EC")) |> 
  mutate(prop_tumors_with = tumors_with/total_samples) |>
  mutate(prop_tumors_without = tumors_without/total_samples)

prevalence_df_w <- prevalence_df %>% 
  pivot_wider(names_from = stages, values_from = c(prop_tumors_with, prop_tumors_without)) %>%
  fill(`prop_tumors_with_Atypical hyperplasia`) %>% 
  fill(`prop_tumors_with_Stage-1 EC`,.direction = "up") %>% 
  select(compound, `prop_tumors_with_Atypical hyperplasia`, `prop_tumors_with_Stage-1 EC`)

prevalence_df_w <- prevalence_df_w[rep(c(T,F),nrow(prevalence_df_w)/2),] # deduplicate

prevalence_df_w <- prevalence_df_w %>% 
  mutate(gene_name = str_remove(compound, pattern = "\\.1") )

# prevalence_df %>% 
#   pivot_longer(cols = c(prop_tumors_with,prop_tumors_without))
#   
library(ggrepel)

font_size <- 12
geom_font_size <- (5/14) * font_size

ggplot(prevalence_df, aes(x=stages, y=prop_tumors_with)) + 
  geom_point() + 
  geom_segment(data = prevalence_df_w, aes(x=1,xend = 2,
                                           y=`prop_tumors_with_Atypical hyperplasia`,
                                           yend=`prop_tumors_with_Stage-1 EC`)) + 
  geom_text_repel(data = prevalence_df_w,aes(x=1, y=`prop_tumors_with_Atypical hyperplasia`,label=gene_name),direction = "y",nudge_x = -0.01,hjust=1,size = geom_font_size) + 
  scale_color_viridis_d() + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,0.45),expand = c(0, 0)) + 
  scale_x_discrete(labels = c("Atypical\nhyperplasia","Stage-1\nendometrial carcinoma")) + 
  theme(text= element_text(size = font_size)) + 
  labs(y="Proportions of tissue samples with\ndriver mutations", x=NULL) -> 
  prev_plot


ggsave(filename = "figures/figure_1_prevalence.png",plot = prev_plot,height = 5,width = 4)

# epistasis ----

# resetting gene and trinucleotide mutation rates for estimation of epistasis over entire data set
cesa_epistasis <- gene_mutation_rates(cesa = cesa, covariates = "UCEC", save_all_dndscv_output = T)
cesa_epistasis <- trinuc_mutation_rates(cesa = cesa_epistasis, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)


RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- cesa_epistasis$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])


# creating data frame of mutation rates
mut_rate_epistasis <- tibble(gene = cesa_epistasis$dNdScv_results$rate_grp_1$genemuts$gene_name,
                             exp_mu = cesa_epistasis$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv)



mut_rate_epistasis$n_syn_sites = nsyn_sites[mut_rate_epistasis$gene]

samples_total <- length(unique(cesa_epistasis$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

mut_rate_epistasis <- mut_rate_epistasis %>% 
  mutate(mut_rate = (exp_mu / n_syn_sites) / samples_total) 

# saving mutation rates
set_rates_epistasis <- mut_rate_epistasis %>%
  select(gene, mut_rate) %>%
  data.table::setDT()

# clearing gene rates 
cesa_epistasis<- clear_gene_rates(cesa = cesa_epistasis)

# resetting gene mutation rates
cesa_epistasis <- set_gene_rates(cesa = cesa_epistasis, rates = set_rates_epistasis, missing_genes_take_nearest = T) 

# epistasis calculation for PIK3CA and PIK3R1

# selecting PIK3CA and PIK3R1 variants with maf prevalence > 1, including stop mutations for PIK3R1
pik3ca_variants <- cesa_epistasis$variants[gene == "PIK3CA" & maf_prevalence > 1]
pik3r1_variants <- cesa_epistasis$variants[gene == "PIK3R1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
# pik3r1_variants <- [gene == "PIK3R1" & (maf_prevalence > 1 | aa_alt == "STOP")]

# combining PIK3CA and PIK3R1 variants into one data set
pik3_variants <- rbind(pik3ca_variants, pik3r1_variants)

# defining compound variants
comp_pik3 <- define_compound_variants(cesa = cesa_epistasis, variant_table = pik3_variants, by = "gene", merge_distance = Inf)


shared_tcga_pik3 <- grep(x = intersect(comp_pik3$samples_with$PIK3CA.1, comp_pik3$samples_with$PIK3R1.1), pattern = "TCGA",value = T)

# cesa_epistasis$maf %>% 
#   filter(Unique_Patient_Identifier == shared_tcga_pik3[2]) %>%
#   filter(genes %in% c("PIK3CA","PIK3R1"))
# 
# cesa_epistasis$variants %>%
#   filter(variant_id == "5:67591246_A>G")
# 
# cesa_epistasis$maf %>% 
#   filter(Unique_Patient_Identifier == shared_tcga_pik3[3]) %>%
#   filter(genes %in% c("PIK3CA","PIK3R1"))



# epistasis calculation for PIK3CA and PIK3R1 relationships
cesa_pik3 <- ces_epistasis(cesa = cesa_epistasis, variants = comp_pik3, run_name = "PIK3R1_vs_PIK3CA")
epistasis_pik3 <- cesa_pik3$epistasis$PIK3R1_vs_PIK3CA


# finding samples where substitutions in PIK3CA and PIK3R1 co-occur

cesa_df <- as.data.frame(cesa_epistasis$maf)

# tumors with substitutions in PIK3CA and PIK3R1
tumors_of_interest_pik3 <- comp_pik3$samples_with$PIK3CA.1[comp_pik3$samples_with$PIK3CA.1 %in% comp_pik3$samples_with$PIK3R1.1]


# selecting only all PIK3CA and PIK3R1 variants with single nucleotide variants in tumors of interest
pik3_data <- cesa_df |>
  filter(genes %in% c("PIK3CA", "PIK3R1")) |>
  filter(variant_type == "snv") |>
  filter(Unique_Patient_Identifier %in% tumors_of_interest_pik3)

# selecting unique patient identifiers who have co-occurring substitutions in PIK3CA and PIK3R1
pik3_data <- pik3_data |>
  count(Unique_Patient_Identifier, genes) |>
  count(Unique_Patient_Identifier) |>
  filter(n > 1)

# effect size analysis for all PIK3R1 and PIK3CA variants 
cesa_pik3 <- ces_variant(cesa = cesa_pik3, variants = pik3_variants, run_name = "PIK3_effect_size")
pik3_selection <- cesa_pik3$selection$PIK3_effect_size 



# find what variants of PIK3R1 and PIK3CA exist in the tumors of interest

# selecting tumors of interest and the two PIK3 variants from the maf data
pik3_maf <- cesa_pik3$maf |>
  filter(Unique_Patient_Identifier %in% tumors_of_interest_pik3) |>
  filter(genes %in% c("PIK3CA", "PIK3R1"))

# combining the selection intensity for each variant with the tumors of interest maf data
pik3_maf_selection <- left_join(pik3_maf, pik3_selection, by = c("top_consequence" = "variant_name")) |>
  relocate(Unique_Patient_Identifier, selection_intensity, genes, top_consequence)



# epistasis calculation for KRAS and FGFR2 


# selecting KRAS and FGFR2 variants with maf prevalence > 1
kras_variants <- cesa_epistasis$variants[gene == "KRAS" & maf_prevalence > 1]
fgfr2_variants <- cesa_epistasis$variants[gene == "FGFR2" & maf_prevalence > 1]

# combining KRAS and FGFR2 variant data set
kras_fgfr2_variants <- rbind(kras_variants, fgfr2_variants)

# defining compound variants
comp_kras_fgfr2 <- define_compound_variants(cesa = cesa_epistasis, variant_table = kras_fgfr2_variants, by = "gene", merge_distance = Inf)

# epistasis calculations for relationship between KRAS and FGFR2
cesa_kras_fgfr2 <- ces_epistasis(cesa = cesa_epistasis, variants = comp_kras_fgfr2, run_name = "KRAS_vs_FGFR2")
epistasis_kras_fgfr2 <- cesa_kras_fgfr2$epistasis$KRAS_vs_FGFR2


# finding samples that have substitutions in both genes

# tumors with substitutions in KRAS and FGFR2
tumors_interest_kras_fgfr2 <- comp_kras_fgfr2$samples_with$KRAS.1[comp_kras_fgfr2$samples_with$KRAS.1 %in% comp_kras_fgfr2$samples_with$FGFR2.1]


# selecting only kras and fgfr genes and sinlge nucleotide variants
kras_fgfr2_data <- cesa_df |>
  filter(genes %in% c("KRAS", "FGFR2")) |>
  filter(variant_type == "snv") |>
  filter(Unique_Patient_Identifier %in% tumors_interest_kras_fgfr2)

# selecting unique patient identifiers who have co-occurring substitutions in KRAS and FGFR2
kras_fgfr2_data <- kras_fgfr2_data |>
  count(Unique_Patient_Identifier, genes) |>
  count(Unique_Patient_Identifier) |>
  filter(n > 1)

# effect size analysis for all kras and fgfr2 variants 
cesa_kras_fgfr2 <- ces_variant(cesa = cesa_kras_fgfr2, variants = kras_fgfr2_variants, run_name = "kras_fgfr2_effect_size")
kras_fgfr2_selection <- cesa_kras_fgfr2$selection$kras_fgfr2_effect_size



# find what variants of KRAS and FGFR2 exist in each of the tumors of interest 

# selecting tumors of interest and samples with KRAS and FGFR2 genes from the maf data
kras_fgfr2_maf <- cesa_kras_fgfr2$maf |>
  filter(Unique_Patient_Identifier %in% tumors_interest_kras_fgfr2) |>
  filter(genes %in% c("KRAS", "FGFR2"))

# combining the selection intensity for each variant with the tumors of interest maf data
kras_maf_selection <- left_join(kras_fgfr2_maf, kras_fgfr2_selection, by = c("top_consequence" = "variant_name")) |>
  relocate(Unique_Patient_Identifier, selection_intensity, genes, top_consequence) |>
  select(Unique_Patient_Identifier, top_consequence, selection_intensity) 









# epistasis plots 



# arrow plot with confidence intervals for PIK3CA vs PIK3R1

current_gene_pik3 <- "PIK3CA"

compare_pik3 <- epistasis_pik3|>
  mutate(variant_A = gsub("\\.1.*","",variant_A))|>
  mutate(variant_B = gsub("\\.1.*","",variant_B))|>
  filter(variant_A == current_gene_pik3 | variant_B == current_gene_pik3) 

compare_pik3 <- compare_pik3|>
  mutate(gene_of_interest = current_gene_pik3)|>
  mutate(other_gene = case_when(
    variant_A == current_gene_pik3 ~ variant_B,
    variant_B == current_gene_pik3 ~ variant_A))|>
  mutate(ces_GOI = case_when(
    variant_A == current_gene_pik3 ~ ces_A0,
    variant_B == current_gene_pik3 ~ ces_B0))|>
  mutate(ces_OTHER = case_when(
    variant_B == current_gene_pik3 ~ ces_A0,
    variant_A == current_gene_pik3 ~ ces_B0))|>
  mutate(ces_GOI_after_OTHER = case_when(
    variant_A == current_gene_pik3 ~ ces_A_on_B,
    variant_B == current_gene_pik3 ~ ces_B_on_A))|>
  mutate(ces_OTHER_after_GOI = case_when(
    variant_B == current_gene_pik3 ~ ces_A_on_B,
    variant_A == current_gene_pik3 ~ ces_B_on_A))|>
  mutate(joint_cov_samples_just_GOI = case_when(
    variant_B == current_gene_pik3 ~ nB0,
    variant_A == current_gene_pik3 ~ nA0))|>
  mutate(joint_cov_samples_just_OTHER = case_when(
    variant_A == current_gene_pik3 ~ nB0,
    variant_B == current_gene_pik3 ~ nA0))|>
  mutate(ci_low_95_ces_GOI = case_when(
    variant_A == current_gene_pik3 ~ ci_low_95_ces_A0,
    variant_B == current_gene_pik3 ~ ci_low_95_ces_B0))|>
  mutate(ci_high_95_ces_GOI = case_when(
    variant_A == current_gene_pik3 ~ ci_high_95_ces_A0,
    variant_B == current_gene_pik3 ~ ci_high_95_ces_B0))|>
  mutate(ci_low_95_ces_OTHER = case_when(
    variant_A == current_gene_pik3 ~ ci_low_95_ces_B0,
    variant_B == current_gene_pik3 ~ ci_low_95_ces_A0))|>
  mutate(ci_high_95_ces_OTHER = case_when(
    variant_A == current_gene_pik3 ~ ci_high_95_ces_B0,
    variant_B == current_gene_pik3 ~ ci_high_95_ces_A0))|>
  mutate(ci_low_95_ces_GOI_after_OTHER = case_when(
    variant_A == current_gene_pik3 ~ ci_low_95_ces_A_on_B,
    variant_B == current_gene_pik3 ~ ci_low_95_ces_B_on_A))|>
  mutate(ci_high_95_ces_GOI_after_OTHER = case_when(
    variant_A == current_gene_pik3 ~ ci_high_95_ces_A_on_B,
    variant_B == current_gene_pik3 ~ ci_high_95_ces_B_on_A))|>
  mutate(ci_low_95_ces_OTHER_after_GOI = case_when(
    variant_A == current_gene_pik3 ~ ci_low_95_ces_B_on_A,
    variant_B == current_gene_pik3 ~ ci_low_95_ces_A_on_B))|>
  mutate(ci_high_95_ces_OTHER_after_GOI = case_when(
    variant_A == current_gene_pik3 ~ ci_high_95_ces_B_on_A,
    variant_B == current_gene_pik3 ~ ci_high_95_ces_A_on_B))|>
  select(gene_of_interest, other_gene, ends_with("OTHER"),ends_with("GOI"), joint_cov_samples_just_GOI, joint_cov_samples_just_OTHER, nAB, n00, ci_low_95_ces_GOI, ci_high_95_ces_GOI, ci_low_95_ces_OTHER, ci_high_95_ces_OTHER, ci_low_95_ces_A_on_B, ci_high_95_ces_A_on_B, ci_low_95_ces_B_on_A, ci_high_95_ces_B_on_A)|>
  mutate(ci_low_95_ces_GOI = if_else(is.na(ci_low_95_ces_GOI), 0, ci_low_95_ces_GOI))|>
  mutate(ci_low_95_ces_OTHER = if_else(is.na(ci_low_95_ces_OTHER), 0, ci_low_95_ces_OTHER))|>
  mutate(ci_low_95_ces_GOI_after_OTHER = if_else(is.na(ci_low_95_ces_GOI_after_OTHER), 0, ci_low_95_ces_GOI_after_OTHER))|>
  mutate(ci_low_95_ces_OTHER_after_GOI = if_else(is.na(ci_low_95_ces_OTHER_after_GOI), 0, ci_low_95_ces_OTHER_after_GOI))|>
  mutate(ci_high_95_ces_GOI = if_else(is.na(ci_high_95_ces_GOI), 50000, ci_high_95_ces_GOI))|>
  mutate(ci_high_95_ces_OTHER = if_else(is.na(ci_high_95_ces_OTHER), 50000, ci_high_95_ces_OTHER))|>
  mutate(ci_high_95_ces_GOI_after_OTHER = if_else(is.na(ci_high_95_ces_GOI_after_OTHER), 50000, ci_high_95_ces_GOI_after_OTHER))|>
  mutate(ci_high_95_ces_OTHER_after_GOI = if_else(is.na(ci_high_95_ces_OTHER_after_GOI), 50000, ci_high_95_ces_OTHER_after_GOI))


compare_pik3$other_gene <- factor(compare_pik3$other_gene, levels = unique(compare_pik3$other_gene))













# new epistasis plot, legends outside

text_size <- 20
geom_text_size <- text_size * (5/14)


# 
# labels for pik3 plot
pik3ca_label <- "PIK3CA"
pik3ca_label_2 <- "PIK3CA in a background of PIK3R1"
pik3r1_label <- "PIK3R1"
pik3r1_label_2 <- "PIK3R1 in a background of PIK3CA"

# creating arrow plot for epistatic relationship between PIK3CA and PIK3R1
ci_plot_pik3 <- ggplot(compare_pik3) +
  geom_segment(data = compare_pik3, aes(x = as.numeric(other_gene)-0.1, xend = as.numeric(other_gene)-0.1, y = ces_GOI,
                                        yend = ces_GOI_after_OTHER), color="#008B8B", 
               arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill = "#008B8B", size = 1) +
  geom_segment(data = compare_pik3, aes(x = as.numeric(other_gene)+0.1, xend = as.numeric(other_gene)+0.1,
                                        y = ces_OTHER, yend = ces_OTHER_after_GOI), color="#C7AA82", 
               arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill = "#C7AA82", size = 1) +
  geom_point(data = compare_pik3, aes(x = as.numeric(other_gene)-0.1, y = ces_GOI), color = "#008B8B", size = 4) +
  geom_point(data = compare_pik3, aes(x = as.numeric(other_gene)+0.1, y = ces_OTHER), color="#C7AA82", size = 4) +
  geom_errorbar(aes(x = as.numeric(other_gene)-0.15, ymin = ci_low_95_ces_GOI, ymax = ci_high_95_ces_GOI), width = 0.05, color = "#008B8B") +
  geom_errorbar(aes(x = as.numeric(other_gene)-0.05, ymin = ci_low_95_ces_GOI_after_OTHER, ymax = ci_high_95_ces_GOI_after_OTHER), width = 0.05, color = "#008B8B") +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.15, ymin = ci_low_95_ces_OTHER_after_GOI, ymax = ci_high_95_ces_OTHER_after_GOI), width = 0.05, color = "#C7AA82") +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.05, ymin = ci_low_95_ces_OTHER, ymax = ci_high_95_ces_OTHER), width = 0.05, color = "#C7AA82") +
  coord_cartesian(xlim = c(0.75, 1.25), ylim = c(0, 4500)) +
  labs(x = NULL, y = "Cancer effect size") +
  theme_classic() +
  scale_y_continuous(expand = c(0.01, 0)) + 
  scale_x_continuous(breaks = c(0.9,1.1),labels = c("PIK3CA","PIK3R1")) +
  # annotate("rect", xmin = .1, xmax = 1.1, ymin = -1100, ymax = -200, alpha = 0.1) +
  # annotate("rect", xmin = 1.15, xmax = 2.1, ymin = -1100, ymax = -200, alpha = 0.1) +
  theme(text = element_text(size = text_size))




ci_plot_pik3_annot <- ggplot(compare_pik3) +
  annotate("text", x = 1.3, y = 4000, hjust = 0, label = pik3ca_label, size = geom_text_size) +
  annotate("point", x = 1.25, y = 4000, size = 3, color = "#008B8B", shape = 16) +
  annotate("text", x = 1.3, y = 3600, hjust = 0, label = pik3ca_label_2, size = geom_text_size) +
  annotate("point", x = 1.25, y = 3600, size = 3, color = "#008B8B", shape = 17) +
  annotate("text", x = 1.3, y = 3200, hjust = 0, label = pik3r1_label, size = geom_text_size) +
  annotate("point", x = 1.25, y = 3200, size = 3, color = "#C7AA82", shape = 16) +
  annotate("text", x = 1.3, y = 2800, hjust = 0, label = pik3r1_label_2, size = geom_text_size) +
  annotate("point", x = 1.25, y = 2800, size = 3, color = "#C7AA82", shape = 17) +
  annotate("rect", xmin = 1.2, xmax = 2.2, ymin = 2600, ymax = 4100, alpha = 0) +
  scale_y_continuous(expand = c(0.01, 0)) + 
  # annotate("rect", xmin = .1, xmax = 1.1, ymin = -1100, ymax = -200, alpha = 0.1) +
  # annotate("rect", xmin = 1.15, xmax = 2.1, ymin = -1100, ymax = -200, alpha = 0.1) +
  theme_void() + 
  theme(text = element_text(size = text_size))

ci_plot_pik3 + ci_plot_pik3_annot




# arrow plot with confidence intervals for KRAS vs FGFR2

current_gene_kras <- "KRAS"

compare_kras <- epistasis_kras_fgfr2|>
  mutate(variant_A = gsub("\\.1.*","",variant_A))|>
  mutate(variant_B = gsub("\\.1.*","",variant_B))|>
  filter(variant_A == current_gene_kras | variant_B == current_gene_kras) 

compare_kras <- compare_kras|>
  mutate(gene_of_interest = current_gene_kras)|>
  mutate(other_gene = case_when(
    variant_A == current_gene_kras ~ variant_B,
    variant_B == current_gene_kras ~ variant_A))|>
  mutate(ces_GOI = case_when(
    variant_A == current_gene_kras ~ ces_A0,
    variant_B == current_gene_kras ~ ces_B0))|>
  mutate(ces_OTHER = case_when(
    variant_B == current_gene_kras ~ ces_A0,
    variant_A == current_gene_kras ~ ces_B0))|>
  mutate(ces_GOI_after_OTHER = case_when(
    variant_A == current_gene_kras ~ ces_A_on_B,
    variant_B == current_gene_kras ~ ces_B_on_A))|>
  mutate(ces_OTHER_after_GOI = case_when(
    variant_B == current_gene_kras ~ ces_A_on_B,
    variant_A == current_gene_kras ~ ces_B_on_A))|>
  mutate(joint_cov_samples_just_GOI = case_when(
    variant_B == current_gene_kras ~ nB0,
    variant_A == current_gene_kras ~ nA0))|>
  mutate(joint_cov_samples_just_OTHER = case_when(
    variant_A == current_gene_kras ~ nB0,
    variant_B == current_gene_kras ~ nA0))|>
  mutate(ci_low_95_ces_GOI = case_when(
    variant_A == current_gene_kras ~ ci_low_95_ces_A0,
    variant_B == current_gene_kras ~ ci_low_95_ces_B0))|>
  mutate(ci_high_95_ces_GOI = case_when(
    variant_A == current_gene_kras ~ ci_high_95_ces_A0,
    variant_B == current_gene_kras ~ ci_high_95_ces_B0))|>
  mutate(ci_low_95_ces_OTHER = case_when(
    variant_A == current_gene_kras ~ ci_low_95_ces_B0,
    variant_B == current_gene_kras ~ ci_low_95_ces_A0))|>
  mutate(ci_high_95_ces_OTHER = case_when(
    variant_A == current_gene_kras ~ ci_high_95_ces_B0,
    variant_B == current_gene_kras ~ ci_high_95_ces_A0))|>
  mutate(ci_low_95_ces_GOI_after_OTHER = case_when(
    variant_A == current_gene_kras ~ ci_low_95_ces_A_on_B,
    variant_B == current_gene_kras ~ ci_low_95_ces_B_on_A))|>
  mutate(ci_high_95_ces_GOI_after_OTHER = case_when(
    variant_A == current_gene_kras ~ ci_high_95_ces_A_on_B,
    variant_B == current_gene_kras ~ ci_high_95_ces_B_on_A))|>
  mutate(ci_low_95_ces_OTHER_after_GOI = case_when(
    variant_A == current_gene_kras ~ ci_low_95_ces_B_on_A,
    variant_B == current_gene_kras ~ ci_low_95_ces_A_on_B))|>
  mutate(ci_high_95_ces_OTHER_after_GOI = case_when(
    variant_A == current_gene_kras ~ ci_high_95_ces_B_on_A,
    variant_B == current_gene_kras ~ ci_high_95_ces_A_on_B))|>
  select(gene_of_interest, other_gene, ends_with("OTHER"),ends_with("GOI"), joint_cov_samples_just_GOI, joint_cov_samples_just_OTHER, nAB, n00, ci_low_95_ces_GOI, ci_high_95_ces_GOI, ci_low_95_ces_OTHER, ci_high_95_ces_OTHER, ci_low_95_ces_A_on_B, ci_high_95_ces_A_on_B, ci_low_95_ces_B_on_A, ci_high_95_ces_B_on_A)|>
  mutate(ci_low_95_ces_GOI = if_else(is.na(ci_low_95_ces_GOI), 0, ci_low_95_ces_GOI))|>
  mutate(ci_low_95_ces_OTHER = if_else(is.na(ci_low_95_ces_OTHER), 0, ci_low_95_ces_OTHER))|>
  mutate(ci_low_95_ces_GOI_after_OTHER = if_else(is.na(ci_low_95_ces_GOI_after_OTHER), 0, ci_low_95_ces_GOI_after_OTHER))|>
  mutate(ci_low_95_ces_OTHER_after_GOI = if_else(is.na(ci_low_95_ces_OTHER_after_GOI), 0, ci_low_95_ces_OTHER_after_GOI))|>
  mutate(ci_high_95_ces_GOI = if_else(is.na(ci_high_95_ces_GOI), 50000, ci_high_95_ces_GOI))|>
  mutate(ci_high_95_ces_OTHER = if_else(is.na(ci_high_95_ces_OTHER), 50000, ci_high_95_ces_OTHER))|>
  mutate(ci_high_95_ces_GOI_after_OTHER = if_else(is.na(ci_high_95_ces_GOI_after_OTHER), 50000, ci_high_95_ces_GOI_after_OTHER))|>
  mutate(ci_high_95_ces_OTHER_after_GOI = if_else(is.na(ci_high_95_ces_OTHER_after_GOI), 50000, ci_high_95_ces_OTHER_after_GOI))


compare_kras$other_gene <- factor(compare_kras$other_gene, levels = unique(compare_kras$other_gene))

kras_label <- "KRAS"
kras_label_2 <- "KRAS in a background of FGFR2"
fgfr2_label <- "FGFR2"
fgfr2_label_2 <- "FGFR2 in a background of KRAS"

ci_plot_kras_fgfr2 <- ggplot(compare_kras) +
  geom_segment(data = compare_kras, aes(x=as.numeric(other_gene)-0.1,xend=as.numeric(other_gene)-0.1,y=ces_GOI,yend=ces_GOI_after_OTHER), color="#008B8B", arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill="#008B8B", size = 1) +
  geom_segment(data = compare_kras, aes(x=as.numeric(other_gene)+0.1,xend=as.numeric(other_gene)+0.1,y=ces_OTHER,yend=ces_OTHER_after_GOI), color="#C7AA82", arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill="#C7AA82", size = 1) +
  geom_point(data = compare_kras, aes(x=as.numeric(other_gene)-0.1,y=ces_GOI), color="#008B8B", size=4) +
  geom_point(data = compare_kras, aes(x=as.numeric(other_gene)+0.1,y=ces_OTHER), color="#C7AA82", size=4) +
  geom_errorbar(aes(x=as.numeric(other_gene)-0.15, ymin = ci_low_95_ces_GOI, ymax = ci_high_95_ces_GOI), width =0.05, color="#008B8B") +
  geom_errorbar(aes(x=as.numeric(other_gene)-0.05, ymin = ci_low_95_ces_GOI_after_OTHER, ymax = ci_high_95_ces_GOI_after_OTHER), width =0.05, color="#008B8B") +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.15, ymin = ci_low_95_ces_OTHER_after_GOI, ymax = ci_high_95_ces_OTHER_after_GOI), width =0.05, color="#C7AA82") +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.05, ymin = ci_low_95_ces_OTHER, ymax = ci_high_95_ces_OTHER), width =0.05, color="#C7AA82") +
  coord_cartesian(ylim = c(0, 4500), xlim = c(0.75, 1.25)) +
  labs(x = NULL, y = "Cancer effect size") +
  scale_x_continuous(breaks = c(0.9,1.1),labels = c("KRAS","FGFR2")) +
  theme_classic() +
  # annotate("rect", xmin = 1.15, xmax = 2, ymin = 2700, ymax = 4100, alpha = 0.1) +
  scale_y_continuous(expand = c(0.01, 0)) + 
  theme(text = element_text(size = text_size))

ci_plot_kras_fgfr2

ci_plot_kras_fgfr2_annot <- ggplot(compare_kras) + 
  annotate("text", x = 1.3, y = 4000, hjust = 0, label = kras_label, size = geom_text_size) +
  annotate("point", x = 1.25, y = 4000, size = 3, color = "#008B8B", shape = 16) +
  annotate("text", x = 1.3, y = 3600, hjust = 0, label = kras_label_2, size = geom_text_size) +
  annotate("point", x = 1.25, y = 3600, size = 3, color = "#008B8B", shape = 17) +
  annotate("text", x = 1.3, y = 3200, hjust = 0, label = fgfr2_label, size = geom_text_size) +
  annotate("point", x = 1.25, y = 3200, size = 3, color = "#C7AA82", shape = 16) +
  annotate("text", x = 1.3, y = 2800, hjust = 0, label = fgfr2_label_2, size = geom_text_size) +
  annotate("point", x = 1.25, y = 2800, size = 3, color = "#C7AA82", shape = 17) +
  annotate("rect", xmin = 1.2, xmax = 2.2, ymin = 2600, ymax = 4100, alpha = 0) +
  theme_void() + 
  theme(text = element_text(size = text_size))


ci_plot_kras_fgfr2 + ci_plot_kras_fgfr2_annot



ci_plot_pik3_annot <- ci_plot_pik3_annot + plot_layout(tag_level = 'new')

ci_plot_kras_fgfr2_annot <- ci_plot_kras_fgfr2_annot + plot_layout(tag_level = 'new')



both_epi_plot <- 
  ci_plot_pik3 + 
  ci_plot_pik3_annot + 
  ci_plot_kras_fgfr2 + 
  ci_plot_kras_fgfr2_annot  + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(nrow = 2) & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(hjust = 0, vjust = 1)) 


ggsave(filename = "figures/figure_3_epistasis.png", width =12, height = 7)


