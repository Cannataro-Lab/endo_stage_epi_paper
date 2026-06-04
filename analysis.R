# Full analysis completed within Summers et al. 
# Driver genes in endometrial hyperplasia and carcinoma exhibit stage-specific mutation rates, changing selection intensities, and antagonistic epistasis

# Code collaboratively built by Mary Summers, Sem Asmelash, J. Nick Fisk, Jeffrey D. Mandell, and Vincent L. Cannataro 
# Code updated to cancereffectsizeR v2.10.2 by Kira A. Glasmacher





# load cancer effect size and necessary packages ----

library(cancereffectsizeR) # v2.10.2 https://townsend-lab-yale.github.io/cancereffectsizeR/
library(ces.refset.hg19) #v1.1.3
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
                   figo_stage) %>%
distinct()

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
# removing samples where column Problem is not equal to NA
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

# removing samples where column Problem is not equal to NA
maf_ucec <- maf_ucec[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
maf_ucec <- maf_ucec[germline_variant_site == F]

# keeping only samples that do not occur in repetitive regions 
maf_ucec <- maf_ucec[(repetitive_region == F | cosmic_site_tier %in% 1:3)]


# combining ucec clinical and maf data
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
cptac.maf <- read_tsv("input_data/cptac.maf", comment = "#")
manifest <- read_tsv("input_data/gdc_manifest.2026-06-01.155053.txt") # using this set of samples will ensure consistency with our results even with newer data releases
cptac_clinical <- read_tsv("input_data/clinical.project-cptac-3.2026-06-01/clinical.tsv")

cptac_clinical <- cptac_clinical %>%
  filter(`cases.primary_site` == "Uterus, NOS") %>%
  select(`cases.case_id`, `diagnoses.ajcc_pathologic_stage`) %>%
  distinct()

# selecting data from CPTAC maf file with tumor sample barcodes matching the samples from manifest
maf_endo_cptac <- cptac.maf |> filter(source_file_id %in% manifest$id)


# joining clinical and sample data sets by case id
maf_endo_cptac <- left_join(x = maf_endo_cptac, y = cptac_clinical, by = c("case_id" = "cases.case_id"))


# pre load maf of clinical and sample data 
maf_cptac <- preload_maf(maf = maf_endo_cptac, 
                         refset = "ces.refset.hg19", 
                         chain_file = "input_data/hg38ToHg19.over.chain", 
                         keep_extra_columns = c("case_id", "diagnoses.ajcc_pathologic_stage"))


# removing samples where column Problem is not equal to NA
maf_cptac <- maf_cptac[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
maf_cptac <- maf_cptac[germline_variant_site == F]

# keeping only samples that do not occur in repetitive regions 
maf_cptac <- maf_cptac[(repetitive_region == F | cosmic_site_tier %in% 1:3)]

# adding stage column and specifying Stage 1 and Rest of Stages
maf_cptac <- maf_cptac |> 
  filter(!is.na(diagnoses.ajcc_pathologic_stage)) |> 
  mutate(Stage = case_when(
    diagnoses.ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage1",
    TRUE ~ "RestOfStages"))





# stage specific analysis of selection intensity ----

# initiate analysis
cesa <- CESAnalysis(refset = "ces.refset.hg19")

# load in data from CPTAC, TCGA and Li projects
cesa <- load_maf(cesa = cesa, maf = maf_cptac, sample_data_cols = "Stage")
cesa <- load_maf(cesa = cesa, maf = VCF_all_maf_data, sample_data_cols = "Stage")
cesa <- load_maf(cesa = cesa, maf = maf_ucec, sample_data_cols = "Stage")


# estimating mutation rates
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "UCEC", samples = cesa$samples[Stage == "AH"],save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "UCEC", samples = cesa$samples[Stage == "Stage1"], save_all_dndscv_output = T)
# cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "UCEC", samples = cesa$samples[Stage == "RestOfStages"], save_all_dndscv_output = T)


# selecting genes of interest
selected_genes <- c("KRAS", "PTEN", "ARID1A", "CTCF", "CTNNB1", "PIK3CA", "CHD4", "FGFR2", "PIK3R1")


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

# check for monotonic increase of expected baseline mutation flux
table(mut_rate_df$cancer_greater)

# defining rate 1 and rate 2 as mutation rates for hyperplasia and stage 1
rate_1 <- mut_rate_df|>
  select(gene, hyp_mu)
rate_2 <- mut_rate_df|>
  select(gene, cancer1_mu)

# change in mutation rate across stages
mut_rates_for_p <- mut_rate_df %>% 
  select(gene, hyp_mu, cancer1_mu) %>% 
  mutate(p_1 = hyp_mu / cancer1_mu) %>% 
  mutate(p_2 = 1 - p_1)

# saving "last" gene mutation rates into separate data frame, "last" rates meaning from last stage cancer1_mu
set_cancer_rates <- mut_rate_df %>%
  select(gene, rate = cancer1_mu) %>%
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
source("R/modified_ces_variant.R")

index_by_state = list()
name_by_state = list()
ordering_col = 'Stage'
ordering = c('AH', 'Stage1')
if(is.null(names(ordering))) {
  if (length(unlist(ordering)) == length(ordering)) {
    names(ordering) = unlist(ordering)
  } else {
    names(ordering) = 1:length(ordering)
  }
}
for (i in 1:length(ordering)) {
  for (j in 1:length(ordering[[i]])) {
    index_by_state[[ordering[[i]][j]]] = i
    name_by_state[[ordering[[i]][j]]] = names(ordering)[i]
  }
}
samples = cancereffectsizeR:::select_samples(cesa_samples_by_groups, samples=cesa_samples_by_groups$samples[Stage != "RestOfStages"])
sample_index_table = samples[, .(Unique_Patient_Identifier = Unique_Patient_Identifier,
                                 group_index = unlist(index_by_state[samples[[ordering_col]]]), 
                                 group_name = unlist(name_by_state[samples[[ordering_col]]]))]

# Calculate non-step-specific selection intensity for LRT----
source("R/new_sequential_lik_const_gamma.R")

for(comp_ind in 1:length(compound)){
  this_comp <- compound[comp_ind, ]
  this_gene <- unlist(unique(this_comp$snv_info$genes))[1]
  these_props <- mut_rates_for_p[mut_rates_for_p$gene %in% this_gene, c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  if(length(this_gene) != 1){
    this_gene <- unlist(str_split(this_gene[1], "\\."))
    this_gene <- this_gene[1]
  }
  
  cat("Running gene:", this_gene, "\n")
  
  if (length(these_props) != 2 || any(is.na(these_props))) {
    stop(paste("Missing or invalid mutation rates for", this_gene))
  }
  
  cesa_samples_by_groups <- modified_ces_variant(cesa = cesa_samples_by_groups,
                                                 samples=cesa_samples_by_groups$samples[Stage != "RestOfStages"], 
                                                 variants = this_comp, 
                                                 model = sequential_lik_dev_const_gamma, 
                                                 lik_args = list(sample_index = sample_index_table, 
                                                                 sequential_mut_prop = these_props), 
                                                 optimizer_args = list(method = 'L-BFGS-B', 
                                                                       lower = 1e-3, 
                                                                       upper = 1e9),
                                                 return_fit = TRUE,
                                                 run_name = paste0(this_gene, "_const_gamma"),
                                                 conf = 0.95)
}

selection_results_constant <- rbind(cesa_samples_by_groups@selection_results$KRAS_const_gamma,
                                    cesa_samples_by_groups@selection_results$PTEN_const_gamma,
                                    cesa_samples_by_groups@selection_results$ARID1A_const_gamma,
                                    cesa_samples_by_groups@selection_results$CTCF_const_gamma,
                                    cesa_samples_by_groups@selection_results$CTNNB1_const_gamma,
                                    cesa_samples_by_groups@selection_results$PIK3CA_const_gamma,
                                    cesa_samples_by_groups@selection_results$CHD4_const_gamma,
                                    cesa_samples_by_groups@selection_results$FGFR2_const_gamma,
                                    cesa_samples_by_groups@selection_results$PIK3R1_const_gamma)

# Step-specific gamma estimation ----
for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))[1]
  these_props <- mut_rates_for_p[mut_rates_for_p$gene %in% this_gene, c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  if(length(this_gene) != 1){
    this_gene <- unlist(str_split(this_gene[1], "\\."))
    this_gene <- this_gene[1]
  }
  
  cat("Running gene:", this_gene, "\n")
  # print(this_comp)
  # print(these_props)
  
  if (length(these_props) != 2 || any(is.na(these_props))) {
    stop(paste("Missing or invalid mutation rates for", this_gene))
  }
  
  cesa_samples_by_groups <- modified_ces_variant(cesa = cesa_samples_by_groups,
                                                 samples=cesa_samples_by_groups$samples[Stage != "RestOfStages"], 
                                                 variants = this_comp, 
                                                 model = sequential_lik_dev, 
                                                 lik_args = list(sample_index = sample_index_table, 
                                                                 sequential_mut_prop = these_props), 
                                                 optimizer_args = list(method = 'L-BFGS-B', 
                                                                       lower = 1e-3, 
                                                                       upper = 1e9),
                                                 return_fit = TRUE,
                                                 run_name = this_gene,
                                                 conf = 0.95)
  }

selection_results_step <- rbind(cesa_samples_by_groups@selection_results$KRAS,
                                cesa_samples_by_groups@selection_results$PTEN,
                                cesa_samples_by_groups@selection_results$ARID1A,
                                cesa_samples_by_groups@selection_results$CTCF,
                                cesa_samples_by_groups@selection_results$CTNNB1,
                                cesa_samples_by_groups@selection_results$PIK3CA,
                                cesa_samples_by_groups@selection_results$CHD4,
                                cesa_samples_by_groups@selection_results$FGFR2,
                                cesa_samples_by_groups@selection_results$PIK3R1)

# creating selection plots with stage specific model ----

scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}

# selecting necessary data

selection_data <- selection_results_step


# computing Likelihood Root 95% CIs
ci_results <- list()

for (gene in selected_genes) {
  fit_list <- attr(cesa_samples_by_groups@selection_results[[gene]], "fit")
  fit <- fit_list[[1]]
  
  lik_fn <- get(paste0("lik_fn_", gene))
  
  ci <- cancereffectsizeR:::univariate_si_conf_ints(
    fit = fit,
    lik_fn = lik_fn,
    min_si = 0.001,
    max_si = 1e9,
    conf = 0.95
  )
  
  coef_vals <- coef(fit)
  
  ci_results[[gene]] <- data.frame(
    gene = gene,
    si_AH = coef_vals["si_AH"],
    si_Stage1 = coef_vals["si_Stage1"],
    ci_AH_low = ifelse(is.null(ci$ci_low_95_si_AH), NA, ci$ci_low_95_si_AH),
    ci_AH_high = ifelse(is.null(ci$ci_high_95_si_AH), NA, ci$ci_high_95_si_AH),
    ci_Stage1_low = ifelse(is.null(ci$ci_low_95_si_Stage1), NA, ci$ci_low_95_si_Stage1),
    ci_Stage1_high = ifelse(is.null(ci$ci_high_95_si_Stage1), NA, ci$ci_high_95_si_Stage1)
  )
}

ci_df <- do.call(rbind, ci_results)

ci_df <- ci_df %>%
  mutate(ci_AH_low = ifelse(ci_AH_low < 0.001 | is.na(ci_AH_low), 0.001, ci_AH_low),
         ci_Stage1_low = ifelse(ci_Stage1_low < 0.001 | is.na(ci_Stage1_low), 0.001, ci_Stage1_low))

ci_df$gene <- factor(ci_df$gene, levels = unique(ci_df$gene))

# Calculate likelihood ratio ----
loglik_step <- selection_results_step$loglikelihood
loglik_simple <- selection_results_constant$loglikelihood

loglik_df <- data.frame(
  gene = selected_genes,
  loglik_step = loglik_step,
  loglik_simple = loglik_simple)

loglik_df <- loglik_df %>%
  mutate(
    loglik_step = ifelse(loglik_step > 1e5, NA, loglik_step), # just making sure that the loglikelihoods are realistic
    loglik_simple = ifelse(loglik_simple > 1e5, NA, loglik_simple), # not some large positive number (happens when convergence failed)
    LRT_stat = -2 * (loglik_simple - loglik_step),
    p_value = pchisq(LRT_stat, df = 1, lower.tail = FALSE), 
    p_less_0.05 = ifelse(p_value < 0.05, TRUE, FALSE))
                    
# arrange genes in plotting df by relative difference in selection intensity between stage 1 and AH, calculated as (si_stage1 - si_AH) / si_AH
ci_df <- ci_df %>% 
  mutate(relative_diff = (si_Stage1 - si_AH) / si_AH) %>% 
  arrange(desc(relative_diff)) %>% 
  # make genes factors with levels determined by relative difference
  mutate(gene = factor(gene, levels = gene))

# reformatting data set
plotting_df <- ci_df |> 
  select(variant_name = gene, starts_with("si"), starts_with("ci")) |>
  pivot_longer(cols = -variant_name, names_to = "data_type") |>
  mutate(stage = ifelse(str_detect(data_type, "AH"), "AH", "Stage1")) |>
  mutate(si_or_ci = case_when(str_detect(data_type, "si") ~"si", 
                              str_detect(data_type, "low") ~"ci_low_95",
                              str_detect(data_type, "high") ~"ci_high_95")) |>
  
  mutate (value = case_when(is.na(value)~0, TRUE~value))

# pivoting data set to create columns for gene, stage, si, and CIs
plotting_df <- plotting_df|> 
  select(-data_type) |>
  pivot_wider(values_from = value, names_from = si_or_ci)

# defining stages to be plotted
plotting_df$stage <- factor(plotting_df$stage, levels = c("AH","Stage1"))

# add LRT result to plotting data frame and arrange by relative diff in selection for plotting
plotting_df <- plotting_df %>% 
  left_join(loglik_df %>% select(gene, p_less_0.05), by = c("variant_name" = "gene")) %>%
  mutate(variant_name = factor(variant_name, levels = ci_df$gene))

# plotting results with LRT significance indicated by shape of points
plotting_df <- plotting_df %>%
  mutate(facet_label = case_when(p_less_0.05 == TRUE ~ paste0(variant_name, "*", sep = ""),
                                 p_less_0.05 == FALSE ~ variant_name)) %>%
  mutate(facet_label = factor(facet_label, levels = c("PIK3R1*", "CTCF", "CHD4", "ARID1A",
                                                      "PIK3CA", "FGFR2", "CTNNB1", "PTEN", "KRAS")))

font_size <- 18

plotting_df |>
  ggplot(aes(x = stage, y = si, color = stage)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .3) +
  facet_wrap(~facet_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("dark green", "purple"),
                     labels = c("Normal through AH", "AH through Stage I")) +
  theme_bw() +
  labs(x = "Progression", y = "Cancer effect size", color = "Progression step") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = font_size),
        text = element_text(size = font_size),
        strip.text = element_text(face = "italic", size = font_size),
        axis.title = element_text(size = font_size),
        axis.text = element_text(size = font_size)) +
  scale_y_continuous(labels = scientific) +
  scale_x_discrete(labels = c("Step 1", "Step 2")) +
  expand_limits(y = 0)


ggsave(filename = "figures/figure_2_stages_selection.png", width = 10, height = 8)


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

ggplot(prevalence_df, aes(x = stages, y = prop_tumors_with)) +
  geom_point() +
  geom_segment(data = prevalence_df_w, aes(x = 1, xend = 2,
                                           y = `prop_tumors_with_Atypical hyperplasia`,
                                           yend = `prop_tumors_with_Stage-1 EC`)) +
  geom_text_repel(data = prevalence_df_w, aes(x = 1, y = `prop_tumors_with_Atypical hyperplasia`, label = gene_name),
                  direction = "y", nudge_x = -0.01, hjust = 1, size = geom_font_size, fontface = "italic") +
  scale_color_viridis_d() +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Atypical\nhyperplasia", "Stage I\nendometrial carcinoma")) +
  theme(text = element_text(size = font_size),
        axis.text.x = element_text(size = font_size),
        axis.text.y = element_text(size = font_size),
        axis.title.y = element_text(size = font_size)) +
  labs(y = "Proportions of tissue samples with\ndriver mutations", x = NULL) ->
  prev_plot


ggsave(filename = "figures/figure_1_prevalence.png",plot = prev_plot, height = 5,width = 4)

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
  select(gene, rate = mut_rate) %>%
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


shared_tcga_pik3 <- grep(x = intersect(comp_pik3$samples_with$PIK3CA, comp_pik3$samples_with$PIK3R1), pattern = "TCGA",value = T)

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
tumors_of_interest_pik3 <- comp_pik3$samples_with$PIK3CA[comp_pik3$samples_with$PIK3CA %in% comp_pik3$samples_with$PIK3R1]


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
tumors_interest_kras_fgfr2 <- comp_kras_fgfr2$samples_with$KRAS[comp_kras_fgfr2$samples_with$KRAS %in% comp_kras_fgfr2$samples_with$FGFR2]


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
  geom_segment(data = compare_pik3, aes(x = as.numeric(other_gene)-0.15, xend = as.numeric(other_gene)-0.05, y = ces_GOI,
                                        yend = ces_GOI_after_OTHER), color="#008B8B", 
               arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill = "#008B8B", linewidth = 1) +
  geom_segment(data = compare_pik3, aes(x = as.numeric(other_gene)+0.05, xend = as.numeric(other_gene)+0.15,
                                        y = ces_OTHER, yend = ces_OTHER_after_GOI), color="#C7AA82", 
               arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill = "#C7AA82", linewidth = 1) +
  geom_point(data = compare_pik3, aes(x = as.numeric(other_gene)-0.15, y = ces_GOI), color = "#008B8B", size = 4) +
  geom_point(data = compare_pik3, aes(x = as.numeric(other_gene)+0.05, y = ces_OTHER), color="#C7AA82", size = 4) +
  geom_errorbar(aes(x = as.numeric(other_gene)-0.15, ymin = ci_low_95_ces_GOI, ymax = ci_high_95_ces_GOI), width = 0.02, color = "#008B8B", alpha = 0.6) +
  geom_errorbar(aes(x = as.numeric(other_gene)-0.05, ymin = ci_low_95_ces_GOI_after_OTHER, ymax = ci_high_95_ces_GOI_after_OTHER), width = 0.02, color = "#008B8B", alpha = 0.6) +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.15, ymin = ci_low_95_ces_OTHER_after_GOI, ymax = ci_high_95_ces_OTHER_after_GOI), width = 0.02, color = "#C7AA82", alpha = 0.6) +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.05, ymin = ci_low_95_ces_OTHER, ymax = ci_high_95_ces_OTHER), width = 0.02, color = "#C7AA82", alpha = 0.6) +
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
  geom_segment(data = compare_kras, aes(x=as.numeric(other_gene)-0.15,xend=as.numeric(other_gene)-0.05,y=ces_GOI,yend=ces_GOI_after_OTHER), color="#008B8B", arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill="#008B8B", linewidth = 1) +
  geom_segment(data = compare_kras, aes(x=as.numeric(other_gene)+0.05,xend=as.numeric(other_gene)+0.15,y=ces_OTHER,yend=ces_OTHER_after_GOI), color="#C7AA82", arrow = arrow(length = unit(0.1, "inches"), type ="closed"), arrow.fill="#C7AA82", linewidth = 1) +
  geom_point(data = compare_kras, aes(x=as.numeric(other_gene)-0.15,y=ces_GOI), color="#008B8B", size=4) +
  geom_point(data = compare_kras, aes(x=as.numeric(other_gene)+0.05,y=ces_OTHER), color="#C7AA82", size=4) +
  geom_errorbar(aes(x=as.numeric(other_gene)-0.15, ymin = ci_low_95_ces_GOI, ymax = ci_high_95_ces_GOI), width =0.02, color="#008B8B", alpha = 0.6) +
  geom_errorbar(aes(x=as.numeric(other_gene)-0.05, ymin = ci_low_95_ces_GOI_after_OTHER, ymax = ci_high_95_ces_GOI_after_OTHER), width =0.02, color="#008B8B", alpha = 0.6) +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.15, ymin = ci_low_95_ces_OTHER_after_GOI, ymax = ci_high_95_ces_OTHER_after_GOI), width =0.02, color="#C7AA82", alpha = 0.6) +
  geom_errorbar(aes(x=as.numeric(other_gene)+0.05, ymin = ci_low_95_ces_OTHER, ymax = ci_high_95_ces_OTHER), width =0.02, color="#C7AA82", alpha = 0.6) +
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


library(patchwork)

both_epi_plot <- 
  ci_plot_pik3 + 
  ci_plot_pik3_annot + 
  ci_plot_kras_fgfr2 + 
  ci_plot_kras_fgfr2_annot +
  plot_layout(nrow = 2) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag.position = c(0, 1),
      plot.tag = element_text(hjust = 0, vjust = 1)
    )
  )


ggsave(both_epi_plot, filename = "figures/figure_3_epistasis.png", width =13, height = 7)

