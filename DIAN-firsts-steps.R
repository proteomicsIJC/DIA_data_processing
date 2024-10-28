###########################################
### R workflow for DIAN data processing ###
###########################################

#### Libraries, functions and data importation
### libraries and WD setting
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(ggplot2)
library(ggfortify)
library(limma)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(janitor)
library(dplyr)
library(tidyr)
library(mice)
library(reshape2)
library(openxlsx)
library(ggrepel)
library(readr) 
library(Peptides)

### Directory preparation
# Check if folders exist, if NO create them
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
file.remove(file.path(wd,"./results/used_parameters.txt"))
file.create(file.path(wd, "./results/used_parameters.txt"))

### functions importation
source("./functions/columns_checker.R")
source("./functions/log2_to_pattern_label_free.R")
source("./functions/make_all_contrasts.R")
source("./functions/remove_batch.R")
source("./functions/remove_samp.R")
source("./functions/remove_contaminants.R")
source("./functions/tim.R")
source("./functions/zero_to_na_label_free.R")
source("./functions/sample_combination.R")

#### Data importation
### Expression matrix
raw_data <- as.data.frame(readr::read_tsv("./raw_data/report.pg_matrix.tsv", locale = readr::locale(encoding = "latin1")))
colnames(raw_data)[1] <- "Protein.Group"
colnames(raw_data)[2] <- "Accession"
colnames(raw_data)[3] <- "Protein.Name"
colnames(raw_data)[4] <- "Gene.names"
colnames(raw_data)[5] <- "Protein.Description"

### Peptide information
peptide_info <- readr::read_tsv("./raw_data/report.pr_matrix.tsv")

### MetaData
meta_data <- as.data.frame(readr::read_tsv("./raw_data/meta_data.tsv", locale = readr::locale(encoding = "latin1")))

### Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

#### Raw Data cleaning 
### raw_data sample name
colnames(raw_data) <- sapply(X = colnames(raw_data), 
                             FUN = function(x) {
                               if (x %in% meta_data$file_name){
                                 meta_data$sample_name[match(x, meta_data$file_name)]}
                               else {
                                 x
                               }
                               })

### Annotation completion with Uniprot.ws
## Accession_1 and Gene.names.1 generation
raw_data$Accession_1 <- sapply(X = raw_data$Accession, FUN = function(x) {unlist(strsplit(x, split = ";"))[1]})
raw_data$Gene.names_1 <- sapply(X = raw_data$Gene.names, FUN = function(x) {unlist(strsplit(x, split = ";"))[1]})

raw_data <- raw_data %>% 
  relocate(Accession_1, .after = Accession) %>%
  relocate(Gene.names_1, .after = Gene.names) %>% 
  relocate(Protein.Name, .after = Protein.Description) %>% 
  ungroup()

## Uniprot.ws petition
id_mapping <- UniProt.ws::mapUniProt("UniProtKB_AC-ID", "UniProtKB", query = raw_data$Accession_1)

raw_data <- raw_data %>% 
  group_by(Accession) %>% 
  mutate(Protein.Description = ifelse(Accession_1 %in% id_mapping$From,
                                        id_mapping$Protein.names[match(Accession_1, id_mapping$From)],NA))
raw_data <- as.data.frame(raw_data)
dian_data <- raw_data

openxlsx::write.xlsx(dian_data, "./results/Raw_data.xlsx")



#### Peptide based analysis
## Missed cleavages
unique_peptides <- unique(peptide_info$Stripped.Sequence)
cleavage_report <- tibble::tibble(
  # The peptides
  peptide = unique_peptides,
  
  # Count Ks and Rs
  k_or_r = sapply(X = unique_peptides, FUN = function(x) {
    length(grep(x = unlist(strsplit(x, split = ""),use.names = F), pattern = "K|R"))  
  }),
  
  # Count Ks and Rs that happen at the end
  final_k_or_r = sapply(X = unique_peptides, FUN = function(x) {
    final_letter <- unlist(strsplit(x, split = ""), use.names = F)
    final_letter <- final_letter[length(final_letter)]
    length(grep(x = final_letter, pattern = "K|R"))
  }),
  # Count KPs and RPs
  kp_or_rp = sapply(X = unique_peptides, FUN = function(x) {
    stringr::str_count(string = x, pattern = "RP|KP")
  })
)

cleavage_report <- cleavage_report %>% 
  mutate(missed_cleavages = k_or_r - final_k_or_r - kp_or_rp)

cleavage_report_per <- cleavage_report %>% 
  subset(select = c(missed_cleavages)) %>% 
  mutate(missed_cleavages_per = as.character(missed_cleavages)) %>% 
  group_by(missed_cleavages_per) %>% 
  mutate(n_of_cleavages = n()) %>% 
  ungroup() %>%
  mutate(per_of_cleavages = n_of_cleavages/length(unique_peptides)*100) %>% 
  subset(select = c(missed_cleavages,per_of_cleavages)) %>% 
  mutate(missed_cleavages = as.character(missed_cleavages))
cleavage_report_per <- cleavage_report_per %>% 
  distinct()

# Plot
missed_cleavages <- ggplot(data = cleavage_report_per, mapping = aes(x = missed_cleavages, y = per_of_cleavages))+
  geom_bar(stat = "identity", mapping = aes(fill = missed_cleavages)) +
  scale_fill_manual(values = c("0" = "skyblue",
                               "1" = "orange",
                               "2" = "red",
                               "3" = "black"))+
  geom_text(data = cleavage_report_per,
            aes(y = per_of_cleavages + 5, x = missed_cleavages, 
                label = paste0(round(x = per_of_cleavages, digits = 2), "%")), size = 8)+
  theme_bw()+
  ggtitle("NAs density per protein groups")+
  xlab("# Missed cleavages")+
  ylab("Peptide %")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))

png(filename = "./plots/0_missed_cleavages.png", width = 1400, height = 1400)
missed_cleavages
dev.off()

## Hydrophobicity
hydrophobicity_report <- tibble::tibble(
  # The peptides
  peptide = unique_peptides,
  hydrophobicity = sapply(X = unique_peptides, 
                          FUN = function(x){Peptides::hydrophobicity(seq = x, scale = "KyteDoolittle")})
)
  
hydro_hgram <- ggplot(data = hydrophobicity_report, aes(x = hydrophobicity)) +
  stat_density(aes(y=after_stat(density)), 
               color="#BCBD22", linewidth = 2, geom="line")+
  theme_bw()+
  ggtitle(label = "Hydrophobicity, KyteDoolittle scale") +
  xlab("Hydrophobicity")+
  ylab("Density")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))


png(filename = "./plots/0_hydrophobicity.png", width = 1000, height = 1000)
hydro_hgram
dev.off()

#### Firsts transformations for the data manipulation
# save the sample_names
sample_names <- meta_data$sample_name
## Transform 0s to NAs
dian_data <- zero_to_NA_label_free(dataset = dian_data, patterns = sample_names)

## Remove contaminants
dian_clean <- remove_contaminants(dataset = dian_data, contaminants = "./raw_data/contaminants.fasta", accession_name = "Accession")

## transform to log2
dian_clean <- log2_to_pattern_label_free(dataset = dian_clean, patterns = sample_names)
openxlsx::write.xlsx(dian_clean, "results/Raw_data_NO_contaminants.xlsx")

# this has become quite useful when then calculating the completeness of the data
to_completeness <- nrow(dian_clean)

#### Data Quality
# Save a palette of colours to work with
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
          "#00B159", "#FCD612", "#FF2F03", "#03D3FF",
          "#2FF923", "#FF03B4", "#f08a0c", "#0d407f",
          "#037A68", "#510840", "#F70D1A", "#e0218a")

### Data to long format
long_format <- dian_clean %>%
  pivot_longer(cols = all_of(sample_names),
               names_to = "sample_name",
               values_to = "intens")

## Add meta_data to long format
long_format <- merge(long_format, meta_data, by = "sample_name")

## Add median calculation (by sample) - Not required
long_format2 <- sample_combination(dataset = long_format, technical_replicates = "sample_name",
                                  expression_unit = "Protein.Group",
                                  sample_names = "sample_name",
                                  intensity_to_combine = "intens",
                                  report_results = T, remove_reps = F)

## Boxplot for raw intensity
intensity_boxplots <- ggplot(long_format, mapping = aes(x = sample_name, y = intens))+
  geom_boxplot(data = long_format, mapping = aes(x = sample_name, y = intens, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Sample name")

png(filename = "./plots/1_raw_intentity.png", width = 1400, height = 1400)
intensity_boxplots
dev.off()


## Completenes
# Calculation
completeness <- long_format %>%
  group_by(sample_name) %>%
  mutate(count_na = 100 - (((sum(is.na(intens))/to_completeness)))*100) %>% 
  subset(select = c(sample_name,count_na)) %>% 
  distinct()
completeness <- merge(completeness, meta_data, by = "sample_name")

# Plot
completeness_barplot <- ggplot(completeness, mapping = aes(y = count_na, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = sample_name))+
  geom_hline(data = completeness ,mapping = aes(yintercept = mean(count_na), col = "red"))+
  theme_bw()+
  xlab("Sample") +
  ylab("% of completeness")+
  ggtitle("Completeness of the data per sample")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Sample name")

png(filename = "./plots/2_completeness.png", width = 1400, height = 1400)
completeness_barplot
dev.off()


## Number of proteins per sample
# Calculation
nprot <- long_format %>%
  group_by(sample_name) %>%
  mutate(count_prot = sum(!is.na(intens))) %>% 
  subset(select = c(sample_name,count_prot)) %>% 
  distinct()
nprot <- merge(nprot, meta_data, by = "sample_name")

# Plot
number_of_prot <- ggplot(nprot, mapping = aes(y = count_prot, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = sample_name))+
  geom_hline(data = nprot ,mapping = aes(yintercept = mean(count_prot), col = "red"))+
  theme_bw()+
  xlab("Sample") +
  ylab("# of proteins")+
  ggtitle("Proteins per sample")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Sample name")

png(filename = "./plots/3_number_of_prots.png", width = 1400, height = 1400)
number_of_prot
dev.off()


## Intensity
# Calculation
int_prot <- long_format %>%
  group_by(sample_name) %>%
  mutate(int_prot = sum(2^intens ,na.rm = T)) %>% 
  subset(select = c(sample_name,int_prot)) %>% 
  distinct()
int_prot <- merge(int_prot, meta_data, by = "sample_name")

# Calculatio of top 10  most intens
int_prot_topn <- long_format %>%
  mutate(intens = 2^intens)%>%
  group_by(sample_name) %>%
  mutate(top_n = rank(-intens)) %>% 
  subset(select = c(sample_name,Protein.Group,intens,top_n)) %>% 
  filter(top_n <= 10) %>% 
  group_by(sample_name) %>% 
  mutate(sum_intens = sum(intens)) %>% 
  subset(select = c(sample_name,sum_intens)) %>% 
  distinct()

# Plot
int_of_prot <- ggplot(int_prot, mapping = aes(y = int_prot, x = sample_name))+
  geom_bar(data = int_prot, stat = "identity", aes(fill = sample_name))+
  geom_hline(data = int_prot ,mapping = aes(yintercept = mean(int_prot), col = "red"))+
  geom_bar(data = int_prot_topn,   
           aes(y = sum_intens, x = sample_name), fill = "black", stat = "identity")+
  
  theme_bw()+
  xlab("Sample") +
  ylab("Intensity")+
  ggtitle("Intensity per sample, top 10 black")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Sample name")

png(filename = "./plots/4_intentity_per_sample.png", width = 1400, height = 1400)
int_of_prot
dev.off()


## Intensity percentage
# Calculation
int_prot_per <- long_format %>%
  group_by(sample_name) %>%
  mutate(int_prot = sum(2^intens ,na.rm = T)) %>% 
  subset(select = c(sample_name,int_prot)) %>% 
  distinct() %>% 
  mutate(int_prot_per = int_prot/int_prot*100)
int_prot_per <- int_prot_per[,c(1,2,3)]

# Calculation of top 10 most intens proteins
int_prot_topn_per <- long_format %>%
  mutate(intens = 2^intens)%>%
  group_by(sample_name) %>%
  mutate(top_n = rank(-intens)) %>% 
  subset(select = c(sample_name,Protein.Group,intens,top_n)) %>% 
  filter(top_n <= 10) %>% 
  group_by(sample_name) %>% 
  mutate(sum_intens = sum(intens)) %>% 
  subset(select = c(sample_name,sum_intens)) %>% 
  distinct()
int_prot_topn_per <- merge(int_prot_per, int_prot_topn_per, 
                           by = "sample_name")
int_prot_topn_per$top_per <- int_prot_topn_per$sum_intens/int_prot_topn_per$int_prot*100

# Plot 
int_of_prot_prot <- ggplot(int_prot_topn_per, mapping = aes(y = int_prot_per, x = sample_name, 
                                                            label = as.character(round(top_per, digits = 2))))+
  geom_bar(data = int_prot_topn_per, stat = "identity", aes(fill = sample_name))+
  geom_hline(data = int_prot_topn_per ,mapping = aes(yintercept = mean(top_per), col = "red"))+
  geom_bar(data = int_prot_topn_per,   
           aes(y = top_per, x = sample_name), fill = "black", stat = "identity")+
  geom_text(data = int_prot_topn_per,
            aes(y = top_per + 5, x = sample_name), size = 8)+
  
  theme_bw()+
  xlab("Sample") +
  ylab("Intensity %")+
  ggtitle("Intensity per sample, top 10 black")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Sample name")

png(filename = "./plots/5_intentity_per_sample_percentage.png", width = 1400, height = 1400)
int_of_prot_prot
dev.off()

## Intensity of Top14 plasma
plasma_proteins <- c("ALBU_HUMAN", "A1AT_HUMAN","APOA1_HUMAN","APOA2_HUMAN","APOB_HUMAN",
                     "A1AG1_HUMAN","HPT_HUMAN","A2MG_HUMAN","TRFE_HUMAN","FIBA_HUMAN", "FIBB_HUMAN",
                     "FIBG_HUMAN","CO3_HUMAN","IGHA1_HUMAN","IGHA2_HUMAN","IGHG1_HUMAN","IGHG2_HUMAN",
                     "IGHG3_HUMAN","IGHG4_HUMAN","IGHM_HUMAN")

int_prot_per_top14 <- long_format %>% 
  filter(Protein.Name %in% plasma_proteins) %>%
  subset(select = c(sample_name,Protein.Name,intens)) %>% 
  mutate(intens = 2^intens) %>% 
  mutate(total_intensity = int_prot_per$int_prot[match(sample_name,int_prot_per$sample_name)]) %>% 
  mutate(total_intensity_per = 100) %>%
  mutate(intens_per = intens/total_intensity*100) %>%
  mutate(intens_per = ifelse(!is.na(intens_per),intens_per,0))
  

int_of_prot_top14 <- ggplot(int_prot_per_top14)+
  geom_bar(data = . %>% 
             subset(select = c(sample_name, total_intensity_per)) %>% 
             distinct(), 
           mapping = aes(x = sample_name, y = total_intensity_per), 
           stat = "identity", fill = "gray", color = "black")+
  geom_bar(data = int_prot_per_top14,
           aes(y = intens_per, x = sample_name, fill = factor(Protein.Name, levels = plasma_proteins)), color = "black", stat = "identity")+
  theme_bw()+
  xlab("Sample") +
  ylab("Intensity %")+
  ggtitle("Intensity per sample, plasma proteins marked")+
  theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Protein")
png(filename = "./plots/6_top14_intensity.png", width = 1000, height = 1000)
int_of_prot_top14
dev.off()

## NA per protein
# Calculation
naprot <- long_format %>%
  group_by(Protein.Group) %>%
  mutate(count_prot = sum(is.na(intens))) %>% 
  subset(select = c(Protein.Group, count_prot)) %>% 
  distinct()

# Plot
na_density <- ggplot(data = naprot, mapping = aes(x = count_prot))+
  geom_histogram(binwidth = 1, bins = 1)+
  theme_bw()+
  ggtitle("NAs density per protein groups")+
  xlab("# NAs")+
  ylab("# proteins")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))

png(filename = "./plots/6_na_density.png", width = 1400, height = 1400)
na_density
dev.off()


## n_prots and total intensity
# Calculations and data modifications
n_prots_total_int <- merge(nprot, int_prot, by = "sample_name")
n_prots_total_int <- n_prots_total_int %>% 
  subset(select = c(sample_name, count_prot, int_prot))

lm_model <- lm(count_prot ~ int_prot, data = n_prots_total_int)
r_squared <- summary(lm_model)$r.squared

# Plot
prots_and_int <- ggplot(data = n_prots_total_int, mapping = aes(x = int_prot, y = count_prot))+
  geom_point(size = 10, colour = "black")+
  geom_smooth(method = "lm", formula = y ~ x)+
  geom_text(x = quantile(n_prots_total_int$int_prot)[4], 
            y = quantile(n_prots_total_int$count_prot)[4], 
            label = paste0("Rsquared: ",as.character(round(r_squared, digits = 5))), color = "red", size= 8)+
  theme_bw()+
  geom_text_repel(nudge_y = 10,
                  size = 6, mapping = aes(label = sample_name))+
  ggtitle("Nº of proteins and Intensity")+
  xlab("Detection intensity")+
  ylab("Nº proteins")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))

png(filename = "./plots/7_number_of_prots_and_intensity.png", width = 1400, height = 1400)
prots_and_int
dev.off()


#### Remove samples
long_format <- remove_samp(dataset = long_format, samples = "A45_ctl_1", reasons = "The sample has a strange chromatogram")
dian_clean <- long_format


#### Impute missing values
dian_clean_imp <- tim(impute = "yes", dataset = dian_clean, NAs_prop = 0.34, 
                      experimental_groups = "cell_line", intensity_to_impute = "intens", 
                      unit_to_impute = "Protein.Group",
                      report_results = T)

# Plot intensity boxplots after imputation
intensity_boxplots_imputation <- ggplot(dian_clean_imp, mapping = aes(x = sample_name, y = imputed_intensity))+
  geom_boxplot(data = dian_clean_imp, mapping = aes(x = sample_name, y = imputed_intensity, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection (Imputed raw intensities)")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Sample name")
png(filename = "./plots/10_intensity_after_tim_application.png", width = 1400, height = 1400)
intensity_boxplots_imputation
dev.off()


#### Data normalization
# median all
median_all <- median(dian_clean_imp$imputed_intensity)
# median by group
dian_clean_imp <- dian_clean_imp %>%
  group_by(sample_name) %>%
  mutate(MED = median(imputed_intensity))
# Do the calculation
dian_clean_imp$normalized_intensity <- (dian_clean_imp$imputed_intensity - dian_clean_imp$MED) + median_all 

# Plot
intensity_boxplots_norm <- ggplot(dian_clean_imp, mapping = aes(x = sample_name, y = normalized_intensity))+
  geom_boxplot(data = dian_clean_imp, mapping = aes(x = sample_name, y = normalized_intensity, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection (Normalized intensities)")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))+
  labs(fill = "Sample name")
png(filename = "./plots/11_normalized_intensity.png", width = 1400, height = 1400)
intensity_boxplots_norm
dev.off()

#### Plots for grouping
### PCA
# matrix extraction
expression_matrix <- reshape2::dcast(dian_clean_imp, Protein.Group ~ sample_name, value.var = "normalized_intensity")
rownames(expression_matrix) <- expression_matrix$Protein.Group 
expression_matrix <- expression_matrix[,-1]

# PCA calculation
pca1 <- prcomp(t(expression_matrix), center = T, scale. = T)
pca1x <- as.data.frame(pca1$x)
pca1x <- pca1x[colnames(expression_matrix),]

pca_meta <- meta_data
rownames(pca_meta) <- pca_meta$sample_name
# data frame to plot
pca1xx <- merge(pca1x, pca_meta, by = "row.names")

# Plot
pca_plot <- ggplot(pca1xx, 
                   aes(x = PC1, y = PC2, label = sample_name))+
  geom_point(size = 8, mapping = aes(color = cell_line))+
  scale_color_manual(values = c("CAV" = "yellow",
                                "A45" = "pink"))+  
  guides(color = guide_legend(title = "Cell line"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 6)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        title = element_text(size = 12),
        legend.text = element_text(size = 10))

png(filename = "./plots/12_first_PCA_same_names_as_first.png", width = 900, height = 900)
pca_plot
dev.off()

### UMAP with clusters
# UMAP calculation
# umap_res <- umap::umap(t(expression_matrix2), n_neighbors = 2, n_components = 2, metric = "euclidean")
# # UMAP data transformations
# umap_df <- as.data.frame(umap_res$layout)
# colnames(umap_df) <- c("UMAP1", "UMAP2")
# umap_df$sample_name <- rownames(umap_df)
# 
# # Plot
# umap_plot <- ggplot(umap_df, 
#                     aes(x = UMAP1, y = UMAP2, label = sample_name))+
#   # geom_point(size = 4, aes(color = Group))+
#   geom_point(size = 4, color = "black")+
#   guides(color = guide_legend(title = "Sample"))+
#   geom_text_repel(nudge_y = 0.01,
#                   size = 4)+
#   xlab("UMAP1") +
#   ylab("UMAP2") +
#   ggtitle(("UMAP"))
# png(filename = "./plots/14_umap_plot.png", width = 1000, height = 1000)
# umap_plot
# dev.off()

#### Extract data
dcast_formula <- as.formula("Protein.Group + Accession + Accession_1 + Gene.names + Gene.names_1 + Protein.Description + Protein.Name ~ sample_name")
expression_matrix_def <- reshape2::dcast(dian_clean_imp, dcast_formula, value.var = "normalized_intensity")
openxlsx::write.xlsx("./results/Expression_matrix.xlsx", x = expression_matrix_def)

