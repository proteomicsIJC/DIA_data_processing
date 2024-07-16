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
library(sva)

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
source("./functions/log2_to_pattern.R")
source("./functions/log2_to_pattern_label_free.R")
source("./functions/make_all_contrasts.R")
source("./functions/remove_batch.R")
source("./functions/remove_samp.R")
source("./functions/remove_contaminants.R")
source("./functions/tim.R")
source("./functions/upsidedown.R")
source("./functions/xlsx_tt.R")
source("./functions/zero_to_na.R")
source("./functions/zero_to_na_label_free.R")
source("./functions/sample_combination.R")

### DIAN data
raw_data <- readr::read_tsv("./raw_data/report.pg_matrix.tsv", locale = readr::locale(encoding = "latin1"))
raw_data <- as.data.frame(raw_data)
colnames(raw_data)[1] <- "Protein.Group"
colnames(raw_data)[2] <- "Accession"
colnames(raw_data)[3] <- "Protein Name"
colnames(raw_data)[4] <- "Gene.names"
colnames(raw_data)[5] <- "Protein Description"

### mutate colnames representing the samples
## This will change at every DIAN run !!
colnames(raw_data)[6:ncol(raw_data)] <- gsub("E.*REC_", "REC_", colnames(raw_data)[6:ncol(raw_data)])
colnames(raw_data)[6:ncol(raw_data)] <- gsub("_DIA", "", colnames(raw_data)[6:ncol(raw_data)])
colnames(raw_data)[6:ncol(raw_data)] <- gsub("\\.raw", "", colnames(raw_data)[6:ncol(raw_data)])
raw_data <- raw_data[, -grep("Phospho", colnames(raw_data))]

# rec_cho_2_16_11 <- apply(raw_data[, 6:7], 1, median, na.rm = TRUE)
# rec_cho_9_11 <- apply(raw_data[, 8:9], 1, median, na.rm = TRUE)
# rec_cho_16_11 <- apply(raw_data[, 10:11], 1, median, na.rm = TRUE)
# rec_ut_9_11 <- apply(raw_data[, 14:15], 1, median, na.rm = TRUE)
# rec_ut_16_11 <- apply(raw_data[, 16:17], 1, median, na.rm = TRUE)
# rec_ut_2_16_11 <- apply(raw_data[, 12:13], 1, median, na.rm = TRUE)
# 
# raw_data <- cbind.data.frame(raw_data[, 1:5], 
#                               "REC_CHO_9_11" = rec_cho_9_11,
#                               "REC_CHO_16_11" = rec_cho_16_11,
#                               "REC_CHO_2_16_11" = rec_cho_2_16_11,
#                               "REC_UT_9_11" = rec_ut_9_11,
#                               "REC_UT_16_11" = rec_ut_16_11,
#                               "REC_UT_2_16_11" = rec_ut_2_16_11)

### Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

#### Annotation completion with Uniprot.ws
### Uniprot.ws work
## Accession_1
raw_data <- raw_data %>% 
  group_by(Accession) %>% 
  mutate(Accession_1 = unlist(strsplit(x = Accession, split = "\\;"), use.names = F)[1]) %>%
  mutate(Gene.names_1 = unlist(strsplit(x = Gene.names, split = "\\;"), use.names = F)[1]) %>% 
  relocate(Accession_1, .after = Accession) %>% 
  relocate(Gene.names_1, .after = Gene.names) %>% 
  ungroup()

## Uniprot.ws petition
id_mapping <- UniProt.ws::mapUniProt("UniProtKB_AC-ID", "UniProtKB", query = raw_data$Accession_1)

raw_data <- raw_data %>% 
  group_by(Accession) %>% 
  mutate(`Protein Description` = ifelse(Accession_1 %in% id_mapping$From,
                                        id_mapping$Protein.names[match(Accession_1, id_mapping$From)],
                                        `Protein Description`))
raw_data <- as.data.frame(raw_data)
dian_data <- raw_data

### Expression data
dian <- raw_data[,8:ncol(raw_data)]
rownames(dian) <- raw_data$Protein.Group

# Meta data
# This is an ungrouped version of the data that helps the imputation
# meta_data <- as.data.frame(tibble::tibble(
#   sample_name = colnames(dian),
#   exp_group = rep("Ungrouped", times = ncol(dian))
# ))

#### Firsts transformations for the data manipulation
# save the sample_names
sample_names <- colnames(dian)
## Transform 0s to NAs
dian_data <- zero_to_NA_label_free(dataset = dian_data, patterns = sample_names)

## Remove contaminants
dian_clean <- remove_contaminants(dataset = dian_data, contaminants = "./raw_data/contaminants.fasta", 
                                  accession_name = "Accession")

## transform to log2
dian_clean <- log2_to_pattern_label_free(dataset = dian_clean, patterns = sample_names)
openxlsx::write.xlsx("./results/protein_list.xlsx", x = dian_clean)

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
# n_prots_total_int <- n_prots_total_int[,c(1,2,4)]

lm_model <- lm(count_prot ~ int_prot, data = n_prots_total_int)
r_squared <- summary(lm_model)$r.squared
# Plot
prots_and_int <- ggplot(data = n_prots_total_int, mapping = aes(x = int_prot, y = count_prot))+
  geom_point(size = 10, colour = "black")+
  geom_smooth(method = "lm", formula = y ~ x)+
  geom_text(x = 1.3e9, y = 550, label = paste0("Rsquared: ",as.character(round(r_squared, digits = 5))), color = "red", size= 8)+
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

## initial sample quantity and number of proteins
# Add the info
# sample_extra_info <- as.data.frame(readr::read_tsv("./raw_data/sample_extra_info.tsv",
#                                                    locale = readr::locale(encoding = "latin1")))
# sample_extra_info <- sample_extra_info[,c(2,3,4)]
# colnames(sample_extra_info) <- c("sample_name","volume","injection_volume")
# n_prots_total_int <- merge(n_prots_total_int, sample_extra_info, by = "sample_name")
# # Calculations
# lm_model <- lm(count_prot ~ volume, data = n_prots_total_int)
# r_squared <- summary(lm_model)$r.squared
# # Plots
# prots_and_vol <- ggplot(data = n_prots_total_int, mapping = aes(x = volume, y = count_prot), label = sample_name)+
#   geom_point(size = 10, colour = "black")+
#   geom_smooth(method = "lm", formula = y ~ x)+
#   geom_text(x = 100, y = 300, label = paste0("Rsquared: ",as.character(round(r_squared, digits = 5))), color = "red", size= 8)+
#   theme_bw()+
#   geom_text_repel(nudge_y = 10,
#                   size = 6, mapping = aes(label = sample_name))+
#   ggtitle("Nº of proteins and Initial sample volume")+
#   xlab("Initial sample volume")+
#   ylab("Nº proteins")+
#   theme(legend.position= "none",
#         axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
#         axis.text.y = element_text(size = 15, face = "bold"),
#         axis.title = element_text(size = 20),
#         title = element_text(size = 22),
#         legend.text = element_text(size = 14, face = "bold"))
# 
# png(filename = "./plots/8_number_of_prots_and_volume.png", width = 1400, height = 1400)
# prots_and_vol
# dev.off()
# 
# ## injected sample and number of proteins
# # Calculations
# lm_model <- lm(count_prot ~ injection_volume, data = n_prots_total_int)
# r_squared <- summary(lm_model)$r.squared
# # Plot
# inj_and_vol <- ggplot(data = n_prots_total_int, mapping = aes(x = injection_volume, y = count_prot), label = sample_name)+
#   geom_point(size = 10, colour = "black")+
#   geom_smooth(method = "lm", formula = y ~ x)+
#   geom_text(x = 0.025, y = 600, label = paste0("Rsquared: ",as.character(round(r_squared, digits = 5))), color = "red", size= 8)+
#   theme_bw()+
#   geom_text_repel(nudge_y = 15,
#                   size = 6, mapping = aes(label = sample_name))+
#   ggtitle("Nº of proteins and Injected volume")+
#   xlab("Injected volume")+
#   ylab("Nº proteins")+
#   theme(legend.position= "none",
#         axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
#         axis.text.y = element_text(size = 15, face = "bold"),
#         axis.title = element_text(size = 20),
#         title = element_text(size = 22),
#         legend.text = element_text(size = 14, face = "bold"))
# 
# png(filename = "./plots/9_number_of_prots_and_injection.png", width = 1400, height = 1400)
# inj_and_vol
# dev.off()

### Remove samples
long_format <- remove_samp(dataset = long_format)
dian_clean <- long_format

### Sample Deconvolution
## Preare the data for the sample deconvolution
sample_deconv <- readr::read_tsv("./raw_data/meta_data_for_deconvolution.tsv")
dian_clean <- merge(dian_clean, sample_deconv, by = "sample_name")
## Do the deconvolution
dian_clean <- sample_combination(dataset = dian_clean, 
                                 technical_replicates = "sample_name", 
                                 sample_names = "real_sample",
                                 expression_unit = "Protein.Group", 
                                 intensity_to_combine = "intens", 
                                 remove_reps = T,report_results = T)
## Changing of colnames after the deconvolution
colnames(dian_clean)[grep("real_sample", colnames(dian_clean))] <- "sample_name" 
colnames(dian_clean)[grep("combined_intensity", colnames(dian_clean))] <- "intens"

#### Impute missing values
## meta_data for the imputation
experimental_groups <- c("REC_CHO_09_11" = "CHO", 
                         "REC_CHO_16_11" = "CHO",
                         "REC_CHO_2_16_11" = "CHO", 
                         "REC_UT_09_11" = "UT",
                         "REC_UT_16_11" = "UT", 
                         "REC_UT_2_16_11" = "UT")

## do the imputation
dian_clean_imp <- tim(impute = "yes", dataset = dian_clean, unit_to_impute = "Protein.Group",
                      NAs_prop = 0.3, intensity_to_impute = "intens", 
                      experimental_groups = experimental_groups,
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
## median all
median_all <- median(dian_clean_imp$imputed_intensity)

## median by group
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
rownames(expression_matrix) <- expression_matrix$protein_group 
expression_matrix <- expression_matrix[,-1]
# PCA calculation
pca1 <- prcomp(t(expression_matrix), center = T, scale. = T)
pca1x <- as.data.frame(pca1$x)
pca1x <- pca1x[colnames(expression_matrix),]
pca1xx <- pca1x
pca1xx$sample_name <- rownames(pca1xx)

# Plot
pca_plot <- ggplot(pca1xx, 
                   aes(x = PC1, y = PC2, label = sample_name))+
  geom_point(size = 4, color = "black")+
  guides(color = guide_legend(title = "Sample"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 4)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))
#### ADD THE SAME AS THE ELOI GARI !

png(filename = "./plots/12_first_PCA.png", width = 1000, height = 1000)
pca_plot
dev.off()

#### Remove batch effect
### Check the ordering of the batches
columns_checker(dataset = dian_clean_imp, 
                unit_of_analysis = "Protein.Group", intensity = "normalized_intensity", 
                tmt = F)

### Remove the batch effect
dian_clean_imp_unbatch <- remove_batch(dataset = dian_clean_imp,
                                       use_combat = T,
                                       remove = "yes", 
                                       intensity = "normalized_intensity", unit_of_work = "Protein.Group", 
                                       tmt = F, where_is_the_batch1 = c("batch1","batch2",
                                                                        "batch2","batch2",
                                                                        "batch1","batch2"))


## Extract the expression_matrix
expression_matrix2 <- reshape2::dcast(dian_clean_imp_unbatch, Protein.Group ~ sample_name, value.var = "unbatched_intensity")
rownames(expression_matrix2) <- expression_matrix2$Protein.Group 
expression_matrix2 <- expression_matrix2[,-1]

# PCA calculation
pca1 <- prcomp(t(expression_matrix2), center = T, scale. = T)
pca1x <- as.data.frame(pca1$x)
pca1x <- pca1x[colnames(expression_matrix2),]
pca1xx <- pca1x
pca1xx$sample_name <- rownames(pca1xx)

# Plot
pca_plot <- ggplot(pca1xx, 
                   aes(x = PC1, y = PC2, label = sample_name))+
  # geom_point(size = 4, aes(color = Group))+
  geom_point(size = 4, color = "black")+
  guides(color = guide_legend(title = "Sample"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 4)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))

png(filename = "./plots/13_PCA_with_groups.png", width = 1000, height = 1000)
pca_plot
dev.off()

### UMAP with clusters
# UMAP calculation
umap_res <- umap::umap(t(expression_matrix2), n_neighbors = 2, n_components = 2, metric = "euclidean")
# UMAP data transformations
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sample_name <- rownames(umap_df)

# Plot
umap_plot <- ggplot(umap_df, 
                    aes(x = UMAP1, y = UMAP2, label = sample_name))+
  # geom_point(size = 4, aes(color = Group))+
  geom_point(size = 4, color = "black")+
  guides(color = guide_legend(title = "Sample"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 4)+
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle(("UMAP"))
png(filename = "./plots/14_umap_plot.png", width = 1000, height = 1000)
umap_plot
dev.off()

#### Extract data
## Annotation
annotation <- dian_clean
annotation <- annotation %>% 
  subset(select = c(Protein.Group, Accession, Accession_1, Gene.names, Gene.names_1,`Protein Description`, `Protein Name`)) %>% 
  distinct()
readr::write_tsv("./results/annotation.tsv", x = annotation)

## expression_matrix (protein_list)
expression_matrix_def <- reshape2::dcast(dian_clean_imp_unbatch, Protein.Group ~ sample_name, value.var="unbatched_intensity" )
expression_matrix_def <- merge(annotation, expression_matrix_def, by = "Protein.Group")
readr::write_tsv("./results/expression_matrix_def.tsv", x = expression_matrix_def)
