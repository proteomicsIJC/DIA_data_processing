################################
### DIA - Data preprocessing ###
################################

### libraries
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
source("./functions/create_meta_data.R")
source("./functions/first_accession.R")
source("./functions/log2_to_pattern.R")
source("./functions/log2_to_pattern_label_free.R")
source("./functions/make_all_contrasts.R")
source("./functions/presence_vs_no_presence.R")
source("./functions/remove_batch.R")
source("./functions/remove_samp.R")
source("./functions/remove_contaminants.R")
source("./functions/tim.R")
source("./functions/tt_extractor.R")
source("./functions/tt_list_cleaner.R")
source("./functions/upsidedown.R")
source("./functions/xlsx_extractor.R")
source("./functions/zero_to_na.R")
source("./functions/zero_to_na_label_free.R")

### data
raw_data <- readr::read_tsv("./raw_data/report.pg_matrix.tsv")
raw_data <- as.data.frame(raw_data)
colnames(raw_data)[1] <- "Protein.Group"
colnames(raw_data)[2] <- "Accession"
colnames(raw_data)[3] <- "Protein Name"
colnames(raw_data)[4] <- "Gene.names"
colnames(raw_data)[5] <- "Protein Description"

for (i in 6:(ncol(raw_data))){
  colnames(raw_data)[i] <- paste0("s",gsub(x = colnames(raw_data)[i], pattern = ".*\\_",
                                           replacement = ""), collapse = "")
  colnames(raw_data)[i] <- gsub(x = colnames(raw_data)[i], pattern = "\\.raw",
                                replacement = "")
}

dian_data <- raw_data

# expression data
dian <- raw_data[,6:ncol(raw_data)]
rownames(dian) <- raw_data$Protein.Group

# Meta data
meta_data <- as.data.frame(tibble::tibble(
  sample_name = colnames(dian),
  exp_group = rep("Ungrouped", times = 27)
))## THIS IS AN UNGROUPED VERSION <3

# Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

### Work DIAN data
# Transform 0s to NAs
sample_names <- colnames(dian)
dian_data <- zero_to_NA_label_free(dataset = dian_data, patterns = sample_names)

# Remove contaminants
dian_clean <- remove_contaminants(dataset = dian_data, contaminants = "./raw_data/contaminants.fasta", 
                                     accession_name = "Accession")

# remove KTRs
removed_keratines <- dian_clean %>% 
  filter(grepl("^KRT", x = Gene.names))

dian_clean <- dian_clean %>% 
  filter(!grepl("^KRT", x = Gene.names))

# transform to log2
dian_clean <- log2_to_pattern_label_free(dataset = dian_clean, patterns = sample_names)

### Data Quality
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
          "#00B159", "#FCD612", "#FF2F03", "#03D3FF",
          "#2FF923", "#FF03B4", "#f08a0c", "#0d407f",
          "#037A68", "#510840", "#F70D1A", "#e0218a")


# Data to long format
long_format <- dian_clean %>%
  pivot_longer(cols = all_of(sample_names),
               names_to = "sample_name",
               values_to = "intens")
# Add metadat to long format
long_format <- merge(long_format, meta_data, by = "sample_name")

# Add median calculation (by sample)
long_format <- long_format %>%
  group_by(sample_name) %>%
  mutate(MED = median(intens, na.rm=T))

# Boxplot
intensity_boxplots <- ggplot(long_format, mapping = aes(x = sample_name, y = intens))+
  geom_boxplot(data = long_format, mapping = aes(x = sample_name, y = intens, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

png(filename = "./plots/intensity_boxplots.png", width = 1000, height = 1000)
intensity_boxplots
dev.off()


# Completenes
completeness <- long_format %>%
  group_by(sample_name) %>%
  mutate(count_na = 100 - (((sum(is.na(intens))/nrow(dian_clean))))*100) %>% 
  subset(select = c(sample_name,count_na)) %>% 
  distinct()
completeness <- merge(completeness, meta_data, by = "sample_name")

completeness_barplot <- ggplot(completeness, mapping = aes(y = count_na, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = sample_name))+
  geom_hline(data = completeness ,mapping = aes(yintercept = mean(count_na), col = "red"))+
  theme_bw()+
  xlab("Sample") +
  ylab("% of completeness")+
  ggtitle("Completeness of the data per sample")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

png(filename = "./plots/completeness.png", width = 1000, height = 1000)
completeness_barplot
dev.off()

# Number of proteins per sample
nprot <- long_format %>%
  group_by(sample_name) %>%
  mutate(count_prot = sum(!is.na(intens))) %>% 
  subset(select = c(sample_name,count_prot)) %>% 
  distinct()
nprot <- merge(nprot, meta_data, by = "sample_name")

number_of_prot <- ggplot(nprot, mapping = aes(y = count_prot, x = sample_name))+
  geom_bar(stat = "identity", aes(fill = sample_name))+
  geom_hline(data = nprot ,mapping = aes(yintercept = mean(count_prot), col = "red"))+
  theme_bw()+
  xlab("Sample") +
  ylab("# of proteins")+
  ggtitle("Proteins per sample")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

png(filename = "./plots/number_of_prot.png", width = 1000, height = 1000)
number_of_prot
dev.off()

# NA per protein
naprot <- long_format %>%
  group_by(Protein.Group) %>%
  mutate(count_prot = sum(is.na(intens))) %>% 
  subset(select = c(Protein.Group, count_prot)) %>% 
  distinct()

# NAs density per protein
na_density <- ggplot(data = naprot, mapping = aes(x = count_prot))+
  geom_histogram(binwidth = 1, bins = 1)+
  theme_bw()+
  ggtitle("NAs density per protein groups")+
  xlab("# NAs")+
  ylab("# proteins")+
  theme(legend.position= "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

png(filename = "./plots/na_density.png", width = 1000, height = 1000)
na_density
dev.off()

### Remove samples
long_format <- remove_samp(dataset = long_format)
dian_clean <- long_format

### Impute missing values
colnames(dian_clean)[3] <- "protein_group"
dian_clean$normalized_intensity <- dian_clean$intens
dian_clean_imp <- tim(impute = "no", dataset = dian_clean)
dian_clean_imp <- dian_clean_imp %>% 
  subset(select = -c(normalized_intensity))

### Normalization by median
# median all
median_all <- median(dian_clean_imp$intens)

# MED
dian_clean_imp <- dian_clean_imp %>%
  group_by(sample_name) %>%
  mutate(MED = median(intens))

# Do the calculation
dian_clean_imp$normalized_intensity <- (dian_clean_imp$intens - dian_clean_imp$MED) + median_all 

# Plot the result
intensity_boxplots_norm <- ggplot(dian_clean_imp, mapping = aes(x = sample_name, y = normalized_intensity))+
  geom_boxplot(data = dian_clean_imp, mapping = aes(x = sample_name, y = normalized_intensity, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection (Normalized intensities)")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
png(filename = "./plots/normalized_intensity.png", width = 1000, height = 1000)
intensity_boxplots_norm
dev.off()

### Extract data
# annotation
colnames(dian_clean_imp)[2] <- "Accession"
colnames(dian_clean_imp)[3] <- "Protein Group"
annotation <- dian_clean_imp %>% 
  subset(select = c(`Protein Group`, Accession, Gene.names, `Protein Description`, `Protein Name`)) %>% 
  distinct()
readr::write_tsv("./results/annotation.tsv", x = annotation)

# meta just in case
readr::write_tsv("./results/meta_data_used_no_clustered.tsv", x = meta_data)

# matrix
expression_matrix <- reshape2::dcast(dian_clean_imp, `Protein Group` ~ sample_name, value.var="normalized_intensity")
expression_matrix_def <- merge(annotation, expression_matrix, by.x = "Protein Group", by.y = "Protein Group")

readr::write_tsv("./results/expression_matrix.tsv", x = expression_matrix_def)


### PCAs
# PCA computation
rownames(expression_matrix) <- expression_matrix$`Protein Group`
expression_matrix <- expression_matrix[,-1]
pca1 <- prcomp(t(expression_matrix), center = TRUE, scale. = TRUE)
pca1x <- as.data.frame(pca1$x)
pca1x <- pca1x[colnames(expression_matrix),]

# meta data for the PCA
meta_pez <- openxlsx::read.xlsx("./raw_data/experimental_design.xlsx")
rownames(meta_pez) <- meta_pez$Sample_name_Ignasi_name
meta_pez <- meta_pez[colnames(expression_matrix),]

pca1xx <- merge(pca1x, meta_pez, by = "row.names")

# Represent the PCA
pca_plot <- ggplot(pca1xx, 
                 aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Tissue, color = Dose),size = 4) +
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 1), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 1), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))

pca_plot

png(filename = "./plots/pca_tissues_dose.png", width = 1000, height = 1000)
pca_plot
dev.off()

######## A BETTER WAY TO DO THE PCAs, check this works well for you !! 
### PCA plotting
# matrix extraction 
expression_matrix <- reshape2::dcast(dian_clean_imp, protein_group ~ sample_name, value.var = "normalized_intensity")
rownames(expression_matrix) <- expression_matrix$protein_group 
expression_matrix <- expression_matrix[,-1]
pca1 <- prcomp(t(expression_matrix), center = T, scale. = T)
pca1x <- as.data.frame(pca1$x)
pca1x <- pca1x[colnames(expression_matrix),]

# meta for pca
pca_meta <- as.data.frame(tibble(
  sample_rep = rownames(pca1x)
))

pca_meta$sample <- pca_meta$sample_rep
for (i in 1:length(pca_meta$sample)){
  pca_meta$sample[i] <- unlist(strsplit(pca_meta$sample_rep[i], "_")[[1]])[1]
}
rownames(pca_meta) <- pca_meta$sample_rep

# data frame to plot
pca1xx <- merge(pca1x, pca_meta, by = "row.names")

# Represent the PCA
pca_plot <- ggplot(pca1xx, 
                   aes(x = PC1, y = PC2, label = sample_rep))+
  geom_point(size = 4, aes(color = sample))+
  guides(color = guide_legend(title = "Sample"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 6)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))+
  geom_polygon(aes(group = sample, fill = sample), alpha = 0.2, show.legend = FALSE)
pca_plot

png(filename = "./plots/6_pca_replicates.png", width = 1400, height = 1400)
pca_plot
dev.off()






