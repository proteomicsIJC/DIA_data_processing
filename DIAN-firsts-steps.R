########################
### Data preparation ###
########################

### libraries
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#source("log2Transf.R")
#source("proteinGroupsCleanner.R")
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
source("./functions/upsidedown.R")
source("./functions/xlsx_tt.R")
source("./functions/zero_to_na.R")
source("./functions/zero_to_na_label_free.R")
library(readr) 

### data
raw_data <- readr::read_tsv("./raw_data/report.pg_matrix.tsv", locale = readr::locale(encoding = "latin1"))
raw_data <- as.data.frame(raw_data)
colnames(raw_data)[1] <- "Protein.Group"
colnames(raw_data)[2] <- "Accession"
colnames(raw_data)[3] <- "Protein Name"
colnames(raw_data)[4] <- "Gene.names"
colnames(raw_data)[5] <- "Protein Description"
for (i in 6:(ncol(raw_data))){
  ## gsub for the general pattern
  colnames(raw_data)[i] <- gsub(pattern = "R:\\\\R_PROTEOMICA\\\\1_PROJECTS_OE480\\\\P094_Francesc Xavier Avilés", replacement = "", x = colnames(raw_data)[i])
  colnames(raw_data)[i] <- gsub(pattern = "\\\\20240426", replacement = "", x = colnames(raw_data)[i])
  ## gsub remove .raw
  colnames(raw_data)[i] <- gsub(pattern = "\\.raw",replacement = "", x = colnames(raw_data)[i])
  ## gsub remove "coletillas"
  colnames(raw_data)[i] <- gsub(pattern = "_B", replacement = "", x = colnames(raw_data)[i])
  colnames(raw_data)[i] <- gsub(pattern = "_15ng", replacement = "", x = colnames(raw_data)[i])
  ## Change PR for user name
  colnames(raw_data)[i] <- gsub(pattern = "_FAVI019", replacement = "HNSCC", x = colnames(raw_data)[i])
}

### Uniprot.ws work
## Accession_1
raw_data <- raw_data %>% 
  group_by(Accession) %>% 
  mutate(Accession_1 = unlist(strsplit(x = Accession, split = "\\;"), use.names = F)[1]) %>% 
  relocate(Accession_1, .after = Accession) %>% 
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

# expression data
dian <- raw_data[,7:ncol(raw_data)]
rownames(dian) <- raw_data$Protein.Group

## barplot(apply(dian,2, sum, na.rm = T), las = 2)

# Meta data
meta_data <- as.data.frame(tibble::tibble(
  sample_name = colnames(dian),
  exp_group = rep("Ungrouped", times = 19)
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

# transform to log2
dian_clean <- log2_to_pattern_label_free(dataset = dian_clean, patterns = sample_names)
to_completeness <- nrow(dian_clean)

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
##long_format <- long_format %>%
##  group_by(sample_name) %>%
##  mutate(MED = median(intens, na.rm=T))

# Boxplot
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

# Completenes
completeness <- long_format %>%
  group_by(sample_name) %>%
  mutate(count_na = 100 - (((sum(is.na(intens))/to_completeness)))*100) %>% 
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
int_prot <- long_format %>%
  group_by(sample_name) %>%
  mutate(int_prot = sum(2^intens ,na.rm = T)) %>% 
  subset(select = c(sample_name,int_prot)) %>% 
  distinct()
int_prot <- merge(int_prot, meta_data, by = "sample_name")

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

## Intensity per 
int_prot_per <- long_format %>%
  group_by(sample_name) %>%
  mutate(int_prot = sum(2^intens ,na.rm = T)) %>% 
  subset(select = c(sample_name,int_prot)) %>% 
  distinct() %>% 
  mutate(int_prot_per = int_prot/int_prot*100)
int_prot_per <- int_prot_per[,c(1,2,3)]

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
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))

png(filename = "./plots/6_na_density.png", width = 1400, height = 1400)
na_density
dev.off()

# n_prots and total intensity
n_prots_total_int <- merge(nprot, int_prot, by = "sample_name")
n_prots_total_int <- n_prots_total_int[,c(1,2,4)]

lm_model <- lm(count_prot ~ int_prot, data = n_prots_total_int)
r_squared <- summary(lm_model)$r.squared

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

# initial sample quantity and number of proteins
sample_extra_info <- as.data.frame(readr::read_tsv("./raw_data/sample_extra_info.tsv",
                                                   locale = readr::locale(encoding = "latin1")))


sample_extra_info <- sample_extra_info[,c(2,3,4)]
colnames(sample_extra_info) <- c("sample_name","volume","injection_volume")
n_prots_total_int <- merge(n_prots_total_int, sample_extra_info, by = "sample_name")

lm_model <- lm(count_prot ~ volume, data = n_prots_total_int)
r_squared <- summary(lm_model)$r.squared

prots_and_vol <- ggplot(data = n_prots_total_int, mapping = aes(x = volume, y = count_prot), label = sample_name)+
  geom_point(size = 10, colour = "black")+
  geom_smooth(method = "lm", formula = y ~ x)+
  geom_text(x = 100, y = 300, label = paste0("Rsquared: ",as.character(round(r_squared, digits = 5))), color = "red", size= 8)+
  theme_bw()+
  geom_text_repel(nudge_y = 10,
                  size = 6, mapping = aes(label = sample_name))+
  ggtitle("Nº of proteins and Initial sample volume")+
  xlab("Initial sample volume")+
  ylab("Nº proteins")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))

png(filename = "./plots/8_number_of_prots_and_volume.png", width = 1400, height = 1400)
prots_and_vol
dev.off()

# injected sample and number of proteins
# n_prots_total_int <- n_prots_total_int %>%
#   filter(!sample_name %in% c("HNSCC_008","HNSCC_009"))

lm_model <- lm(count_prot ~ injection_volume, data = n_prots_total_int)
r_squared <- summary(lm_model)$r.squared

inj_and_vol <- ggplot(data = n_prots_total_int, mapping = aes(x = injection_volume, y = count_prot), label = sample_name)+
  geom_point(size = 10, colour = "black")+
  geom_smooth(method = "lm", formula = y ~ x)+
  geom_text(x = 0.025, y = 600, label = paste0("Rsquared: ",as.character(round(r_squared, digits = 5))), color = "red", size= 8)+
  theme_bw()+
  geom_text_repel(nudge_y = 15,
                  size = 6, mapping = aes(label = sample_name))+
  ggtitle("Nº of proteins and Injected volume")+
  xlab("Injected volume")+
  ylab("Nº proteins")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))

png(filename = "./plots/9_number_of_prots_and_injection.png", width = 1400, height = 1400)
inj_and_vol
dev.off()

### Remove samples
long_format <- remove_samp(dataset = long_format, samples = c(""))
dian_clean <- long_format

### Sample Deconvolution
############################ SOURCE AND MAKE THE FUNCTION WORK HERE !!

### Impute missing values
colnames(dian_clean)[3] <- "protein_group"
dian_clean_imp <- tim(impute = "yes", dataset = dian_clean, NAs_prop = 0.1, intensity_to_impute = "intens")

intensity_boxplots_norm <- ggplot(dian_clean_imp, mapping = aes(x = sample_name, y = imputed_intensity))+
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
intensity_boxplots_norm
dev.off()

### RENORMALIZE !
###dian_clean_imp <- dian_clean_imp %>% 
###  subset(select = -c(normalized_intensity))

# median all
median_all <- median(dian_clean_imp$imputed_intensity)

# MED
dian_clean_imp <- dian_clean_imp %>%
  group_by(sample_name) %>%
  mutate(MED = median(imputed_intensity))

# Do the calculation
dian_clean_imp$normalized_intensity <- (dian_clean_imp$imputed_intensity - dian_clean_imp$MED) + median_all 

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
pca_meta <- pca_meta %>%
  mutate(batch = ifelse((str_detect(pattern = "post", string = sample_rep)) | 
                          (str_detect(pattern = "bC", string = sample_rep)),"batch_2","batch_1"))

rownames(pca_meta) <- pca_meta$sample_rep


# data frame to plot
pca1xx <- merge(pca1x, pca_meta, by = "row.names")

# Represent the PCA
pca_plot <- ggplot(pca1xx, 
                   aes(x = PC1, y = PC2, label = sample_rep))+
  geom_point(size = 4, color = "black")+
  guides(color = guide_legend(title = "Sample"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 4)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))
####  geom_polygon(aes(group = sample, fill = sample), alpha = 0.2, show.legend = FALSE)

png(filename = "./plots/12_first_PCA.png", width = 1000, height = 1000)
pca_plot
dev.off()

### Kmeans
number_of_clusts <- factoextra::fviz_nbclust(x = t(expression_matrix), kmeans, method = "silhouette", k.max = 17)

png(filename = "./plots/13_number_of_clusters.png", width = 1000, height = 1000)
number_of_clusts
dev.off()

kmeans_results <- kmeans(x = t(expression_matrix), 2, nstart = 25)
clusters <- as.data.frame(tibble::tibble(
  sample = names(kmeans_results$cluster),
  group = paste0("group_",kmeans_results$cluster)
))

## PCA with clusters
pca1xx <- pca1xx %>% 
  mutate(Group = clusters$group[match(sample_rep, clusters$sample)])

# Represent the PCA
pca_plot <- ggplot(pca1xx, 
                   aes(x = PC1, y = PC2, label = sample_rep))+
  geom_point(size = 4, aes(color = Group))+
  guides(color = guide_legend(title = "Sample"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 4)+
  xlab(paste("PC1(", round(100*(pca1$sdev[1]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ylab(paste("PC2(", round(100*(pca1$sdev[2]^2/sum(pca1$sdev^2)), 2), "%)", sep = "")) +
  ggtitle(("Principal Component Analysis (PCA)"))

png(filename = "./plots/13_PCA_with_groups.png", width = 1000, height = 1000)
pca_plot
dev.off()

## UMAP with clusters
umap_res <- umap::umap(t(expression_matrix), n_neighbors = 9, n_components = 2, metric = "euclidean")
umap_df <- umap_res$layout
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df <- merge(umap_df, clusters, by.x = "row.names", by.y = "sample")
colnames(umap_df)[1] <- "sample"
colnames(umap_df)[4] <- "Group"

umap_plot <- ggplot(umap_df, 
                    aes(x = UMAP1, y = UMAP2, label = sample))+
  geom_point(size = 4, aes(color = Group))+
  guides(color = guide_legend(title = "Sample"))+
  geom_text_repel(nudge_y = 0.01,
                  size = 4)+
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle(("UMAP"))

png(filename = "./plots/14_umap_plot.png", width = 1000, height = 1000)
umap_plot
dev.off()

### Extract data
# dian clean column name change
colnames(dian_clean_imp)[2] <- "Accession"
colnames(dian_clean_imp)[3] <- "Protein Group"
annotation <- dian_clean_imp %>% 
  subset(select = c(`Protein Group`, Accession, Gene.names, `Protein Description`, `Protein Name`)) %>% 
  distinct()

# annotation
annotation <- dian_clean
colnames(annotation)[2] <- "Accession"
colnames(annotation)[3] <- "Protein Group"
annotation <- annotation %>% 
  subset(select = c(`Protein Group`, Accession, Gene.names, `Protein Description`, `Protein Name`)) %>% 
  distinct()
readr::write_tsv("./results/annotation.tsv", x = annotation)

# meta just in case
readr::write_tsv("./results/meta_data_unclustered.tsv", x = meta_data)

# matrix
expression_matrix <- reshape2::dcast(dian_clean_imp, `Protein Group` ~ sample_name, value.var="normalized_intensity" )
expression_matrix_def <- merge(annotation, expression_matrix, by.x = "Accession", by.y = "Protein Group")
readr::write_tsv("./results/expression_matrix_def.tsv", x = expression_matrix_def)
