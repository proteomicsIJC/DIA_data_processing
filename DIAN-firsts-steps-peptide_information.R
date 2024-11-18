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
### Peptide information
peptide_info <- readr::read_tsv("./raw_data/report.pr_matrix.tsv")

### Hela information
hela_peptides1 <- readr::read_tsv("./raw_data/Hela_raw_files/report.pr_matrix_1.tsv")
colnames(hela_peptides1)[ncol(hela_peptides1)] <- "intens"
hela_peptides1$sample_name <- "Hela_1"

hela_peptides2 <- readr::read_tsv("./raw_data/Hela_raw_files/report.pr_matrix_2.tsv")
colnames(hela_peptides2)[ncol(hela_peptides2)] <- "intens"
hela_peptides2$sample_name <- "Hela_2"

hela_peptides3 <- readr::read_tsv("./raw_data/Hela_raw_files/report.pr_matrix_3_1.tsv")
colnames(hela_peptides3)[ncol(hela_peptides3)] <- "intens"
hela_peptides3$sample_name <- "Hela_3"

hela_peptides4 <- readr::read_tsv("./raw_data/Hela_raw_files/report.pr_matrix_3_2.tsv")
colnames(hela_peptides4)[ncol(hela_peptides4)] <- "intens"
hela_peptides4$sample_name <- "Hela_4"

hela_peptides5 <- readr::read_tsv("./raw_data/Hela_raw_files/report.pr_matrix_3_3.tsv")
colnames(hela_peptides5)[ncol(hela_peptides5)] <- "intens"
hela_peptides5$sample_name <- "Hela_5"

hela_peptides6 <- readr::read_tsv("./raw_data/Hela_raw_files/report.pr_matrix_3_4.tsv")
colnames(hela_peptides6)[ncol(hela_peptides6)] <- "intens"
hela_peptides6$sample_name <- "Hela_6"


hela_peptides <- rbind(hela_peptides1, hela_peptides2, hela_peptides3, hela_peptides4, hela_peptides5, hela_peptides6)
remove(hela_peptides1);remove(hela_peptides2) 
remove(hela_peptides3);remove(hela_peptides4) 
remove(hela_peptides5);remove(hela_peptides6)

### MetaData
meta_data <- as.data.frame(readr::read_tsv("./raw_data/meta_data.tsv", locale = readr::locale(encoding = "latin1")))
# save the sample_names
sample_names <- meta_data$sample_name

#### Raw Data cleaning 
colnames(peptide_info) <- sapply(X = colnames(peptide_info), 
                             FUN = function(x) {
                               if (x %in% meta_data$file_name){
                                 meta_data$sample_name[match(x, meta_data$file_name)]}
                               else {
                                 x
                               }
                             })

#### Peptide based analysis
peptide_info <- peptide_info %>% 
  filter(!duplicated(Stripped.Sequence))
hela_peptides <- hela_peptides %>% 
  filter(!duplicated(Stripped.Sequence)) 

peptide_info_long <- peptide_info %>% 
  pivot_longer(cols = all_of(sample_names),
               names_to = "sample_name",
               values_to = "intens") %>% 
  subset(select = c(Stripped.Sequence,sample_name, intens)) %>% 
  filter(!is.na(intens)) %>% 
  subset(select = -c(intens))
peptide_info_long$type_of_sample <- "sample"

hela_peptides_long <- hela_peptides %>% 
  subset(select = c(Stripped.Sequence,sample_name, intens)) %>% 
  filter(!is.na(intens)) %>% 
  subset(select = -c(intens))
hela_peptides_long$type_of_sample <- "Hela"
peptide_info_long <- rbind(peptide_info_long, hela_peptides_long)

## Missed cleavages
cleavage_report <- tibble::tibble(
  # The peptides
  peptide = peptide_info_long$Stripped.Sequence,
  # The samples
  sample_name = peptide_info_long$sample_name,
  
  # Count Ks and Rs
  k_or_r = sapply(X = peptide_info_long$Stripped.Sequence, FUN = function(x) {
    length(grep(x = unlist(strsplit(x, split = ""),use.names = F), pattern = "K|R"))  
  }),
  
  # Count Ks and Rs that happen at the end
  final_k_or_r = sapply(X = peptide_info_long$Stripped.Sequence, FUN = function(x) {
    final_letter <- unlist(strsplit(x, split = ""), use.names = F)
    final_letter <- final_letter[length(final_letter)]
    length(grep(x = final_letter, pattern = "K|R"))
  }),
  # Count KPs and RPs
  kp_or_rp = sapply(X = peptide_info_long$Stripped.Sequence, FUN = function(x) {
    stringr::str_count(string = x, pattern = "RP|KP")
  })
)

cleavage_report <- cleavage_report %>% 
  mutate(missed_cleavages = k_or_r - final_k_or_r - kp_or_rp)

cleavage_report_per <- cleavage_report %>% 
  mutate(missed_cleavages = as.character(missed_cleavages)) %>% 
  group_by(sample_name) %>% 
  mutate(number_of_peptides = n()) %>% 
  ungroup() %>% 
  group_by(sample_name, missed_cleavages) %>% 
  mutate(number_of_cleavages = n()) %>% 
  ungroup() %>% 
  subset(select = c(sample_name,missed_cleavages,number_of_peptides,number_of_cleavages)) %>% 
  distinct() %>% 
  mutate(cleavage_per = number_of_cleavages/number_of_peptides*100)

cleavage_report_per$sample_type <- sapply(X = cleavage_report_per$sample_name, 
                                          FUN = function(x) {
                                            if (grepl(x = x, pattern = "Hela")){
                                              "Hela"
                                            } else {
                                              "Sample"
                                            }
                                          })

# Plot
missed_cleavages <- ggplot(data = cleavage_report_per)+
  # geom_bar(stat = "identity", aes(x = missed_cleavages, y = cleavage_per, fill = missed_cleavages)) + # in case we want the non-stacked model again
  geom_bar(stat = "identity", aes(x = sample_name, y = cleavage_per, 
                                  fill = factor(missed_cleavages, levels = c("3","2","1","0")))) +
  coord_cartesian(clip = "off", ylim = c(0,100)) +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "orange", "2" = "red", "3" = "black")) +
  # geom_text(data = cleavage_report_per,
  #           aes(y = cleavage_per + 5, x = missed_cleavages,
  #               label = paste0(round(x = cleavage_per, digits = 1), "%")), # in case we want the non-stacked model again
  #           size = 6) +
  geom_text_repel(data = cleavage_report_per,
            aes(y = 100 - cleavage_per/2, 
                x = sample_name,
                label = paste0(round(x = cleavage_per, digits = 1), "%")),
            size = 6) +
  ggnewscale::new_scale_fill() +
  # geom_rect(data = cleavage_report_per,
  #           aes(xmin = -Inf, xmax = Inf, ymin = 105, ymax = 110, fill = sample_type), 
  #           alpha = 0.2) +
  # scale_fill_manual(values = c("Hela" = "skyblue", "Sample" = "black")) +
  theme_bw() +
  ggtitle("Missed cleavages per sample") +
  xlab("# Missed cleavages") +
  ylab("Peptide %") +
  # facet_wrap(~factor(sample_name, levels = c(
  #   c(sort(unique(cleavage_report$sample_name)[-grep(pattern = "Hela_", x = unique(cleavage_report$sample_name))]))   # in case we want the non-stacked model again
  #   c(sort(unique(cleavage_report$sample_name)[grep(pattern = "Hela_", x = unique(cleavage_report$sample_name))])))),  
  #   scales = "free") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        strip.background = element_rect(fill = NA))

png(filename = "./plots/0_missed_cleavages.png", width = 1400, height = 1400)
missed_cleavages
dev.off()

## Hydrophobicity
# unique_peptides <- unique(peptide_info_long$Stripped.Sequence)
hydrophobicity_report <- tibble::tibble(
  # The peptides
  peptide = peptide_info_long$Stripped.Sequence,
  sample_name = peptide_info_long$sample_name,
  type_of_sample = peptide_info_long$type_of_sample,
  hydrophobicity = sapply(X = peptide_info_long$Stripped.Sequence, 
                          FUN = function(x){Peptides::hydrophobicity(seq = x, scale = "KyteDoolittle")})
)

hydro_hgram <- ggplot(data = hydrophobicity_report , aes(x = hydrophobicity)) +
  # geom_density(mapping = aes(color = type_of_sample), linewidth = 2)+
  geom_density(mapping = aes(color = sample_name), linewidth = 2)+
  theme_bw()+
  ggtitle(label = "Hydrophobicity, KyteDoolittle scale") +
  xlab("Hydrophobicity")+
  ylab("Density")+
  theme(legend.position= "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20),
        title = element_text(size = 22),
        legend.text = element_text(size = 14, face = "bold"))

png(filename = "./plots/0_hydrophobicity.png", width = 1000, height = 1000)
hydro_hgram
dev.off()

