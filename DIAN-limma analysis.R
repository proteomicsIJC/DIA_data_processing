############################################
### R for limma after data preprocessing ###
############################################

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

### functions importation
source("./functions/make_all_contrasts.R")
source("./functions/upsidedown.R")
source("./functions/xlsx_tt.R")

### data
## expression matrix
expression_matrix <- as.data.frame(readr::read_tsv("./results/expression_matrix_def.tsv"))
## annotation
annotation <- as.data.frame(readr::read_tsv("./results/annotation.tsv"))
## meta_data
meta_data <- as.data.frame(readr::read_tsv("./results/meta_data.tsv"))


#### data reshape
## expression matrix
expression_mat <- expression_matrix[,c(c(1), ## Accession
                                       c(8:ncol(expression_matrix)))] ## samples
row.names(expression_mat) <- expression_mat$Protein.Group
expression_mat <- expression_mat[,-1]
expression_mat <- as.matrix(expression_mat)

## annotation
rownames(annotation) <- annotation$Protein.Group

## meta_data
rownames(meta_data) <- meta_data$sample_name
meta_data <- meta_data[colnames(expression_mat),]
table(colnames(expression_mat) == rownames(meta_data))

##################
### limma work ###
##################

## Desing matrix
groups <- meta_data$group
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("^groups", "", colnames(design))
colnames(design) <- gsub(" ","_", colnames(design))
design

## Model fit
fit <- lmFit(expression_mat, design = design)

## Prepare the contrasts
contrasts_all <- make_all_contrasts(design = design, differentiating_element = " vs. ")

## Fit the contrasts
fit1 <- contrasts.fit(fit = fit, contrasts = contrasts_all)
fit1 <- eBayes(fit = fit1)

## Toptable xunga
tt <- topTable(fit = fit1, number = Inf)


## do the xlsx
xlsx_work <- xlsx_tt(fit__1 = fit1, 
                     meta_data = meta_data, meta_sample_column = "sample_name", meta_data_column = "group", 
                     annotation = annotation,
                     expression_matrix =  expression_mat,
                     differentation_element = " vs. ",
                     color_samples = c("CHO" = "darkorange",
                                       "UT" = "skyblue"),
                     filename = "./results/limma_example.xlsx")

