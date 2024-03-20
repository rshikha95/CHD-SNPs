# A computational approach to study SNPs associated with congenital heart disease
## Author: "Shikha Vashisht"
## date: "2024-01-05"

## Introduction: 
## This R code serves the purpose of systemetic SNP compilation, filtering based on Congenital Heart Disease (CHD) - associated traits, annotation, and integration with human heart organogenesis epigenetic data. The workflow includes the identification of potential CHD-associated enhancers, subsequent analysis for the identification of transcription factor binding sites (TFBS) within these enhancers, and conservation analysis. The script further aims to pinpoint regulatory CHD-SNPs (rCHD-SNPs) using human heart eQTL data. To enhance the robustness of the findings, a cross-validation step with enhancers from the human enhancer disease database is incorporated. The final step involves the visualization of the results, offering a holistic view of the potential regulatory elements associated with Congenital Heart Disease.

```{r}
# Install required packages

required_packages <- c("dplyr", "tidyr", "readr", "stringr", "data.table", "ggplot2", "GenomicRanges", "rtracklayer")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]


if(length(missing_packages) > 0) {
  install.packages(missing_packages)
  
}

print(missing_packages)

## Load required packages
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

## Step1: Download the SNPs data from GWAS-catalog and Clinvar database 
 
gwascatalog_url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
file_name <- "gwas_catalog_v1.0.2-associations_e110_r2023-07-29.tsv"
file_path <- "C:/Users/vashs/OneDrive/Documents/"
download.file(gwascatalog_url, destfile =paste0(file_path, file_name, sep=""))

clinvar_url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
file_name <- "variant_summary.txt.gz"
file_path <- "C:/Users/vashs/OneDrive/Documents/"
download.file(clinvar_url, paste0(file_path, file_name, sep=""))

# Step 2: Load the GWAS-catalog and ClinVar genetic variation files change to the directory where files are stored

setwd("~/CHD_Project/POSTDOC_PROJECT/new_analysis")
gwas<- read_tsv("gwas_catalog_v1.0.2-associations_e110_r2023-07-29.tsv", col_names = T)
head(gwas)
dim(gwas)

clinvar_variantSummary<- fread("variant_summary_clinvar.txt", header=T, sep="\t", stringsAsFactors = F)

head(clinvar_variantSummary)
dim(clinvar_variantSummary)

## Filter clinvar variations for hg38 assembly

clinvar_variantSummary_hg38<- clinvar_variantSummary[clinvar_variantSummary$Assembly == "GRCh38", ]

# Download CHD-traits list from Supplementary Table #1

setwd("~/CHD_Project/POSTDOC_PROJECT/new_analysis/Paper_final_things_for_submissions/supp_mat/Supplementary_Tables")

chd_traits<- read_xlsx("Supp_Table1.xlsx", sheet = 1, skip = 1)
nrow(chd_traits)
head(chd_traits)
chd_traits<- as.data.frame(chd_traits)
chd_traits<- chd_traits$Term
chd_traits_strings<- paste(chd_traits, collapse = "|")

## Step 3: Exclude CHD-unrelated traits.

exclusion_regex <- paste(c("Asthma", "down syndrome", "beta blocker", "body fat", "treatment", "hypertension", "diabetes", "treated", "exercise", "lipoprotein", "postoperative","Perioperative", "resting"), collapse = "|")

## Step 4: SNPs filtration based on CHD-traits

filt_gwas_chd<- gwas %>%
filter(stringr::str_detect(gwas$`DISEASE/TRAIT`, regex(chd_traits_strings, ignore_case = TRUE)) & !stringr::str_detect(`DISEASE/TRAIT`, regex(paste0("(?i)", exclusion_regex))))
filt_clinvar_chd<- clinvar_variantSummary_hg38 %>%
filter(stringr::str_detect(clinvar_variantSummary_hg38$PhenotypeList, regex(chd_traits_strings, ignore_case = TRUE)) & !stringr::str_detect(PhenotypeList, regex(paste0("(?i)", exclusion_regex))))

# Get unique CHD-SNPs from GWAS-catalog and ClinVar Database

 ## GWAS-catalog

filt_gwas_chd_uniq<- filt_gwas_chd %>% distinct(SNPS, .keep_all = TRUE)
nrow(filt_gwas_chd_uniq)

## ClinVar Database

filt_clinvar_chd_uniq <- filt_clinvar_chd %>% distinct(`RS# (dbSNP)`, .keep_all = TRUE)
nrow(filt_clinvar_chd_uniq)
table(filt_clinvar_chd_uniq$Type)

filt_clinvar_chd_SNV<- filt_clinvar_chd[filt_clinvar_chd$Type == "single nucleotide variant", ]
nrow(filt_clinvar_chd_SNV)

filt_clinvar_chd_SNV_uniq <- filt_clinvar_chd_SNV %>% distinct(`RS# (dbSNP)`, .keep_all = TRUE)
nrow( filt_clinvar_chd_SNV_uniq)

## Save files with unique CHD-SNPs from GWAS-catalog and ClinVar Database

write_tsv(filt_gwas_chd_uniq,"chd_filtered_SNPs/gwas_filtered_chd_uniq.txt")

write.table(file="clinvar_filtered_chd_SNV_uniq.txt", filt_clinvar_chd_SNV_uniq, col.names = T, row.names = F, quote=F)

## Prepare the input for ANNOVAR annotation 

# GWAS-catalog CHD-SNPs ANNOVAR Input

col14<- as.numeric(filt_gwas_chd_uniq$CHR_POS)
filt_gwas_chd_uniq$CHR_END<- col14
annovarInput_filt_gwas_chd_uniq<-filt_gwas_chd_uniq 
annovarInput_filt_gwas_chd_uniq$empty1<-rep(0,nrow(filt_gwas_chd_uniq))
annovarInput_filt_gwas_chd_uniq$empty2<-rep(0,nrow(filt_gwas_chd_uniq))

annovarInput_filt_gwas_chd_uniq_nw<-annovarInput_filt_gwas_chd_uniq[c(12:13,39:41,1:11,14:38)]
nrow(annovarInput_filt_gwas_chd_uniq_nw)

 # ClinVar CHD-SNPs ANNOVAR Input

annovar_filt_clinvar_chd_snv_uniq<- filt_clinvar_chd_SNV_uniq[,c(19:21,33:34, 1:18,22:32)]

# Save ANNOVAR input files 

write.table(file="annovarInput_gwas_filtered_chd_uniq.txt",annovarInput_filt_gwas_chd_uniq_nw,col.names = F, row.names = F, quote = F,sep = "\t")

write.table(file="annovarInput_clinvar_filtered_chd_snv_uniq.txt", annovar_filt_clinvar_chd_snv_uniq, col.names = F, row.names = F, quote = F,sep = "\t")

```

```{shell ANNOVAR annonation of CHD-SNPs from GWAS-catalog and ClinVar}
# Linux command to annotate CHD-SNPs using ANNOVAR
# Using ensembl annotations for hg38 genome assembly
# Install Annovar program 

annotate_variation.pl -out gwas_chd_filt_annovar_uniq -build hg38 annovarInput_gwas_filtered_chd_uniq.txt humandb/ -dbtype ensGene

annotate_variation.pl -out clinvar_chd_filt_SNV_annovar_engGene_uniq -build hg38 annovarInput_clinvar_filtered_chd_snv_uniq.txt humandb/ -dbtype ensGene
```

```{r ANNOVAR CHD-SNPs Annotations}
# Load ANNOVAR annotation results for GWAS-catalog and Clinvar

annovar_ensGene <- read_tsv("chd_filtered_SNPs/annovar_res/gwas_chd_filt_annovar_uniq.variant_function", col_names=F )
annovar_ensgene_context<-data.frame(table(annovar_ensGene$X1))
annovar_clinvar_chd_SNV_uniq_ensGene <- read_tsv("clinvar_chd_filt_SNV_annovar_engGene_uniq.variant_function", col_names=F )
annovar_clinvar_chd_SNV_uniq_ensGene_type<-data.frame(table(annovar_clinvar_chd_SNV_uniq_ensGene$X1))

```
```{r Segregate CHD-SNPs into noncoding and coding}

# Segregate CHD-SNPs into noncoding and coding 

library(dplyr)

# Define noncoding 

noncoding<- c("downstream", "upstream","intergenic", "intronic","ncRNA_intronic","upstream;downstream")

noncoding_gwas_uniq<- annovar_ensGene %>% filter(X1 %in% noncoding)

nrow(noncoding_gwas_uniq)
noncoding_gwas_uniq<- noncoding_gwas_uniq[,c(3:5,2,1,8:43)]

# Expand GWAS-catalog noncoding CHD-SNPs coordinates

gwas_noncoding_uniq_expand<- noncoding_gwas_uniq
gwas_noncoding_uniq_expand$X3<- paste("chr",gwas_noncoding_uniq_expand$X3, sep="")
gwas_noncoding_uniq_expand$snp_pos<-gwas_noncoding_uniq_expand$X4
gwas_noncoding_uniq_expand$X4<- as.numeric(gwas_noncoding_uniq_expand$X4)- 75
gwas_noncoding_uniq_expand$X5<- as.numeric(gwas_noncoding_uniq_expand$X5)+ 75
gwas_noncoding_uniq_expand<- gwas_noncoding_uniq_expand[,c(1:3,44,4:43)]

write_tsv(noncoding_gwas_uniq,"noncoding_gwas_chd_uniq")

write.table(file="noncoding_gwas_chd_uniq_expand_coords.bed",gwas_noncoding_uniq_expand, col.names = F, row.names = F, quote=F, sep="\t" )

# ClinVar Database

noncoding_clinvar_uniq<- annovar_clinvar_chd_SNV_uniq_ensGene %>% filter(X1 %in% noncoding)
nrow(noncoding_clinvar_uniq) 
noncoding_clinvar_uniq<- noncoding_clinvar_uniq[,c(3:5,2,1,8:43)]

# Expand ClinVar noncoding CHD-SNPs coordinates
clinvar_noncoding_uniq_expand<- noncoding_clinvar_uniq
clinvar_noncoding_uniq_expand$X3<- paste("chr",clinvar_noncoding_uniq_expand$X3, sep="")
clinvar_noncoding_uniq_expand$X4<- as.numeric(clinvar_noncoding_uniq_expand$X4)- 75
clinvar_noncoding_uniq_expand$X5<- as.numeric(clinvar_noncoding_uniq_expand$X5)+ 75 

write_tsv(noncoding_clinvar_uniq,"chd_filtered_SNPs/noncoding_clinvar_uniq")

write.table(file="chd_filtered_SNPs/clinvar_noncoding_uniq_expand_coords.bed",clinvar_noncoding_uniq_expand, col.names = F, row.names = F, quote=F, sep="\t" )

# Define coding

coding<- c("exonic","exonic;splicing", "ncRNA_exonic", "ncRNA_exonic;splicing")

# GWAS-catalog

coding_gwas_uniq<- annovar_ensGene %>% filter(X1 %in% coding)
nrow(coding_gwas_uniq)
genes_list_gwas<- unique(coding_gwas_uniq$X2)
split_genes_gwas <- unlist(strsplit(genes_list_gwas, ","))

write.table(file="exonicCoding_genesList_uniqGWAScat.txt",split_genes_gwas, col.names = F, row.names = F, quote = F )

# ClinVar database

coding_clinvar_uniq<- annovar_clinvar_chd_SNV_uniq_ensGene %>% filter(X1 %in% coding)
nrow(coding_clinvar_uniq)
genes_list_clinvar<-unique(coding_clinvar_uniq$X2)
split_genes_clinvar <- unlist(strsplit(genes_list_clinvar, ","))

write.table(file="exonic_coding_genes_list_uniqClinvar.txt",split_genes_clinvar, col.names = F, row.names = F, quote = F )

```

```{r Visualization: Distribution Plots of CHD-SNPs}

# GWAS-catalog CHD-SNPs

dist_gwas<- annovar_ensgene_context
dist_gwas$perc<- (dist_gwas$Freq/sum(dist_gwas$Freq)*100)
class(dist_gwas$perc)
dist_gwas$perc<- round(dist_gwas$perc,2)
dist_gwas$fill<- paste(dist_gwas$Var1,"=", dist_gwas$perc,"%", sep = "" )
dist_gwas <- dist_gwas %>% rename("Gwas-catalog CHD-SNPs" ="fill")

gwas_distPlot<-ggplot(dist_gwas, aes(x = 2, y = Freq, fill = `Gwas-catalog CHD-SNPs`)) + 
  geom_bar(stat = "identity", color = "white") + 
  coord_polar(theta = "y", start = 0, clip = "off") + 
  theme_void() +
  theme(
    text = element_text(size = 12 ),  # Adjust the font size as needed
    plot.title = element_text(hjust = 0.1)
  )

ggsave("gwas_CHD_SNPs_distPlot.tiff", plot = gwas_distPlot,
    dpi = 300, width = 7, height = 6, units = "in",type = "cairo", compression = "lzw")

# Clinvar CHD-SNPs

annovar_clinvar_chd_SNV_uniq_ensGene_type<-data.frame(table(annovar_clinvar_chd_SNV_uniq_ensGene$X1))
dist_clinvar<- annovar_clinvar_chd_SNV_uniq_ensGene_type
dist_clinvar$perc<- (dist_clinvar$Freq/sum(dist_clinvar$Freq)*100)
class(dist_clinvar$perc)
dist_clinvar$perc<- round(dist_clinvar$perc,2)
dist_clinvar$fill<- paste(dist_clinvar$Var1,"=", dist_clinvar$perc,"%", sep = "" )
dist_clinvar <- dist_clinvar %>%
rename("ClinVar CHD-SNPs" ="fill")

clinvar_distPlot <- ggplot(dist_clinvar, aes(x = 2, y = Freq, fill = `ClinVar CHD-SNPs`)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0, clip = "off") +
  theme_void() +
  theme(
    text = element_text(size = 12), 
    plot.title = element_text(hjust = 0.1)
  )
ggsave("clinvar_CHD_SNPs_distPlot.tiff", plot = clinvar_distPlot,
    dpi = 300, width = 7, height = 6, units = "in",type = "cairo", compression = "lzw")

```

```{shell Integration of Heart organogenesis epigenetic data from Cotney Lab}

# Step 1: Download the data files from Cotney's lab ftp server

dataLink="https://cotneylab.cam.uchc.edu/~jcotney/HEART_HUB/hg38/primary"

for sample in CS13_12383 CS13_12690 CS14_12408 CS14_14135 CS16_12997 CS16_14315 CS17_12291 CS17_12331 CS18_12059 CS18_12456 CS19_11914 CS19_12135 CS20_12448 CS20_12451 CS21_11849 CS21_12093 CS23_12058 CS23_12151; do
    for mark in H3K4me1 H3K27ac H3K4me3; do
        wget "${dataLink}/${sample}_${mark}.pval.signal.bigWig"
    done
done

# Step 2: Convert the file format of downloaded data from .bigWig to .wig and then convert .wig files to .bed

# .bigWig to .wig

for file in *.bigWig; do
    bigWigToWig "$file" "${file%.bigWig}.wig"
done


# .wig files to .bed

for file in *.wig; do
   convert2bed --input=wig --output=bed --zero-indexed < "$file" > "${file%.wig}.bed"
done  


# Step 3: Perform intersections of noncoding CHD-SNPs containing 150bp elements with  Coteny's data in Step 2

# Step 3a: Intersection with Noncoding elements with CHD-SNPs from GWAS-Catalog

for file in *.bed; do # process all the heart developmental stages with histone modification marks in Cotney's data
    main_part=$(basename "$file" | cut -d'.' -f1)
    output_name="${main_part}_gwasNonCod_intersect.bed"
    bedtools intersect -wo -a "$file" -b noncoding_gwas_chd_uniq_expand_coords.bed > "$output_name"
done

 # Step 3a: Intersection with Noncoding elements with CHD-SNPs from ClinVar Database
for file in *.bed; do
    main_part=$(basename "$file" | cut -d'.' -f1)
    output_name="${main_part}_clinvarNonCod_intersect.bed"
    bedtools intersect -wo -a "$file" -b clinvar_noncoding_uniq_expand_coords.bed > "$output_name"
done

```
```{r Correlation analysis}

# Step 4: Compute correlation between biological replicates from Cotney's data

stages <- c("CS13", "CS14", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")
marks <- c("H3K27ac", "H3K4me1", "H3K4me3")
rep1_files_gwas<- list.files("../../rep1/gwas", pattern="*_intersect.bed", full.names = T)
rep2_files_gwas<- list.files("../../rep2/gwas", pattern="*_intersect.bed", full.names = T)


# Step 4a: Initialize the correlation matrix

cor_matrix <- matrix(NA, nrow=length(stages), ncol=length(marks), 
                     dimnames=list(stages, marks))

# Step 4b: Iterate through each stage and mark, and compute the correlation

for (s in stages) {
    for (m in marks) {
        # Identify the matching files
        file1 <- grep(paste0(s, ".*", m), rep1_files_gwas, value=TRUE)
        file2 <- grep(paste0(s, ".*", m), rep2_files_gwas, value=TRUE)
        
        # Read data only if both files exist
        if (length(file1) > 0 && length(file2) > 0) {
            data1 <- fread(file1, sep="\t", stringsAsFactors=FALSE)
            data2 <- fread(file2, sep="\t", stringsAsFactors=FALSE)
            
            # Identify the common genomic coordinates
            gr1 <- GRanges(seqnames = data1$V1, ranges = IRanges(start = data1$V2, 
                   end =data1$V3), score = data1$V5)
            gr2 <- GRanges(seqnames = data2$V1, ranges = IRanges(start = data2$V2, 
                   end =  data2$V3), score = data2$V5)
            
            olaps <- findOverlaps(gr1, gr2)
            
 # Extract overlapping regions and their pvalSignal (or score) from both gr1 and gr2
            olap_gr1 <- gr1[queryHits(olaps)]
            olap_gr2 <- gr2[subjectHits(olaps)]
            
 # Compute the correlation for the overlapping regions
            if (length(olap_gr1) > 0) {
                cor_val <- cor(mcols(olap_gr1)$score, mcols(olap_gr2)$score)
                cor_matrix[s, m] <- cor_val
            }
        }
    }
}

```

```{r Visualization of correlated replicates from cotney data}

library(pheatmap)

colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Plot the heatmap

pheatmap(cor_matrix, 
         color = colors,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Correlation between rep1 and rep2",
         annotation_col = NULL, 
         border_color = "black")
```

```{r Merging files}

# Step 5: Combine all files resulted after all intersections of GWAS-catalog noncoding CHD-SNPs containing elements with cotney's data files and generate a combined file

library(data.table)

# List of all files to process resulted from GWAS-catalog noncoding elements with CHD-SNPs

rep1_files_gwas<- list.files("/rep1/gwas", pattern="*_intersect.bed", full.names = T)
rep2_files_gwas<- list.files("/rep2/gwas", pattern="*_intersect.bed", full.names = T)

# Function to extract stage and mark from each filename
get_stage_mark <- function(filename) {
  parts <- unlist(strsplit(basename(filename), "_"))
  list(stage = parts[1], mark = parts[3])
}

# Initialize an empty list to store data tables from each file

all_data <- list()

# Loop through each file, filter, augment and store data

for(file in rep1_files_gwas) { # similarly run this code for "rep2_files_gwas"
  data <- fread(file, stringsAsFactors = FALSE)
  
  # Filter based on threshold for column V5
  filtered_data <- data[V5 >= 9, ]
  
  # Extract stage and mark from the filename
  info <- get_stage_mark(file)
  
  # Add stage and mark columns to the data table
  filtered_data[, stage := info$stage]
  filtered_data[, mark := info$mark]
  
 # Append to the list
  all_data[[length(all_data) + 1]] <- filtered_data
}

# Combine all data tables into one

combined_data <- rbindlist(all_data)
combined_data_sub<- combined_data[,c(1:3,5,48:50, 6:10,23, 24, 19, 30,31, 36,37)]

# Write the combined data to a file

write.table(file="combined_files_gwas.bed",combined_data_sub,row.names = FALSE,
            col.names = FALSE, sep = "\t", quote = FALSE, ) 

# Process the combined file for 3 histone modifications marks (H3K27ac, H3K4me1 
# and H3K4me3)

data<- combined_data_sub

subsets <- list()

for (m in marks) {
    subsets[[m]] <- data[data$mark == m, ]
}

H3K27ac_data <- subsets[["H3K27ac"]]
H3K4me1_data <- subsets[["H3K4me1"]]
H3K4me3_data <- subsets[["H3K4me3"]]
table(H3K27ac_data$stage)
table(H3K4me1_data$stage)
table(H3K4me3_data$stage)

H3K27ac_subset_stages <-list()

   for (s in stages) {
        H3K27ac_subset_stages[[s]] <- H3K27ac_data[H3K27ac_data$stage == s,]
                 }

H3K4me1_subset_stages <- list()

   for (s in stages) {
        H3K4me1_subset_stages[[s]] <- H3K4me1_data[H3K4me1_data$stage == s,]
                 }

 H3K4me3_subset_stages <- list()

   for (s in stages) {
   H3K4me3_subset_stages[[s]] <- H3K4me3_data[H3K4me3_data$stage == s,]
                 }

# Step 6: Merge combined files in step5 for Biological rep 1 and rep 2 

# GWAS catalog Intersections with Cotney's data all stages and 3 histone modification Marks

# Combining rep1 and rep2 files

setwd("/Users/vashs/OneDrive/Documents/CHD_Project/POSTDOC_PROJECT/new_analysis/chd_filtered_SNPs/cotney_heart_epigenetic_data_intersections/intersection_with_unmergedBed_files/")

combined_data_sub_rep1<- fread("rep1/gwas/combined_files_gwas_rep1.bed", stringsAsFactors = F, header = F)
combined_data_sub_rep2<- fread("rep2/gwas/combined_files_gwas_rep2.bed", stringsAsFactors = F, header = F)

combined_data_gwas_rep1_2<- rbind(combined_data_sub_rep1, combined_data_sub_rep2)

enhancers <- unique(combined_data_gwas_rep1_2[, .(V1, V2, V3, V8, V9, V10, V11)])

# Stages list

stages<- c("CS13","CS14","CS16","CS17","CS18","CS19","CS20","CS21","CS23")
marks <- c("H3K27ac", "H3K4me1", "H3K4me3")

get_scores_and_marks <- function(enhan, data) {
  scores <- c()
  marks <- c()
  for(stage in stages) {
    stage_data <- data[V1 == enhan$V1 & V2 == enhan$V2 & V3 == enhan$V3 & V6 == stage]
    if(nrow(stage_data) > 0) {
      scores <- c(scores, unique(stage_data$V4))
      marks <- c(marks, unique(stage_data$V7))
    }
  }
  list(scores = scores, marks = marks)
}

# Initialize a storage box for GWAS-catalog 

storage_box <- data.table(enhancers)

 for(stage in stages) {
   storage_box[[stage]] <- 0
 }
   storage_box[, score := NA_character_]
   storage_box[, mark := NA_character_]

# Populate the storage box
   
for(i in 1:nrow(storage_box)) {
  enhancer <- storage_box[i]
  subset_data <- combined_data[V1 == enhancer$V1 & V2 == enhancer$V2 & V3 == enhancer$V3]
  
  sm <- get_scores_and_marks(enhancer, subset_data)
  storage_box[i, score := paste(sm$scores, collapse = ",")]
  storage_box[i, mark := paste(sm$marks, collapse = ",")]
  
  for(stage in stages) {
    if(nrow(subset_data[V6 == stage]) > 0) {
      storage_box[i, (stage) := 1]
    }
  }
}

nrow(storage_box)

write.table(file="combined_reps_all_stages_summary_gwas.txt", storage_box, col.names = T, row.names = F, quote = F, sep = "\t")

# Further refine and find unique putative enhancer entries

enhancer_stage_info_rep1_rep2_gwas <- storage_box[, 
                                   .(
                                     CS13 = max(CS13, na.rm = TRUE),
                                     CS14 = max(CS14, na.rm = TRUE),
                                     CS16 = max(CS16, na.rm = TRUE),
                                     CS17 = max(CS17, na.rm = TRUE),
                                     CS18 = max(CS18, na.rm = TRUE),
                                     CS19 = max(CS19, na.rm = TRUE),
                                     CS20 = max(CS20, na.rm = TRUE),
                                     CS21 = max(CS21, na.rm = TRUE),
                                     CS23 = max(CS23, na.rm = TRUE),
                                     scores = paste(score, collapse=","),
                                     marks = paste(mark, collapse=",")
                                   ), 
                                   by = .(V8, V9, V10, V11)
]

# View the result file
print(enhancer_stage_info_rep1_rep2_gwas)

# Calculate the number of stages the putative CHD-enhancer is present in

enhancer_stage_info_rep1_rep2_gwas[, num_stages := rowSums(.SD), .SDcols = c("CS13", "CS14", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")]

# Calculate the number of different marks associated with putative CHD-enhancer
enhancer_stage_info_rep1_rep2_gwas[, num_marks := lengths(lapply(strsplit(marks, ","), unique))]

# Extract annotations
annotation_data <- enhancer_stage_info_rep1_rep2_gwas[, .(num_stages, num_marks)]

# Save file
write.table(file="putCardiacEnh_allStagesSummary_reps.txt", enhancer_stage_info_rep1_rep2_gwas, col.names = T, row.names = F, quote = F, sep = "\t")

# Visualization

binary_matrix_rep1_rep2_gwas <- enhancer_stage_info_rep1_rep2_gwas[, c(5:13)]

sets_rep1_rep2_gwas<- colnames(enhancer_stage_info_rep1_rep2_gwas)[c(5:13)]

# color the intersection based on stages

library(ComplexUpset)
library(ggplot2)
library(RColorBrewer)

color_palette <- c(
  CS13 = "#e41a1c",
  CS14 = "#377eb8",
  CS16 = "#4daf4a",
  CS17 = "#984ea3",
  CS18 = "#ff7f00",
  CS19 = "#ffff33",
  CS20 = "#a65628",
  CS21 = "#f781bf",
  CS23 = "#999999"
)

png(filename = "putCardiacEnh_repsUpset_usorted.png", width = 5000, height = 2200, res = 300)
ComplexUpset::upset(
  binary_matrix_rep1_rep2_gwas,
  sets_rep1_rep2_gwas,
  queries=list(
    upset_query(set='CS13', fill=color_palette['CS13']),
    upset_query(set='CS14', fill=color_palette['CS14']),
    upset_query(set='CS16', fill=color_palette['CS16']),
    upset_query(set='CS17', fill=color_palette['CS17']),
    upset_query(set='CS18', fill=color_palette['CS18']),
    upset_query(set='CS19', fill=color_palette['CS19']),
    upset_query(set='CS20', fill=color_palette['CS20']),
    upset_query(set='CS21', fill=color_palette['CS21']),
    upset_query(set='CS23', fill=color_palette['CS23'])
  ),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        #counts=FALSE,
        bar_number_threshold=1,  # to show all numbers on top of bars
        width=0.5, text=list(
          vjust=-1.0001,
          hjust=-0.101, angle=10), text_colors = "blue"   
      )
      + scale_y_continuous(expand=expansion(mult=c(0, 0.1)))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=14),  
    colors=c('grey95', 'white')
  ),

  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.6))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
  sort_sets=FALSE, sort_intersections=FALSE, width_ratio=0.2)

dev.off()

# Heatmap

heatmap_data_rep1_re2 <- as.matrix(enhancer_stage_info_rep1_rep2_gwas[, 5:13])
png(filename = "putCardEnh_stageWiseHeatmap_gwas.png", width = 3000, height = 3000, res = 300)
heatmap(heatmap_data_rep1_re2, 
        main="Putative Cardiac Enhancers by human developmental Stage",
        xlab="Stages", 
        ylab="Enhancers", 
        Rowv=NA, Colv=NA,  # Prevents reordering of rows and columns
        col = c("white", "blue"), 
        scale="none", 
        margins=c(5,3))
dev.off()

# Group unique enhancers based on stage time points such as early, intermediate, late and their combinations

early_enh_gwas <- enhancer_stage_info_rep1_rep2_gwas[CS13 == 1 & CS14 == 1 & CS16 == 1]
nrow(early_enh)

intermed_enh_gwas <- enhancer_stage_info_rep1_rep2_gwas[CS17 == 1 & CS18 == 1 & CS19 == 1]
nrow(intermed_enh_gwas)

late_enh_gwas <- enhancer_stage_info_rep1_rep2_gwas[CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(late_enh_gwas)

always_active_gwas <- enhancer_stage_info_rep1_rep2_gwas[CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(always_active_gwas)

early_interm_enh_gwas <- enhancer_stage_info_rep1_rep2_gwas[CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 0 & CS21 == 0 & CS23 == 0]
nrow(early_interm_enh_gwas)

early_late_enh_gwas <- enhancer_stage_info_rep1_rep2_gwas[CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 0 & CS18 == 0 & CS19 == 0 & CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(early_late_enh_gwas)

early_late_enh_gwas
late_interm_enh_gwas <- enhancer_stage_info_rep1_rep2_gwas[CS13 == 0 & CS14 == 0 & CS16 == 0 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(late_interm_enh_gwas)

nrow(enhancer_stage_info_rep1_rep2_gwas)

list_of_best_enhancers_gwas <- list(early_enh = early_enh_gwas,
intermed_enh = intermed_enh_gwas,
late_enh = late_enh_gwas,
always_active = always_active_gwas,
early_interm_enh = early_interm_enh_gwas,
early_late_enh = early_late_enh_gwas,
late_interm_enh = late_interm_enh_gwas)
early_enh_gwas$group <- "early_enh"
intermed_enh_gwas$group <- "intermed_enh"
late_enh_gwas$group <- "late_enh"
always_active_gwas$group <- "always_active"
early_interm_enh_gwas$group <- "early_interm_enh"
early_late_enh_gwas$group <- "early_late_enh"
late_interm_enh_gwas$group <- "late_interm_enh"
best_enhancers_gwas_df <- rbind(early_enh_gwas, intermed_enh_gwas, late_enh_gwas, always_active_gwas, early_interm_enh_gwas, early_late_enh_gwas, late_interm_enh_gwas)

selected_enhancers_gwas <- enhancer_stage_info_rep1_rep2_gwas[num_stages >= 5]

palette <- brewer.pal(n = length(unique(summary_bestEnh_gwas$group)), name = "Set3")

png(filename = "topEnh_gwas_by_stageGrp.png", width = 1700, height = 1500, res = 300)
ggplot(summary_bestEnh_gwas, aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = palette) +
  labs(title = "Top Enhancers Distribution by Group", x = "Group", y = "Number of Enhancers") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

dev.off()

# GWAS-Catalog final putative cardiac enhancer table

enhancer_summaryTable_gwas <- left_join(enhancer_stage_info_rep1_rep2_gwas, combined_data_gwas_rep1_2, by = c("V8" = "V8", "V9" = "V9", "V10" = "V10", "V11" = "V11")) %>%
  group_by(V8, V9, V10, V11) %>%
  slice(1) %>%
  ungroup() %>%
  select(V8, V9, V10, V11, V12, V14, V15, V16, V17, V18, V19, num_stages, num_marks)

write.table(file="enhancerSumaryTab_gwas_final.txt", enhancer_summaryTable_gwas,col.names = F, row.names = F, quote = F, sep="\t")

```

```{r}
# Step 7: Combine all files resulted after all intersections of Clinvar noncoding CHD-SNPs containing elements with cotney's data files and generate a combined file

rep1_files_clinvar<-list.files("/rep1/clinvar", pattern="*_intersect.bed", 
                               full.names = T)

rep2_files_clinvar<-list.files("/rep2/clinvar", pattern="*_intersect.bed", 
                               full.names = T)

get_stage_mark <- function(filename) {
    parts <- unlist(strsplit(basename(filename), "_"))
    list(stage = parts[1], mark = parts[3])
}

all_data <- list()
for(file in rep1_files_clinvar) { # similarily run this loop for rep2 intersection files
    data <- fread(file, stringsAsFactors = FALSE)
    
  # Filter based on threshold for column V5 (MACS2 pvalSignal)
    
    filtered_data <- data[V5 >= 9, ]
    
  # Extract stage and mark from the filename
    
    info <- get_stage_mark(file)
    
  # Add stage and mark columns to the data table
    
    filtered_data[, stage := info$stage]
    filtered_data[, mark := info$mark]
    
  # Append to the list
    
    all_data[[length(all_data) + 1]] <- filtered_data
}

combined_data_rep1_clinvar <- rbindlist(all_data)
combined_data_rep1_clinvar_sub<- combined_data_rep1_clinvar[,c(1:3, 5, 42:44, 6:8, 41, 22, 11:12, 17, 26)]

combined_data_rep1_clinvar_sub$V22<- paste("rs",combined_data_rep1_clinvar_sub$V22, sep="" )

write.table(file="combined_files_clinvarRep1.bed",combined_data_rep1_clinvar,row.names = FALSE, sep = "\t", quote = FALSE, col.names = F) # Similarly save for rep2 files
write.table(file="combined_files_clinvarRep1_sub.bed",combined_data_rep1_clinvar_sub,row.names = FALSE, sep = "\t", quote = FALSE, col.names = F) # Similarly save for rep2 files

# Combine rep1 and rep2 files

combined_data_clinvar_rep1_2<- rbind(combined_data_rep1_clinvar_sub, combined_data_rep2_clinvar_sub)

colnames(combined_data_clinvar_rep1_2)<- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12","V13", "V14", "V15", "V16" )

# Create a unique list of enhancers

enhancers <- unique(combined_data_clinvar_rep1_2[, .(V1, V2, V3, V8, V9, V10, V11)])

# Stages list

stages<- c("CS13","CS14","CS16","CS17","CS18","CS19","CS20","CS21","CS23")
marks<-("H3K27ac", "H3K4me1", "H3K4me3")

get_scores_and_marks <- function(enhan, data) {
  scores <- c()
  marks <- c()
  for(stage in stages) {
    stage_data <- data[V1 == enhan$V1 & V2 == enhan$V2 & V3 == enhan$V3 & V6 == stage]
    if(nrow(stage_data) > 0) {
      scores <- c(scores, unique(stage_data$V4))
      marks <- c(marks, unique(stage_data$V7))
    }
  }
  list(scores = scores, marks = marks)
}

# Initialize storage box for clinvar

storage_box_clinvar <- data.table(enhancers)
for(stage in stages) {
  storage_box_clinvar[[stage]] <- 0
}
storage_box_clinvar[, score := NA_character_]
storage_box_clinvar[, mark := NA_character_]

# Populate the storage box

for(i in 1:nrow(storage_box_clinvar)) {
  enhancer <- storage_box_clinvar[i]
  subset_data <- combined_data_clinvar_rep1_2[V1 == enhancer$V1 & V2 == enhancer$V2 & V3 ==     enhancer$V3]
  
  sm <- get_scores_and_marks(enhancer, subset_data)
  storage_box_clinvar[i, score := paste(sm$scores, collapse = ",")]
  storage_box_clinvar[i, mark := paste(sm$marks, collapse = ",")]
  
  for(stage in stages) {
    if(nrow(subset_data[V6 == stage]) > 0) {
      storage_box_clinvar[i, (stage) := 1]
    }
  }
}

nrow(storage_box_clinvar)
write.table(file="combined_reps_all_stages_summ_clinvar.txt", storage_box_clinvar, col.names = T, row.names = F, quote = F, sep = "\t")


enhancer_stage_info_rep1_rep2_clinvar <- storage_box_clinvar[, 
                                             .(
                                                 CS13 = max(CS13, na.rm = TRUE),
                                                 CS14 = max(CS14, na.rm = TRUE),
                                                 CS16 = max(CS16, na.rm = TRUE),
                                                 CS17 = max(CS17, na.rm = TRUE),
                                                 CS18 = max(CS18, na.rm = TRUE),
                                                 CS19 = max(CS19, na.rm = TRUE),
                                                 CS20 = max(CS20, na.rm = TRUE),
                                                 CS21 = max(CS21, na.rm = TRUE),
                                                 CS23 = max(CS23, na.rm = TRUE),
                                                 scores = paste(score, collapse=","),
                                                 marks = paste(mark, collapse=",")
                                             ), 
                                             by = .(V6, V7, V8, V41)
]

enhancer_stage_info_rep1_rep2_clinvar[, num_stages := rowSums(.SD), .SDcols = c("CS13", "CS14", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")]

# Calculate the number of different marks associated with the enhancer

enhancer_stage_info_rep1_rep2_clinvar[, num_marks := lengths(lapply(strsplit(marks, ","), unique))]

# Extracting annotations

annotation_data_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[, .(num_stages, num_marks)]

# Save File

write.table(file="putCardiacEnh_allStagesSummaryclinvr.txt", enhancer_stage_info_rep1_rep2_clinvar, col.names = T, row.names = F, quote = F, sep = "\t")

# Visulaization 

binary_matrix_rep1_rep2_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[, c(5:13)]
sets_rep1_rep2_clinvar<- colnames(enhancer_stage_info_rep1_rep2_clinvar)[c(5:13)]

png(filename = "putativeCardiacEnh_plot_reps_clinvr.png", width = 2200, height = 2000, res = 300)

upset(binary_matrix_rep1_rep2_clinvar, sets = sets_rep1_rep2_clinvar,  matrix.color = "#56B4E9", keep.order = T)

dev.off()

png(filename = "putativeCardiacEnh_comUpset_reps.png", width = 4000, height = 2200, res = 600)
ComplexUpset::upset(
  binary_matrix_rep1_rep2_clinvar,
  sets_rep1_rep2_clinvar,
  queries=list(
    upset_query(set='CS13', fill=color_palette['CS13']),
    upset_query(set='CS14', fill=color_palette['CS14']),
    upset_query(set='CS16', fill=color_palette['CS16']),
    upset_query(set='CS17', fill=color_palette['CS17']),
    upset_query(set='CS18', fill=color_palette['CS18']),
    upset_query(set='CS19', fill=color_palette['CS19']),
    upset_query(set='CS20', fill=color_palette['CS20']),
    upset_query(set='CS21', fill=color_palette['CS21']),
    upset_query(set='CS23', fill=color_palette['CS23'])
  ),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,  # to show all numbers on top of bars
        width=0.4, text=list(
          vjust=-0.00001,
          hjust=-0.0001, angle=30)   # to adjust width of the bars
      )
      + scale_y_continuous(expand=expansion(mult=c(0, 0.01)))
      + theme(
        # to hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # to show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=14),  
    colors=c('grey95', 'white')
  ),
  
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.4))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
  sort_sets=FALSE, sort_intersections_by='cardinality',n_intersections=80, warn_when_dropping_groups=TRUE, width_ratio=0.08)
dev.off()

# Heatmap

heatmap_data_rep1_re2_clinvar <- as.matrix(enhancer_stage_info_rep1_rep2_clinvar[, 5:13])
png(filename = "../../clinvar/combined_rep1_ep2/putCardEnh_stageWiseHeatmap.png", width = 3000, height = 3000, res = 300)
heatmap(heatmap_data_rep1_re2_clinvar, 
        main="Putative Cardiac Enhancers by human developmental Stage",
        xlab="Stages", 
        ylab="Enhancers", 
        Rowv=NA, Colv=NA,  # Prevents reordering of rows/columns
        col = c("white", "blue"), 
        scale="none", 
        margins=c(5,5))
dev.off()

# Grouping

early_enh_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[CS13 == 1 & CS14 == 1 & CS16 == 1]
nrow(early_enh_clinvar)

intermed_enh_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[CS17 == 1 & CS18 == 1 & CS19 == 1]
nrow(intermed_enh_clinvar)
late_enh_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(late_enh_clinvar)
always_active_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(always_active_clinvar)
early_interm_enh_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 0 & CS21 == 0 & CS23 == 0]
nrow(early_interm_enh_clinvar)
early_late_enh_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 0 & CS18 == 0 & CS19 == 0 & CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(early_late_enh_clinvar)
early_late_enh_clinvar
late_interm_enh_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[CS13 == 0 & CS14 == 0 & CS16 == 0 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 1 & CS21 == 1 & CS23 == 1]
nrow(late_interm_enh_clinvar)
nrow(enhancer_stage_info_rep1_rep2_clinvar)
list_of_best_enhancers_clinvar <- list(early_enh = early_enh_clinvar,
                                    intermed_enh = intermed_enh_clinvar,
                                    late_enh = late_enh_clinvar,
                                    always_active = always_active_clinvar,
                                    early_interm_enh = early_interm_enh_clinvar,
                                    early_late_enh = early_late_enh_clinvar,
                                    late_interm_enh = late_interm_enh_clinvar)
early_enh_clinvar$group <- "early_enh"
intermed_enh_clinvar$group <- "intermed_enh"
late_enh_clinvar$group <- "late_enh"
always_active_clinvar$group <- "always_active"
early_interm_enh_clinvar$group <- "early_interm_enh"
early_late_enh_clinvar$group <- "early_late_enh"
late_interm_enh_clinvar$group <- "late_interm_enh"
best_enhancers_clinvar_df <- rbind(early_enh_clinvar, intermed_enh_clinvar, late_enh_clinvar, always_active_clinvar, early_interm_enh_clinvar, early_late_enh_clinvar, late_interm_enh_clinvar)

selected_enhancers_clinvar <- enhancer_stage_info_rep1_rep2_clinvar[num_stages >= 5]

palette <- brewer.pal(n = length(unique(summary_bestEnh_$group)), name = "Set3")

png(filename = "topEnh_gwas_by_stageGrp.png", width = 1700, height = 1500, res = 300)
ggplot(summary_bestEnh_gwas, aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = palette) +
  labs(title = "Top Enhancers Distribution by Group", x = "Group", y = "Number of Enhancers") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

dev.off()

# Create the summary table

enhancer_summaryTable_clinvar <- left_join(enhancer_stage_info_rep1_rep2_clinvar, combined_data_clinvar_rep1_2, by = c("V6" = "V8", "V7" = "V9", "V8" = "V10", "V41" = "V11")) %>% group_by(V6, V7, V8, V41) %>% slice(1) %>% ungroup() %>%
select(V6, V7, V8, V41, V12, V13, V14, V15, V16, num_stages, num_marks)

write.table(file="enhancerSumaryTab_clinvar_final.txt", enhancer_summaryTable_clinvar,col.names = F, row.names = F, quote = F, sep="\t")
```

```{r Putative CHD_enhancers file}

# Step 9: Create a combined "ONE MASTER CHD_enhancers file" by merging:
# 1. enhancer_stage_info_rep1_rep2_gwas and 2. enhancer_stage_info_rep1_rep2_clinvar 

nrow(enhancer_stage_info_rep1_rep2_gwas) # Final CHD-enhancers obtained from GWAS-Catalog
[1] 1139
nrow(enhancer_stage_info_rep1_rep2_clinvar) # Final CHD-enhancers obtained from ClinVar
[1] 917

enhancers_final_gwas<- enhancer_stage_info_rep1_rep2_gwas
enhancers_final_clinvar<- enhancer_stage_info_rep1_rep2_clinvar

colnames(enhancers_final_gwas)<- c("chr_enh", "start_enh","end_enh", "SNP_pos","CS13", "CS14", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23", "pval_Signal","histoneMod_Marks" , "num_stages", "num_marks")

colnames(enhancers_final_clinvar)<- c("chr_enh", "start_enh","end_enh", "SNP_pos","CS13", "CS14", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23", "pval_Signal","histoneMod_Marks" , "num_stages", "num_marks")

Master_enhancer_list_gwas_clinvar<- rbind(enhancers_final_gwas, enhancers_final_clinvar)
write.table(file="Master_enhancer_list_gwas_clinvar.txt", Master_enhancer_list_gwas_clinvar, col.names = T, row.names = F, quote = F, sep="\t")

early_enh_all<- Master_enhancer_list_gwas_clinvar %>% 
  filter(CS13 == 1 & CS14 == 1 & CS16 == 1)
nrow(early_enh_all)

intermed_enh_all <- Master_enhancer_list_gwas_clinvar %>% 
  filter (CS17 == 1 & CS18 == 1 & CS19 == 1)
nrow(intermed_enh_all)

late_enh_all <- Master_enhancer_list_gwas_clinvar %>% 
  filter(CS20 == 1 & CS21 == 1 & CS23 == 1)
nrow(late_enh_all)

always_active_all <- Master_enhancer_list_gwas_clinvar %>% filter(CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 1 & CS21 == 1 & CS23 == 1)
nrow(always_active_all)

early_interm_enh_all <- Master_enhancer_list_gwas_clinvar %>% filter (CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 0 & CS21 == 0 & CS23 == 0)
nrow(early_interm_enh_all)

early_late_enh_all <- Master_enhancer_list_gwas_clinvar %>% filter(CS13 == 1 & CS14 == 1 & CS16 == 1 & CS17 == 0 & CS18 == 0 & CS19 == 0 & CS20 == 1 & CS21 == 1 & CS23 == 1)
nrow(early_late_enh_all)

late_interm_enh_all <- Master_enhancer_list_gwas_clinvar %>% filter (CS13 == 0 & CS14 == 0 & CS16 == 0 & CS17 == 1 & CS18 == 1 & CS19 == 1 & CS20 == 1 & CS21 == 1 & CS23 == 1)
nrow(late_interm_enh_all)

list_of_best_enhancers_all <- list(early_enh = early_enh_all,
                                       intermed_enh = intermed_enh_all,
                                       late_enh = late_enh_all,
                                       always_active = always_active_all,
                                       early_interm_enh = early_interm_enh_all,
                                       early_late_enh = early_late_enh_all,
                                       late_interm_enh = late_interm_enh_all)
early_enh_all$group <- "early_enh"
intermed_enh_all$group <- "intermed_enh"
late_enh_all$group <- "late_enh"
always_active_all$group <- "always_active"
early_late_enh_all$group <- "early_late_enh"
best_enhancers_all_df <- rbind(early_enh_all, intermed_enh_all, late_enh_all, always_active_all, early_late_enh_all)

selected_enhancers_all<- Master_enhancer_list_gwas_clinvar[num_stages >= 5]

# Summarize the data by group

summary_bestEnh_all <- best_enhancers_all_df[, .(count = .N), by = group]
palette <- brewer.pal(n = length(unique(summary_bestEnh_all$group)), name = "Set3")

stageWise_enh_final <- ggplot(summary_bestEnh_all, aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, position = position_dodge(width = 0.9), size = 6) + scale_fill_manual(values = palette) + labs(title = "Enhancers Distribution by Developmental Group", x = "Group", y = "Number of Enhancers") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),  
    axis.title.x = element_text(face = "bold", size = 16), 
    axis.title.y = element_text(face = "bold", size = 16), 
    plot.title = element_text(face = "bold", size = 18, hjust = 0.3),
    plot.background = element_rect(fill = "white")  
  )

# Tiff Format
ggsave("histogramPlot_stageWise_FinalEnh.tiff", plot = stageWise_enh_final, dpi = 300, width = 12, height = 12, units = "in", type = "cairo", compression = "lzw")


# Make col names consistent

colnames(enhancer_stage_info_rep1_rep2)[1:4]<- c("chr","start","end","SNP_pos") 
colnames(enhancer_stage_info_rep1_rep2_clinvar)[1:4]<- c("chr","start","end","SNP_pos")

enhancer_stage_info_rep1_rep2_gwas_clinvar<- rbind(enhancer_stage_info_rep1_rep2,enhancer_stage_info_rep1_rep2_clinvar)
duplicated(enhancer_stage_info_rep1_rep2_gwas_clinvar)

# Check if any duplicated entries

which(!duplicated(enhancer_stage_info_rep1_rep2_gwas_clinvar)=="FALSE")
length(which(!duplicated(enhancer_stage_info_rep1_rep2_gwas_clinvar$SNP_pos)=="FALSE"))
nrow(enhancer_stage_info_rep1_rep2_gwas_clinvar)

# Visualization of all putative CHD-enhancers: GWAS-catalog and ClinVar

library(ComplexUpset)
library(UpSetR)

binary_matrix_rep1_rep2_gwas_clinvar <- enhancer_stage_info_rep1_rep2_gwas_clinvar[, c(5:13)]
sets_rep1_rep2_gwas_clinvar<- colnames(enhancer_stage_info_rep1_rep2_gwas_clinvar)[c(5:13)]

combined_upsetPlot_finalEnh <- ComplexUpset::upset(
  binary_matrix_rep1_rep2_gwas_clinvar,
  sets_rep1_rep2_gwas_clinvar,
  queries = list(
    upset_query(set = 'CS13', fill = color_palette['CS13']),
    upset_query(set = 'CS14', fill = color_palette['CS14']),
    upset_query(set = 'CS16', fill = color_palette['CS16']),
    upset_query(set = 'CS17', fill = color_palette['CS17']),
    upset_query(set = 'CS18', fill = color_palette['CS18']),
    upset_query(set = 'CS19', fill = color_palette['CS19']),
    upset_query(set = 'CS20', fill = color_palette['CS20']),
    upset_query(set = 'CS21', fill = color_palette['CS21']),
    upset_query(set = 'CS23', fill = color_palette['CS23'])
  ),
  base_annotations = list(
    'Intersection size' = (
      intersection_size(
        counts = FALSE, bar_number_threshold = 1,  # make counts=TRUE to show all numbers on top of bars
        width = 0.5, text = list(
          vjust = -0.01,
          hjust = -0.1101, angle = 15
        )   # to adjust width of the bars
      ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = 'black'),
          text = element_text(size = 19, face = "bold"),  
        )
    )
  ),
  stripes = upset_stripes(
    geom = geom_segment(size = 14),  
    colors = c('grey95', 'white')
  ),
  
  matrix = intersection_matrix(
    geom = geom_point(
      shape = 'circle filled',
      size = 3.5,
      stroke = 0.45
    )
  ),
  set_sizes = (
    upset_set_size(geom = geom_bar(width = 0.6)) +
      theme(
        axis.line.x = element_line(colour = 'black'),
        axis.ticks.x = element_line()
      )
  ),
  sort_sets = FALSE, sort_intersections = FALSE, width_ratio = 0.2
)

# Save plots 

ggsave("combinedUpsetPlot_finalEnh.tiff", plot = combined_upsetPlot_finalEnh, dpi = 300, width = 16, height = 12, units = "in", type = "cairo", compression = "lzw")

# Combined heatmap

par(bg = "white", cex.main = 5.5, cex.lab = 6, cex.axis = 6.5, font.axis=2)

tiff(filename = "putCardEnh_Heatmap_gwasClinvar.tif",width = 18,height = 14,unit="in",compression = "lzw",res = 300,type = "cairo")

heatmap(
  heatmap_data_rep1_re2_gwas_clinvar,
  main = "Presence of Putative Cardiac Enhancers in Human Developmental Stage",
  xlab = "Stages",
  ylab = "Enhancers",
  Rowv = NA,
  Colv = NA,
  col = c("grey", "blue"),
  scale = "none",
  margins = c(6, 3)
)

dev.off()

```
# Save the image of R session

save.image("combine_intersection_files.RData")

```{r}
# Retrieve DNA sequences for each putative CHD-enhancer for motif analysis

# Load required packages

library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(JASPAR2022)
library(TFBSTools)
library(parallel)

# Build background for the motif enrichment analysis

# Load the hg38 genome

genome <- Hsapiens

get_genome_frequencies <- function(genome, window_length=150) {
  # Extract individual chromosome sequences
  all_sequences <- lapply(seqnames(genome), function(chr) {
    as.character(genome[[chr]])
  })
}

genome_freq <- get_genome_frequencies(Hsapiens)  
write.table(file="background_model.txt", genome_freq, header=F, quote=F)

# Load Master list of CHD-enhancers

enh<- read.table("Master_PutativeCardiacEnhancers.txt",header=T,sep="\t", stringsAsFactors = F )
# Convert dataframe into a GRanges object

gr <- GRanges(seqnames=enh$chr_enh, ranges=IRanges(start=enh$start_enh, end=enh$end_enh))

make_header <- function(row) {
  paste0(">", enh$chr_enh[row], ":", enh$start_enh[row], "-", enh$end_enh[row])
}

sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
headers <- lapply(seq_len(nrow(enh)), make_header)

fasta_sequences <- mapply(function(header, sequence) {
  c(header, as.character(sequence))
}, headers, as.list(sequences), SIMPLIFY = FALSE)

# Save DNA sequences

writeLines(unlist(fasta_sequences), "all_Cardiac_putEnhancers.fasta")

```

```{shell}

# FIMO command

fimo --o fimo_analysis --bfile background_model.txt vertebrate_2022_CORE_JASPAR2_meme.txt all_Cardiac_putEnhancers.fasta

```

```{r}
# Load motif analysis results from FIMO

motifs<- read_tsv("fimo.tsv", col_names = TRUE)
filtered_motifs <- motifs[motifs$`q-value` < 0.05, ]
unique_Sigmotifs <- filtered_motifs[!duplicated(filtered_motifs$motif_alt_id), ]
nrow(unique_Sigmotifs)

# Save the siginificant motifs

write.table(file="Sig_Motifs_qval_0.05.txt", col.names = T, row.names = F, quote=F,
            sep="\t")


```

```{r Human Enhancer Disease Database (HEDD) intersection with CHD-enhancers}

library(GenomicRanges)

hedd_enh<- read.table("../../HEDD_intersec/HEDD_Enhancer_hg19.txt", header = T, sep="\t", stringsAsFactors = F)

nrow(hedd_enh)

hedd_enh$Length <- hedd_enh$End - hedd_enh$Start + 1

summary(hedd_enh$Length)

hist(hedd_enh$Length, main="Distribution of enhancers in HED Database", xlab="HEDD Enhancer Length", border="blue", col="lightgreen")

boxplot(hedd_enh$Length, ylab="Length")

hedd_enh1<- hedd_enh[,c(2:4,6,5,1)]
hedd_enh1$Chr<- paste("chr", hedd_enh1$Chr, sep="")
hedd_enh1_lftovr_input<- hedd_enh1[,c(1:3)]

# Download chain files for performing liftover of genomic coordinates of CHD-enhancers
# from hg38 to hg19 to perform the intersection with HEDD which has its data in hg19 genomic assembly. 

chainUrl <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

chainFile<- import.chain("hg19ToHg38.over.chain")

# Create a GRanges object for HEDD data

hedd_gr <- GRanges(seqnames = hedd_enh1_lftovr_input$Chr, 
              ranges = IRanges(start = hedd_enh1_lftovr_input$Start, end = hedd_enh1_lftovr_input$End))

# Perform LiftOver from hg38 to hg19 genomic assembly

hedd_gr_lifted <- liftOver(gr, chain)

# Filter out NULL entries (those that could not be lifted over)

gr_lifted <- gr_lifted[!sapply(gr_lifted, is.null)]

# Convert back to a data frame

df_lifted <- data.frame(Chr = seqnames(gr_lifted),
                        Start = start(gr_lifted),
                        End = end(gr_lifted),
                        Length = width(gr_lifted),
                        Source = mcols(gr_lifted)$Source,
                        EnhID = mcols(gr_lifted)$EnhID)

# Check the results

head(df_lifted)

# Save file

write.table(file="hedd_enh_lifted_hg19_hg38.txt", hedd_gr_lifted_df[,c(3:5)], col.names = T, row.names = F, sep="\t", quote=F)

```
```{shell}

# Bedtools intersect to find the overlapping enahncers with HEDD

bedtools intersect -wo -a hedd_enh_lifted_hg19_hg38.txt -b Master_PutativeCardiacEnhancers_withSNPs_Genes.txt > hedd_intersect_masterPutCardiacEnh_with_details.txt

```

```{r Conservation analysis}
# Conservation analysis of CHD-enhancers across 99 vertebrate species

# Download 99-way vertebrate phastCons conservation scores file from UCSC server using:
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/
  
# Download 99-way vertebrate phyloP conservation scores file from UCSC server using:
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw

# CHD-enhancers input

master_putCardiacENh<- fread("Master_PutativeCardiacEnhancers.txt", header = T, stringsAsFactors = F)

nrow(master_putCardiacENh)

master_putCardiacENh_withDetails<- read.table("Master_PutativeCardiacEnhancers_withSNPs_Genes.txt", header=T, stringsAsFactors = F, sep = "\t")

master_putCardiacENh_withDetails$metadata<- apply(master_putCardiacENh_withDetails, 1, function(row) paste(row, collapse = "-", sep=""))

cleaned_strings <- gsub("\\s+", "", master_putCardiacENh_withDetails$metadata)
PhastCons_conserAnalysis_input<- master_putCardiacENh_withDetails[,c(1:3)]

PhastCons_conserAnalysis_input$name<- cleaned_strings

# Save file

write.table(file="masterPutEnhInput.txt", PhastCons_conserAnalysis_input, col.names = F, row.names = F, quote = F, sep = "\t")

```
```{shell}
  
# Install "bigWigAverageOverBed" to compute 
# average conservation scores within putative cardiac enhancer sequences from # hg38.phastCons100way.bw
  
conda install -c bioconda ucsc-bigwigaverageoverbed 

bigWigAverageOverBed -minMax hg38.phastCons100way.bw masterPutEnhInput.txt conserved_putCardiacEnhancers.bed
```

```{r}

# SNPs are single base variations: we are also checking the conservation 
# of SNPs locations using PhyloP_100way vertebrate conservation data


SNP_positions_conservnAna<- data.frame(SNP_chr=master_putCardiacENh_withDetails$Chr, SNP_start=master_putCardiacENh_withDetails$SNP_pos, SNP_end=master_putCardiacENh_withDetails$SNP_pos+1, rsid=master_putCardiacENh_withDetails$rs_ID)

SNP_positions_conservnAna[1:10,]

write.table(file="SNPs_positionsInput.bed", SNP_positions_conservnAna, col.names = F, row.names = F, quote = F, sep = "\t")

```

```{r Random nucleotide sequence generation}

# 10000 random sequences of 150bp are generated from noncoding human genome and 
# compute conservation scores. The aim is to compare the conservation scores 
# computed for CHD-enhancers with these random DNA sequences scores.

# To generate the random sequences, chr1:22, chrM and chrX and chrY are utilized

# Load required packages

library(RMariaDB)
library(Biostrings)  
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)

getRandomSequences <- function(genome, txdb, numSequences, lengthSequence) {
  
  # Specify chromosomes
  
  chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
  
  allExons <- exons(txdb, columns="TXNAME", filter=list(TXCHROM=chroms))
  flankedExons <- flank(allExons, width=lengthSequence, both=TRUE)
  blockedRegions <- reduce(c(allExons, flankedExons))
  
  complementRegions <- setdiff(GRanges(seqnames = chroms, 
                                ranges = IRanges(start = 1, 
                                end = seqlengths(genome)[chroms])), blockedRegions)
  
  # Filter out regions smaller than lengthSequence
  
  validRegions <- complementRegions[width(complementRegions) >= lengthSequence]
  
  randomSeqDF <- data.frame(seqnames=character(0), 
                       start=integer(0), end=integer(0), name=character(0), 
                       stringsAsFactors=FALSE)
  
  for(i in 1:numSequences){
    sampledRegion <- sample(validRegions, 1)
    maxStart <- start(sampledRegion) + width(sampledRegion) - lengthSequence
    randStart <- sample(start(sampledRegion):maxStart, 1)
    randEnd <- randStart + lengthSequence - 1
    randomSeqDF <- rbind(randomSeqDF, data.frame(seqnames=seqnames(sampledRegion),
                      start=randStart, end=randEnd, name=paste0("randomSeq", i)))
  }
  
  # Write BED file
  
  write.table(randomSeqDF, file="randomSequences.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
  return(randomSeqDF)
}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Call function to generate sequences randomly

randomSeqs <- getRandomSequences(BSgenome.Hsapiens.UCSC.hg38, txdb, 10000, 150)


```

```{shell Computation of PhastCons scores for random noncoding sequences}

cd /mnt/c/Users/vashs/OneDrive/Documents/CHD_Project/POSTDOC_PROJECT/new_analysis/chd_filtered_SNPs/cotney_heart_epigenetic_data_intersections/intersection_with_unmergedBed_files/final_putativeCardiacEnh_gwasClinvar/conservation_analysis/random_nonCoding

bigWigAverageOverBed -minMax ../../../final_putativeCardiacEnh_gwasClinvar/conservation_analysis/hg38.phastCons100way.bw randomSequences.bed conserved_randomSeq.bed

```

```{r Statistical significance analysis on differences between CHD-enhancers and random noncoding sequences conservation scores}

randomSeq_cons<- read.delim("conserved_randomSeq.bed", header = F, stringsAsFactors = F, sep="\t")
cardiacEnh_cons<- read.delim("conserved_putCardiacEnhancers.bed", header = F, stringsAsFactors = F, sep="\t")

mean_cardEnh<-mean(cardiacEnh_cons$V5)
mean_randomSeq<- mean(randomSeq_cons$V5)
sd_cardEnh <- sd(cardiacEnh_cons$V5)
sd_randomSeq <- sd(randomSeq_cons$V5)

data <- data.frame(
PhastCons_Score = c(cardiacEnh_cons$V5, randomSeq_cons$V5),
Type = factor(rep(c("CardiacEnhancer", "Random"), times=c(length(cardiacEnh_cons$V5), length(randomSeq_cons$V5))))
)

conservationBoxPlot <- ggplot(data, aes(y = PhastCons_Score, x = Type, fill = Type)) +
  geom_boxplot() +
  theme_minimal()

conservationBoxPlot <- conservationBoxPlot +
  theme(
    axis.text.x = element_text(face = "bold", size = 16,family="arial"),
    axis.text.y = element_text(face = "bold", size = 16,family="arial"),  
    axis.title.x = element_text(face = "bold", size = 22,family="arial"), 
    axis.title.y = element_text(face = "bold", size = 22,family="arial"), 
    plot.background = element_rect(fill = "white"), 
    legend.key = element_rect(fill = "white")     
    )

# Save the plot

ggsave("boxplot_conservationAnalysis.tiff",plot = conservationBoxPlot,dpi = 300,
  width = 8, height = 10, units = "in", compression = "lzw")

```

```{r Integration of eQTL data from human heart tissue}

# Download the eQTL data files from GTEx database for human heart tissue 
# including "Artery_Aorta.v8.egenes", "Artery_Coronary.v8.egenes", "Heart_Atrial_Appendage.v8.egenes" and "Heart_Left_Ventricle.v8.egenes"

# Load required packages

library(dplyr)
library(tidyverse)

# Load CHD-enhancers file

enhancer_df<- read.delim("Master_PutativeCardiacEnhancers_withSNPs_Genes.txt", header = T, stringsAsFactors = F, sep="\t")

# Define a function to perform overlap with each eQTL file

overlap_eqtl <- function(eqtl_file_path, enhancer_df, file_name) {
  eqtl_df <- read.delim(eqtl_file_path, header = TRUE, stringsAsFactors = FALSE)
  
  merged_df <- merge(enhancer_df, eqtl_df, by.x="SNP_pos", by.y="variant_pos", all.x=TRUE)
  
  # Filter only the rows where there is a match
  
  matched_df <- merged_df[!is.na(merged_df$gene_id), ]
  
  # Add a column for the eQTL file where it matches
  
  matched_df$file_matched <- file_name
  
  return(matched_df)
}

# eQTL files

eqtl_files<- list.files(pattern="*.egenes.txt", full.names = T)

# Perform overlap for each file and store the results in a list

results <- lapply(eqtl_files, function(file) overlap_eqtl(file, enhancer_df, file))

# Combine all results

final_enh_heart_eQTLdf <- bind_rows(results)

# Save to a file

write.table(final_enh_heart_eQTLdf, "heart_eQTL_results.txt", sep="\t", row.names=FALSE)

```
