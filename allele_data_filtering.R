library(dplyr)
library(optparse)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Define and parse command-line arguments
option_list <- list(
  make_option(c("--allele_table"), type = "character", default = NULL, help = "Path to the allele table"),
  make_option(c("--resmarkers_table"), type = "character", default = NULL, help = "Path to the resmarkers table"),
  make_option(c("--CFilteringMethod"), type = "character", default = "global_max", help = "Contaminants filtering method: global_max, global_q95, amp_max, amp_q95"),
  make_option(c("--MAF"), type = "numeric", default = 0, help = "Minimum allele frequency; default 0"),
  make_option(c("--exclude_file"), type = "character", default = NULL, help = "Path to the file containing sampleIDs to exclude"),
  make_option(c("--use_case_amps"), type = "character", default = NULL, help = "Path to the file amplicons of your use case")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

allele_table <- opt$allele_table
resmarkers_table <- opt$resmarkers_table
CFilteringMethod <- opt$CFilteringMethod
MAF <- opt$MAF
exclude_file<- opt$exclude_file
use_case_amps<- opt$use_case_amps

#import allele table
allele.data<-read.csv(allele_table, sep ="\t")
amp_cats<-read.csv("amplicon_categories.csv", header =T)

# calculate initial read counts and allele freqs
allele.data <- allele.data %>%
  group_by(sampleID,locus) %>%
  mutate(norm.reads.locus = reads/sum(reads)) %>%
  mutate(n.alleles = n())

#add categories to allele.data
allele.data <- merge(allele.data, amp_cats[, c("locus.pool", "Category")], by.x = "locus", by.y = "locus.pool")

## 0) Exclude samples based on sampleIDs provided in a text file if the 'exclude' argument is provided
if (!is.null(exclude_file)) {
    remove_samples <- read.csv(exclude_file, sep = "\t", header = FALSE)
    allele.data <- subset(allele.data, !(sampleID %in% remove_samples$V1))
} 


# 1) calculate contaminants filtering thresholds from negative controls
if (any(grepl("(?i)Neg", allele.data$sampleID))) {
  
  neg_controls_index <- grepl("(?i)Neg", allele.data$sampleID)
  neg_controls <- allele.data[neg_controls_index, ]
  NEG_threshold_max <- max(neg_controls$reads) # max threshold across amplicons
  NEG_threshold_q95 <- quantile(neg_controls$reads, 0.95) # q95 threshold across amplicons
  
  all.loci<-unique(allele.data$locus)
  
  NEG_thresholds_max <- aggregate(reads ~ locus, data = neg_controls, FUN = function(x) max(x)) # max thresholds for each amplicon
  missing_loci <- setdiff(all.loci, NEG_thresholds_max$locus)
  missing_data <- data.frame(locus = missing_loci, reads = 0)
  NEG_thresholds_max <- rbind(NEG_thresholds_max, missing_data)
  
  NEG_thresholds_q95 <- aggregate(reads ~ locus, data = neg_controls, FUN = function(x) quantile(x, probs = 0.95)) # q95 thresholds for each amplicon
  missing_loci <- setdiff(all.loci, NEG_thresholds_q95$locus)
  missing_data <- data.frame(locus = missing_loci, reads = 0)
  NEG_thresholds_q95 <- rbind(NEG_thresholds_q95, missing_data)
  
} else {
  
  NEG_threshold_max <- 0
  NEG_threshold_q95 <- 0
  NEG_thresholds_max <- 0
  NEG_thresholds_q95 <- 0
}

#1) identify false positives
pos_controls_before_index <- grepl("(?i)3D7", allele.data$sampleID) & !grepl("(?i)(Dd2|HB3|PM)", allele.data$sampleID)
pos_controls_before <- allele.data[pos_controls_before_index, ]
multiple_alleles_before<-pos_controls_before[pos_controls_before$n.alleles > 1,] # all alleles from 3D7 samples must be 1 (single copy), more than that can be considered false positives

false_positives_initial <- suppressWarnings({
  multiple_alleles_before %>%
    group_by(sampleID, locus) %>%
    filter(norm.reads.locus != max(norm.reads.locus))}) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

#print(paste("There were", dim(allele.data)[1]-dim(false_positives_initial)[1], "alleles and", dim(false_positives_initial)[1], "false positive alleles across", length(unique(false_positives_initial$locus)), "amplicons from", length(unique(pos_controls_before$sampleID)), "`3D7 (monoclonal, single copy)` positive controls before applying any filter"))

initial_alleles<-dim(allele.data)[1]-dim(false_positives_initial)[1]
initial_false_positives_div<-sum(false_positives_initial$Category == "Diversity")
initial_false_positives_res<-sum(false_positives_initial$Category == "Resistance")
initial_false_positives_imm<-sum(false_positives_initial$Category == "Immune")
initial_false_positives_diag<-sum(false_positives_initial$Category == "Diagnostic")
initial_amplicons_w_false_positives<-length(unique(false_positives_initial$locus))
initial_controls_w_false_positives<-length(unique(false_positives_initial$sampleID))

initial_alleles_positive<-unique(false_positives_initial$asv)
initial_alleles_negative<-unique(neg_controls$asv)
initial_shared_contaminants_alleles<-length(intersect(initial_alleles_positive, initial_alleles_negative))

initial_amps_positive<-unique(false_positives_initial$locus)
initial_amps_negative<-unique(neg_controls$locus)
initial_shared_contaminants_amps<-length(intersect(initial_amps_positive, initial_amps_negative))


#### CONTAMINATS FILTER ####
CFilteringMethod_ <- CFilteringMethod #for exports

# apply contaminants filter set by user
if (CFilteringMethod == "global_max"){
  CFilteringMethod <- NEG_threshold_max
}else if (CFilteringMethod == "global_q95"){
  CFilteringMethod <- NEG_threshold_q95
}else if (CFilteringMethod == "amp_max"){
  CFilteringMethod <- NEG_thresholds_max
}else{
  CFilteringMethod <- NEG_thresholds_q95}

if (is.null(CFilteringMethod)) {
  print("No negative controls found. Skipping contaminants filter.")
  filtered_allele.data <- allele.data
} else {
  if (class(CFilteringMethod) =="data.frame") {
    allele.data_2 <- merge(allele.data, CFilteringMethod, by = "locus", suffixes = c("", ".NEG_threshold"))
    filtered_allele.data <- allele.data[allele.data_2$reads > allele.data_2$reads.NEG_threshold, ]
    filtered_allele.data <- filtered_allele.data[, !(names(filtered_allele.data) %in% c("norm.reads.locus", "n.alleles"))] #remove old allele freqs and counts
  } else {
    filtered_allele.data <- allele.data[allele.data$reads > CFilteringMethod, ]
    filtered_allele.data <- filtered_allele.data[, !(names(filtered_allele.data) %in% c("norm.reads.locus", "n.alleles"))] #remove old allele freqs and counts
  }
}

# recalculate allele freqs for each sample based on remaining read counts & allele counts based on remaining alleles
filtered_allele.data <- filtered_allele.data %>%
  group_by(sampleID,locus) %>%
  mutate(norm.reads.locus = reads/sum(reads))%>%
  mutate(n.alleles = n())

# identify positive controls and false positive alleles from 3D7 (single copy) positive controls (step not needed)
pos_controls_index <- grepl("(?i)3D7", filtered_allele.data$sampleID) & !grepl("(?i)(Dd2|HB3|PM)", filtered_allele.data$sampleID)
pos_controls <- filtered_allele.data[pos_controls_index, ]
multiple_alleles<-pos_controls[pos_controls$n.alleles > 1,] # all alleles from 3D7 samples must be 1 (single copy), more than that can be considered false positives

false_positives_1 <- suppressWarnings({
  multiple_alleles %>%
    group_by(sampleID, locus) %>%
    filter(norm.reads.locus != max(norm.reads.locus))}) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

#print(paste("There were", dim(filtered_allele.data)[1]-dim(false_positives_1)[1], "alleles and", dim(false_positives_1)[1], "false positive alleles across", length(unique(false_positives_1$locus)), "amplicons from", length(unique(false_positives_1$sampleID)), "`3D7 (monoclonal, single copy)` positive controls after applying the contaminants filter"))

contaminants_filter_alleles<-dim(filtered_allele.data)[1]-dim(false_positives_1)[1]
contaminants_filter_false_positives_div<-sum(false_positives_1$Category == "Diversity")
contaminants_filter_false_positives_res<-sum(false_positives_1$Category == "Resistance")
contaminants_filter_false_positives_imm<-sum(false_positives_1$Category == "Immune")
contaminants_filter_false_positives_diag<-sum(false_positives_1$Category == "Diagnostic")
contaminants_filter_amplicons_w_false_positives<-length(unique(false_positives_1$locus))
contaminants_filter_controls_w_false_positives<-length(unique(false_positives_1$sampleID))

contaminants_filter_alleles_positive<-unique(false_positives_1$asv)
contaminants_filter_alleles_negative<-unique(neg_controls$asv)
contaminants_filter_shared_contaminants_alleles<-length(intersect(contaminants_filter_alleles_positive, contaminants_filter_alleles_negative))

contaminants_filter_amps_positive<-unique(false_positives_1$locus)
contaminants_filter_amps_negative<-unique(neg_controls$locus)
contaminants_filter_shared_contaminants_amps<-length(intersect(contaminants_filter_amps_positive, contaminants_filter_amps_negative))


#### MAF FILTER ####

# apply MAF filter to remove potential false positives
filtered_allele.data <- filtered_allele.data[filtered_allele.data$norm.reads.locus > MAF, ]
filtered_allele.data <- filtered_allele.data[, !(names(filtered_allele.data) %in% c("n.alleles"))] #remove old allele counts

# recalculate allele counts based on remaining alleles
filtered_allele.data <- filtered_allele.data %>%
  group_by(sampleID,locus) %>%
  #   mutate(norm.reads.locus = reads/sum(reads))%>%
  mutate(n.alleles = n())

# identify positive controls and false positive alleles from 3D7 (single copy) positive controls identify false positives (step not needed)
pos_controls_index <- grepl("(?i)3D7", filtered_allele.data$sampleID) & !grepl("(?i)(Dd2|HB3|PM)", filtered_allele.data$sampleID)
pos_controls <- filtered_allele.data[pos_controls_index, ]
multiple_alleles<-pos_controls[pos_controls$n.alleles > 1,] # all alleles from 3D7 samples must be 1 (single copy), more than that can be considered false positives

false_positives_2 <- suppressWarnings({
  multiple_alleles %>%
    group_by(sampleID, locus) %>%
    filter(norm.reads.locus != max(norm.reads.locus))}) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

#print(paste("There are", dim(filtered_allele.data)[1]-dim(false_positives_2)[1], "alleles and", dim(false_positives_2)[1], "false positive alleles across", length(unique(false_positives_2$locus)), "amplicons from", length(unique(false_positives_2$sampleID)), "`3D7 (monoclonal, single copy)` positive controls after applying the contaminats and MAF filters"))

frequency_filter_alleles<-dim(filtered_allele.data)[1]-dim(false_positives_2)[1]
frequency_filter_false_positives_div<-sum(false_positives_2$Category == "Diversity")
frequency_filter_false_positives_res<-sum(false_positives_2$Category == "Resistance")
frequency_filter_false_positives_imm<-sum(false_positives_2$Category == "Immune")
frequency_filter_false_positives_diag<-sum(false_positives_2$Category == "Diagnostic")
frequency_filter_amplicons_w_false_positives<-length(unique(false_positives_2$locus))
frequency_filter_controls_w_false_positives<-length(unique(false_positives_2$sampleID))

frequency_filter_alleles_positive<-unique(false_positives_2$asv)
frequency_filter_alleles_negative<-unique(neg_controls$asv)
frequency_filter_shared_contaminants_alleles<-length(intersect(frequency_filter_alleles_positive, frequency_filter_alleles_negative))

frequency_filter_amps_positive<-unique(false_positives_2$locus)
frequency_filter_amps_negative<-unique(neg_controls$locus)
frequency_filter_shared_contaminants_amps<-length(intersect(frequency_filter_amps_positive, frequency_filter_amps_negative))


#### REPORT ####

report <- data.frame(
  initial = numeric(8),
  contaminants_filter = numeric(8),
  frequency_filter = numeric(8),
  row.names = c(
    "total_alleles",
    "false_positive_diversity_alleles",
    "false_positive_resistance_alleles",
    "false_positive_diagnostic_alleles",
    "false_positive_immunity_alleles",
    "amplicons_w_false_positive_alleles",
    "positive_controls_w_false_positives_alleles",
    "shared_contaminant_alleles"
  ),
  check.names = FALSE
)

report$initial <- c(initial_alleles, initial_false_positives_div, initial_false_positives_res, initial_false_positives_diag, initial_false_positives_imm, initial_amplicons_w_false_positives, initial_controls_w_false_positives, initial_shared_contaminants_alleles)
report$contaminants_filter <- c(contaminants_filter_alleles, contaminants_filter_false_positives_div, contaminants_filter_false_positives_res, contaminants_filter_false_positives_diag, contaminants_filter_false_positives_imm, contaminants_filter_amplicons_w_false_positives, contaminants_filter_controls_w_false_positives, contaminants_filter_shared_contaminants_alleles)
report$frequency_filter <- c(frequency_filter_alleles, frequency_filter_false_positives_div, frequency_filter_false_positives_res, frequency_filter_false_positives_diag, frequency_filter_false_positives_imm, frequency_filter_amplicons_w_false_positives, frequency_filter_controls_w_false_positives, frequency_filter_shared_contaminants_alleles)
report<-cbind(rownames(report), report)
colnames(report)[1] <- ""
colnames(report)[3]<-paste0("contaminants_filter_", CFilteringMethod_)
colnames(report)[4]<-paste0("frequency_filter_", MAF)

base_filename <- basename(allele_table)
filename <- tools::file_path_sans_ext(base_filename)

filtered_allele.data <- filtered_allele.data[, c(2, 1, 3:ncol(filtered_allele.data))]

write.table(filtered_allele.data,file=paste0(filename, "_", CFilteringMethod_, "_", as.character(MAF), "_filtered.csv"),quote=F,sep="\t",col.names=T,row.names=F)



## VISUALIZATION ##

report_2 <- report[1:5,]
colnames(report_2)[1]<-"alleles"
report_2$alleles<-c("total_alleles", "diversity", "resistance", "diagnostic", "immunity")

plot_data <- gather(report_2, key = "Category", value = "Value", -alleles)

# Define the desired order of columns for the x-axis
desired_order <- c(unique(plot_data$Category)[1], unique(plot_data$Category)[2], unique(plot_data$Category)[3])

# Reorder the rows based on the desired order of columns
plot_data$Category <- factor(plot_data$Category, levels = desired_order)

plot_total <- ggplot(plot_data[plot_data$alleles == "total_alleles", ], aes(x = Category, y = Value, group = 1)) +
  geom_line(linewidth = 1.5) +
  labs(title = "",
       x = "",
       y = "Total Alleles") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max(plot_data$Value) * 1.02))+  # Set the lower limit to 0
  theme(axis.text.x = element_text(angle = -10, hjust = 0.33))

plot_rest <- ggplot(plot_data[plot_data$alleles != "total_alleles", ], aes(x = Category, y = Value, group = alleles, color = alleles)) +
  geom_line(linewidth = 1.5) +
  labs(title = "",
       x = "",
       y = "False Positive Alleles") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, max(plot_data[plot_data$alleles != "total_alleles", ][,3]) + 1)) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = -10, hjust = 0.33))


### MICROHAPS FILTERING ###
if (!is.null(resmarkers_table)){
  microhaps<-read.csv(resmarkers_table, sep ="\t")
  
  # calculate read initial counts and allele freqs
  microhaps <- microhaps %>%
    group_by(sampleID,locus) %>%
    mutate(norm.reads.locus = Reads/sum(Reads)) %>%
    mutate(n.alleles = n())
  
  
  ## contaminants filtering
  if (class(CFilteringMethod)=="data.frame"){
    colnames(CFilteringMethod)[2]<-"Reads"
  }
  
  
  if (is.null(CFilteringMethod)) {
    print("No negative controls found. Skipping contaminants filter.")
    microhaps_filtered <- microhaps
  } else {
    if (class(CFilteringMethod) =="data.frame") {
      microhaps_2 <- merge(microhaps, CFilteringMethod, by = "locus", suffixes = c("", ".NEG_threshold"))
      microhaps_filtered <- microhaps[microhaps$Reads > microhaps_2$Reads.NEG_threshold, ]
      microhaps_filtered <- microhaps_filtered[, !(names(microhaps_filtered) %in% c("norm.reads.locus", "n.alleles"))] #remove old allele freqs and counts
    } else {
      microhaps_filtered <- microhaps[microhaps$Reads > CFilteringMethod, ]
      microhaps_filtered <- microhaps_filtered[, !(names(microhaps_filtered) %in% c("norm.reads.locus", "n.alleles"))] #remove old allele freqs and counts
    }
  }
  
  # recalculate read counts and allele freqs
  microhaps_filtered <- microhaps_filtered %>%
    group_by(sampleID,locus) %>%
    mutate(norm.reads.locus = Reads/sum(Reads)) %>%
    mutate(n.alleles = n())
  
  #frequency filtering
  microhaps_filtered <- microhaps_filtered[microhaps_filtered$norm.reads.locus > MAF, ]
  microhaps_filtered <- microhaps_filtered[, !(names(microhaps_filtered) %in% c("n.alleles"))] #remove old allele counts
  
  # recalculate allele counts based on remaining alleles
  microhaps_filtered <- microhaps_filtered %>%
    group_by(sampleID,locus) %>%
    #   mutate(norm.reads.locus = reads/sum(reads))%>%
    mutate(n.alleles = n())
  
  # EXPORT
  
  base_filename2 <- basename(resmarkers_table)
  filename2 <- tools::file_path_sans_ext(base_filename2)
  
  write.table(microhaps_filtered,file=paste0(filename2, "_", CFilteringMethod_, "_", as.character(MAF), "_filtered.csv"),quote=F,sep="\t",col.names=T,row.names=F)
}


### Specific use case report Are there amplicons from use_case on false_positives_initial, false_positives_1, false_positives_2?

if (!is.null(use_case_amps)){
  uc_amps <- read.csv("use_case.csv") #cambiar por use_case_amps!!!
  
  use_case_false_positives_initial <- sum(uc_amps$locus.pool %in% false_positives_initial$sampleID)
  use_case_false_positives_1 <- sum(uc_amps$locus.pool %in% false_positives_1$sampleID)
  use_case_false_positives_2 <- sum(uc_amps$locus.pool %in% false_positives_2$sampleID)
  
  use_case_row<-c("false_positives_use_case_amplicons", use_case_false_positives_initial, use_case_false_positives_1, use_case_false_positives_2)
  use_case_report <- rbind(report, use_case_row)
  rownames(use_case_report)[9]<-"false_positives_use_case_amplicons"
  
  #report with use case
  write.table(use_case_report,file=paste0(filename, "_", CFilteringMethod_, "_", as.character(MAF), "_filter_USE_CASE_report.csv"),quote=F,sep=",",col.names=T,row.names=F)
  
  use_case_report_2 <- use_case_report[9,]
  colnames(use_case_report_2)[1]<-"amps"
  use_case_report_2$amps<-c("use_case_amps")
  
  plot_data_2 <- gather(use_case_report_2, key = "Category", value = "Value", -amps)
  
  # Define the desired order of columns for the x-axis
  desired_order <- c(unique(plot_data_2$Category)[1], unique(plot_data_2$Category)[2], unique(plot_data_2$Category)[3])
  
  # Reorder the rows based on the desired order of columns
  plot_data_2$Category <- factor(plot_data_2$Category, levels = desired_order)
  plot_data_2<-plot_data_2[-4,]
  
  plot_use_case <- ggplot(plot_data_2, aes(x = Category, y = as.numeric(Value), group = amps, color = amps)) +
    geom_line(linewidth = 1.5) +
    labs(title = "",
         x = "",
         y = "False Positive Amplicons") +
    theme_minimal()+
    scale_y_continuous(limits = c(0, max(as.numeric(plot_data_2[,3])) + 1)) +
    theme(legend.title=element_blank())+
    theme(axis.text.x = element_text(angle = -10, hjust = 0.33))
  
  #plot with use case
  fig<-grid.arrange(plot_total, plot_rest, plot_use_case, ncol = 3)
  ggsave(paste0(filename, "_", CFilteringMethod_, "_", as.character(MAF), "_filter_USE_CASE_report.png"), fig, width = 21, height = 10.5, dpi = 300)
  
}else{
  #### EXPORTS WITHOUT USE CASE####
  
  #report
  write.table(report,file=paste0(filename, "_", CFilteringMethod_, "_", as.character(MAF), "_filter_report.csv"),quote=F,sep=",",col.names=T,row.names=F)
  
  #plot
  fig<-grid.arrange(plot_total, plot_rest, ncol = 2)
  ggsave(paste0(filename, "_", CFilteringMethod_, "_", as.character(MAF), "_filter_report.png"), fig, width = 14, height = 7, dpi = 300)
}
