library(dplyr)


args = commandArgs(trailingOnly=T)
path=args[1]
CFilteringMethod=args[2] # OPTIONS: max, q95, amp_max, amp_q95 
MAF=args[3] #minimum allele frequency; default 0.02

allele.data<-read.csv(path, sep ="\t")
amp_cats<-read.csv("amplicon_categories.csv", header =T)

#add categories to allele.data
allele.data <- merge(allele.data, amp_cats[, c("locus.pool", "Category")], by.x = "locus", by.y = "locus.pool", all.x = TRUE)

## 0) Exclude samples based on sampleIDs provided in a text file if the 'exclude' argument is provided
if (!is.null(args[4])) {
  exclude_file <- args[4]
  if (file.exists(exclude_file)) {
    remove_samples <- read.csv(exclude_file, sep = "\t", header = FALSE)
    allele.data <- subset(allele.data, !(sampleID %in% remove_samples$V1))
  }
}


# 1) calculate contaminants filtering thresholds from negative controls
if (any(grepl("(?i)Neg", allele.data$sampleID))) {
  neg_controls_index <- grepl("(?i)Neg", allele.data$sampleID)
  neg_controls <- allele.data[neg_controls_index, ]
  NEG_threshold_max <- max(neg_controls$reads) # max threshold across amplicons
  NEG_threshold_q95 <- quantile(neg_controls$reads, 0.95) # q95 threshold across amplicons
  NEG_thresholds_max <- aggregate(reads ~ locus, data = neg_controls, FUN = function(x) max(x)) # max thresholds for each amplicon
  NEG_thresholds_q95 <- aggregate(reads ~ locus, data = neg_controls, FUN = function(x) quantile(x, probs = 0.95)) # q95 thresholds for each amplicon
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

false_positives_before <- multiple_alleles_before %>%
  group_by(sampleID, locus) %>%
  filter(norm.reads.locus != max(norm.reads.locus)) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

#print(paste("There were", dim(allele.data)[1]-dim(false_positives_before)[1], "alleles and", dim(false_positives_before)[1], "false positive alleles across", length(unique(false_positives_before$locus)), "amplicons from", length(unique(pos_controls_before$sampleID)), "`3D7 (monoclonal, single copy)` positive controls before applying any filter"))

initial_alleles<-dim(allele.data)[1]-dim(false_positives_before)[1]
initial_false_positives_div<-sum(false_positives_before$Category == "Diversity")
initial_false_positives_res<-sum(false_positives_before$Category == "Resistance")
initial_false_positives_imm<-sum(false_positives_before$Category == "Immune")
initial_false_positives_diag<-sum(false_positives_before$Category == "Diagnostic")
initial_amplicons_w_false_positives<-length(unique(false_positives_before$locus))
initial_controls_w_false_positives<-length(unique(pos_controls_before$sampleID))

initial_alleles_positive<-unique(false_positives_before$asv)
initial_alleles_negative<-unique(neg_controls$asv)
initial_shared_contaminants_alleles<-length(intersect(initial_alleles_positive, initial_alleles_negative))

initial_amps_positive<-unique(false_positives_before$locus)
initial_amps_negative<-unique(neg_controls$locus)
initial_shared_contaminants_amps<-length(intersect(initial_amps_positive, initial_amps_negative))


#### CONTAMINATS FILTER ####

# apply contaminants filter set by user
if (CFilteringMethod == "max"){
  CFilteringMethod <- NEG_threshold_max
}else if (CFilteringMethod == "q95"){
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

false_positives <- multiple_alleles %>%
  group_by(sampleID, locus) %>%
  filter(norm.reads.locus != max(norm.reads.locus)) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

#print(paste("There were", dim(filtered_allele.data)[1]-dim(false_positives)[1], "alleles and", dim(false_positives)[1], "false positive alleles across", length(unique(false_positives$locus)), "amplicons from", length(unique(false_positives$sampleID)), "`3D7 (monoclonal, single copy)` positive controls after applying the contaminants filter"))

contaminants_filter_alleles<-dim(filtered_allele.data)[1]-dim(false_positives_before)[1]
contaminants_filter_false_positives_div<-sum(false_positives$Category == "Diversity")
contaminants_filter_false_positives_res<-sum(false_positives$Category == "Resistance")
contaminants_filter_false_positives_imm<-sum(false_positives$Category == "Immune")
contaminants_filter_false_positives_diag<-sum(false_positives$Category == "Diagnostic")
contaminants_filter_amplicons_w_false_positives<-length(unique(false_positives$locus))
contaminants_filter_controls_w_false_positives<-length(unique(false_positives$sampleID))

contaminants_filter_alleles_positive<-unique(false_positives$asv)
contaminants_filter_alleles_negative<-unique(neg_controls$asv)
contaminants_filter_shared_contaminants_alleles<-length(intersect(contaminants_filter_alleles_positive, contaminants_filter_alleles_negative))

contaminants_filter_amps_positive<-unique(false_positives$locus)
contaminants_filter_amps_negative<-unique(neg_controls$locus)
contaminants_filter_shared_contaminants_amps<-length(intersect(contaminants_filter_amps_positive, contaminants_filter_amps_negative))


#### MAF FILTER ####

# set default MAF filter to 2% if no value is provided (https://link.springer.com/content/pdf/10.1038/srep41108.pdf)
if (is.null(MAF) == TRUE || is.na(MAF) == TRUE){
  MAF <- 0.02
}

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

false_positives <- multiple_alleles %>%
  group_by(sampleID, locus) %>%
  filter(norm.reads.locus != max(norm.reads.locus)) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

#print(paste("There are", dim(filtered_allele.data)[1]-dim(false_positives)[1], "alleles and", dim(false_positives)[1], "false positive alleles across", length(unique(false_positives$locus)), "amplicons from", length(unique(false_positives$sampleID)), "`3D7 (monoclonal, single copy)` positive controls after applying the contaminats and MAF filters"))

frequency_filter_alleles<-dim(filtered_allele.data)[1]-dim(false_positives)[1]
frequency_filter_false_positives_div<-sum(false_positives$Category == "Diversity")
frequency_filter_false_positives_res<-sum(false_positives$Category == "Resistance")
frequency_filter_false_positives_imm<-sum(false_positives$Category == "Immune")
frequency_filter_false_positives_diag<-sum(false_positives$Category == "Diagnostic")
frequency_filter_amplicons_w_false_positives<-length(unique(false_positives$locus))
frequency_filter_controls_w_false_positives<-length(unique(false_positives$sampleID))

frequency_filter_alleles_positive<-unique(false_positives$asv)
frequency_filter_alleles_negative<-unique(neg_controls$asv)
frequency_filter_shared_contaminants_alleles<-length(intersect(frequency_filter_alleles_positive, frequency_filter_alleles_negative))

frequency_filter_amps_positive<-unique(false_positives$locus)
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
colnames(report)[4]<-paste0("frequency_filter_", MAF)

#### EXPORTS ####

base_filename <- basename(path)
filename <- tools::file_path_sans_ext(base_filename)

write.table(filtered_allele.data,file=paste0(filename, "_filtered.csv"),quote=F,sep="\t",col.names=T,row.names=F)
write.table(report,file=paste0(filename, "_filter_report.csv"),quote=F,sep=",",col.names=T,row.names=F)
