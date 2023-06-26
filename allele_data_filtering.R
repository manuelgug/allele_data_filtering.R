library(dplyr)


args = commandArgs(trailingOnly=T)
path=args[1]
CFilteringMethod=args[2] # OPTIONS: max, q95, amp_max, amp_q95 
MAF=args[3] #minimum allele frequency; default 0.02


allele.data<-read.csv(path, sep ="\t")

#0) identify false positives (step not needed)
pos_controls_before_index <- grepl("(?i)3D7", allele.data$sampleID) & !grepl("(?i)(Dd2|HB3|PM)", allele.data$sampleID)
pos_controls_before <- allele.data[pos_controls_before_index, ]
multiple_alleles_before<-pos_controls_before[pos_controls_before$n.alleles > 1,] # all alleles from 3D7 samples must be 1 (single copy), more than that can be considered false positives

false_positives_before <- multiple_alleles_before %>%
  group_by(sampleID, locus) %>%
  filter(norm.reads.locus != max(norm.reads.locus)) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

print(paste("There were", dim(allele.data)[1]-dim(false_positives_before)[1], "alleles and", dim(false_positives_before)[1], "false positive alleles across", length(unique(false_positives_before$locus)), "amplicons from", length(unique(pos_controls_before$sampleID)), "`3D7 (monoclonal, single copy)` positive controls before applying any filter"))


#### CONTAMINATS FILTER ####

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

# apply contaminants filter
# set filter selected by user
if (CFilteringMethod == "max"){
  CFilteringMethod <- NEG_threshold_max
  }else if (CFilteringMethod == "q95"){
  CFilteringMethod <- NEG_threshold_q95
  }else if (CFilteringMethod == "amp_max"){
  CFilteringMethod <- NEG_thresholds_max
  }else{
    CFilteringMethod <- NEG_thresholds_q95}

# actually apply contaminants filter depending on user input
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

print(paste("There were", dim(filtered_allele.data)[1]-dim(false_positives)[1], "alleles and", dim(false_positives)[1], "false positive alleles across", length(unique(false_positives$locus)), "amplicons from", length(unique(false_positives$sampleID)), "`3D7 (monoclonal, single copy)` positive controls after applying the contaminants filter"))


#### MAF FILTER ####

# set default MAF filter to 2% if no value is provided (https://link.springer.com/content/pdf/10.1038/srep41108.pdf)
if (is.null(MAF) == TRUE || is.na(MAF) == TRUE){
  MAF <- 0.02
}

# apply MAF filter to remove potential false positives
filtered_allele.data <- filtered_allele.data[filtered_allele.data$norm.reads.locus > MAF, ]
# filtered_allele.data <- filtered_allele.data[, !(names(filtered_allele.data) %in% c("norm.reads.locus", "n.alleles"))] #remove old allele freqs and counts

# recalculate allele freqs for each sample based on remaining read counts & allele counts based on remaining alleles
# filtered_allele.data <- filtered_allele.data %>%
#   group_by(sampleID,locus) %>%
#   mutate(norm.reads.locus = reads/sum(reads))%>%
#   mutate(n.alleles = n())


# identify positive controls and false positive alleles from 3D7 (single copy) positive controls identify false positives (step not needed)
pos_controls_index <- grepl("(?i)3D7", filtered_allele.data$sampleID) & !grepl("(?i)(Dd2|HB3|PM)", filtered_allele.data$sampleID)
pos_controls <- filtered_allele.data[pos_controls_index, ]
multiple_alleles<-pos_controls[pos_controls$n.alleles > 1,] # all alleles from 3D7 samples must be 1 (single copy), more than that can be considered false positives

false_positives <- multiple_alleles %>%
  group_by(sampleID, locus) %>%
  filter(norm.reads.locus != max(norm.reads.locus)) # remove the most frequent allele of each amplicon from each sample and keep the others as false positives for filtering

print(paste("There are", dim(filtered_allele.data)[1]-dim(false_positives)[1], "alleles and", dim(false_positives)[1], "false positive alleles across", length(unique(false_positives$locus)), "amplicons from", length(unique(false_positives$sampleID)), "`3D7 (monoclonal, single copy)` positive controls after applying the contaminats and MAF filters"))


# export filtered allele data
write.table(filtered_allele.data,file="allele_data_filtered.txt",quote=F,sep="\t",col.names=T,row.names=F)
