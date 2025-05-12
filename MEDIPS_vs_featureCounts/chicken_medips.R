# Load MEDIPS library
library(MEDIPS)
library("BSgenome.Ggallus.UCSC.galGal5")
library(doParallel)
library(foreach)

################################################################################################
################################################################################################
################################################################################################
##############define wour workspace where your bams are located#################################
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/")

##############define wour workspace################################################
# Load the list of filenames
filenames <- list.files(pattern = ".MQ10.bam$")
##############define wour output directory################################################has to be your wd plus "merged"
bam_dir= file.path("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/")
output_dir <- file.path("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/chicken_RBC/")
#manually
#output_dir <- file.path("/proj/gbs_medip/GEroNIMO/analysis/VI-3422/Aligned/gal7/HPT/merged")

##############define wour workspace################################################
# Define the maximum number of samples to run per subset
max_samples_per_subset <- 10
# Set the number of cores to use
# numCores <- 2
# cl <- makeCluster(2)
# registerDoParallel(cl)

################################################################################################
################################################################################################
################################################################################################

# Divide filenames into subsets
filename_subsets <- split(filenames, ceiling(seq_along(filenames) / max_samples_per_subset))
# Specify the input file, ROI, and other parameters

genome="BSgenome.Ggallus.UCSC.galGal5"
chr_all=paste("chr", c(1:28,30:33, "Z", "W"), sep="")#eliminating chr32, avoide the erro in XStringSet in evaluating the argument 'x' in selecting a method for function 'XStringSet': Error in ans[] <- x : replace
uniq=0
ROI <- read.delim(file.path(bam_dir, "merged/final_sorted.saf"), sep = "\t", header = T, stringsAsFactors = FALSE)

newname<-list("chr", "start", "end", "name")
# Get the number of rows in the ROI data frame
n_rows <- nrow(ROI)
# Create a new data frame with the desired column names and numbers from 1 to n_rows
ROI <- data.frame(ROI$Chr, ROI$Start, ROI$End, 1:n_rows, stringsAsFactors = FALSE)
names(ROI)<-newname
ROI <- subset(ROI, chr %in% chr_all)

  subset_results <- list()

  # Process each file in the subset
  for (file in filenames) {
    result <- MEDIPS.createROIset(file = file, BSgenome = genome, uniq = uniq, ROI = ROI, chr.select = chr_all, paired = TRUE, isSecondaryAlignment = FALSE, simpleCigar = TRUE)
    subset_results[[file]] <- result
  }


# Merge the results from all subsets in the same order
r <- unlist(subset_results, recursive = FALSE)

#save(r, file = "r_test.rda")
save(r, file = file.path(output_dir, "r_SecondaryAlignment_chicken.rda"))
rm(list = ls())


###run the necessary libraries and files again
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/")

##############define wour workspace################################################
# Load the list of filenames
filenames <- list.files(pattern = ".MQ10.bam$")
##############define wour output directory################################################has to be your wd plus "merged"
bam_dir= file.path("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/")
output_dir <- file.path("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/chicken_RBC/")
load(file = file.path(output_dir, "r_SecondaryAlignment_chicken.rda"))
library(MEDIPS)
library("BSgenome.Ggallus.UCSC.galGal5")


# Continue with the rest of your code

CS = MEDIPS.couplingVector(pattern = "CG", refObj = r[[1]])


methList <- list()  # Create an empty list to store the read assignments
for (i in seq_along(r)) {
  # Perform the MEDIPS.meth analysis
meth_result <- MEDIPS.meth(MSet1 = r[[i]], MSet2 = NULL, CSet = NULL, ISet1 = NULL, ISet2 = NULL, p.adj = NULL, diff.method = NULL, CNV = FALSE, MeDIP = FALSE)
#meth_result <- MEDIPS.meth(r[[i]])
  # Remove the statistics from the result
 # meth_result$statistics <- NULL
  # Store the read assignments in the list
  methList[[i]] <- meth_result
}

#save(methList, file="methList.rda")
save(methList, file = file.path(output_dir, "methListr_chicken_SecondaryAlignment.rda"))
#load(file="methList_SecondaryAlignment.rda")

###########matrix of counts################################
meth_countmatrix <- methList[[1]][, c("chr", "start", "stop")]

for (i in seq_along(methList)) {
  counts_column <- methList[[i]][[4]]
  column_name <- colnames(methList[[i]])[4]
  meth_countmatrix <- cbind(meth_countmatrix, counts_column)
  colnames(meth_countmatrix)[ncol(meth_countmatrix)] <- column_name
}

#save(meth_countmatrix, file="meth_countmatrix_SecondaryAlignment.rda")
#write.table(meth_countmatrix, "meth_countmatrix.txt" , sep="\t" , row.names = F , quote = F)
#save(meth_countmatrix, file = "meth_countmatrix_SecondaryAlignment.rda")
save(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix_chicken.rda"))

write.table(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix_chicken.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(meth_countmatrix, file = meth_countmatrix.txt, sep = "\t", row.names = FALSE, quote = FALSE)


################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
