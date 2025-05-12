# Load MEDIPS library
library(MEDIPS)
library("BSgenome.Sscrofa.UCSC.susScr11")
library(doParallel)
library(foreach)

################################################################################################
################################################################################################
################################################################################################
##############define wour workspace where your bams are located#################################
#setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/pig_sperm/")
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/pig_sperm/")

##############define wour workspace################################################
# Load the list of filenames
filenames <- list.files(pattern = ".MQ10.bam$")
##############define wour output directory################################################has to be your wd plus "merged"
bam_dir= file.path("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/pig_sperm/")
output_dir <- file.path("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/pig_sperm/")


##############define wour workspace################################################
# Define the maximum number of samples to run per subset
max_samples_per_subset <- 10
# Set the number of cores to use
numCores <- 2
cl <- makeCluster(2)
registerDoParallel(cl)

################################################################################################
################################################################################################
################################################################################################

# Divide filenames into subsets
# filename_subsets <- split(filenames, ceiling(seq_along(filenames) / max_samples_per_subset))
# Specify the input file, ROI, and other parameters

genome="BSgenome.Sscrofa.UCSC.susScr11"

ROI <- read.delim(file.path(bam_dir, "merged/final_sorted.saf"), sep = "\t", header = T, stringsAsFactors = FALSE)
newname<-list("chr", "start", "end", "name")
# Get the number of rows in the ROI data frame
n_rows <- nrow(ROI)
# Create a new data frame with the desired column names and numbers from 1 to n_rows
ROI <- data.frame(ROI$Chr, ROI$Start, ROI$End, 1:n_rows, stringsAsFactors = FALSE)
names(ROI)<-newname
uniq=0
chr_all=paste("chr", c(1:18, "X", "Y"), sep="")
ROI <- subset(ROI, chr %in% chr_all)
subset_results <- list()
# Process each file in the subset
for (file in filenames) {
  result <- MEDIPS.createROIset(file = file, BSgenome = genome, uniq = uniq, ROI = ROI, paired = TRUE, isSecondaryAlignment = FALSE, simpleCigar = TRUE)
  subset_results[[file]] <- result
}


# Merge the results from all subsets in the same order
r <- unlist(subset_results, recursive = FALSE)

#save(r, file = "r_test.rda")
save(r, file = file.path(output_dir, "r_SecondaryAlignment_pig.rda"))
rm(list = ls())


###run the necessary libraries and files again
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/pig_sperm/")
output_dir <- file.path(getwd())
load(file = file.path(output_dir, "r_SecondaryAlignment_pig.rda"))
library(MEDIPS)
library("BSgenome.Sscrofa.UCSC.susScr11")


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
save(methList, file = file.path(output_dir, "methListr_pig_SecondaryAlignment.rda"))
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
save(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix_pig.rda"))

write.table(meth_countmatrix, file = file.path(output_dir, "meth_countmatrix_pig.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(meth_countmatrix, file = meth_countmatrix.txt, sep = "\t", row.names = FALSE, quote = FALSE)
