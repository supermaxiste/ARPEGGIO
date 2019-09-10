########################### DMR Seq ###########################
# This script takes as input cov files from Bismark and outputs DMRs

library(dmrseq)
library(data.table)
library(BiocParallel)

# Four command line arguements are needed: first is the number of samples
# for the first species analyzed, second is the number of samples for the
# second species analyzed, third is the output name with extension and 
# fourth is the number of cores. After those four, the cov files from the
# two species you want to compare need to be added (as many as you have).
# The order of the cov files MUST be diploid species first and polyploid species second

comm_args <- commandArgs(trailingOnly = TRUE)

# First argument: number of samples for one species
samples1 <- comm_args[1]

# Second argument: number of samples for the other species
samples2 <- comm_args[2]

# Second argument: output name (with extension)
output <- comm_args[3]

# Third argument: number of cores
cores <- comm_args[4]

# All other arguments: cov files (in the right order)
samples <- as.integer(samples1) + as.integer(samples2)
sample_counter <- 0
for (i in 1:samples){
  if (!is.na(comm_args[i+4])){
    assign(paste0("file_", i), comm_args[i+4])
    sample_counter <- sample_counter + 1
  }
}

# Create vector of cov files

cov_files <- c()
for (i in 1:sample_counter){
  cov_files[i] <- get(paste0("file_", i))
}

# Read cov files

bismarkBSseq <- read.bismark(files = c(cov_files),
                             rmZeroCov = TRUE, 
                             strandCollapse = FALSE,
                             verbose = TRUE)

# Specify conditions
sampleNames = c(rep("par", as.integer(samples1)), rep("kam", as.integer(samples2)))
pData(bismarkBSseq)$Species <- sampleNames

# Filtering step
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")==0) == 0)
sample.idx <- which(pData(bismarkBSseq)$Species %in% c("par", "kam"))

bs.filtered <- bismarkBSseq[loci.idx, sample.idx]

register(MulticoreParam(cores))

# DMRseq function, normally takes around 1.5h
regions <- dmrseq(bs = bs.filtered, testCovariate = "Species")

# Save Robject

outputR <- paste0(substr(output, 1, nchar(output)-3), "Rdata")
save(regions, file = outputR)

#This took about 1.5 h
regions_dataframe <- as.data.frame(regions)

write.csv(regions_dataframe, file=output, quote = FALSE, row.names = FALSE, col.names = TRUE)
