#This script generates three cov files: one with CpG only, one with
#CHG only and one with CHH only. The input should be a Bismark coverage file
#with the follwing format:
## <scaffold> <position> <strand> <methylatedReads> <unmethylatedReads>
##  <context> <specificContext>
# First argument is the coverage file, second is output folder and third
# is output name

## Load data.table to read the file quickly
library(data.table)

comm_args <- commandArgs(trailingOnly = TRUE)
full_path <- normalizePath(comm_args[1])
## Read coverage file
coverage_file <- fread(full_path)
## Set working directory
out <- comm_args[2]
out <- normalizePath(out)
setwd(out)
# Save output name
output <- comm_args[3]
colnames(coverage_file) <- c("scaffold", "position", "strand", "mC", 
                             "uC", "context", "spec_context")
print("File has been read")

# Replace all NAs with 0s

coverage_file[is.na(coverage_file$mC),]$mC <- 0
coverage_file[is.na(coverage_file$uC),]$uC <- 0

# Save all positions with reads

with_reads <- (coverage_file$mC+coverage_file$uC)!=0

## We vectorize conditions to simplify our problem
## First check only contexts of interest: in our case CG, CHG and CHH (with reads)

CG_contexts <- coverage_file$context=="CG" & with_reads
CHG_contexts <- coverage_file$context=="CHG" & with_reads
CHH_contexts <- coverage_file$context=="CHH" & with_reads

## We now generate our cov files with the following format:
## <scaffold> <start_position> <end_position> <% methylation>
# <count_methylated> <count_unmethylated>

cov_CG <- cbind(coverage_file$scaffold[CG_contexts],
                coverage_file$position[CG_contexts],
                coverage_file$position[CG_contexts],
                (coverage_file$mC[CG_contexts]/(coverage_file$mC[CG_contexts]+coverage_file$uC[CG_contexts]))*100,
                coverage_file$mC[CG_contexts],
                coverage_file$uC[CG_contexts])

#Replace infinite values with 0 (not needed anymore)
#infinite <- cov_CG[,4]==Inf
#cov_CG[infinite,4] <- 0

cov_CHG <- cbind(coverage_file$scaffold[CHG_contexts],
                coverage_file$position[CHG_contexts],
                coverage_file$position[CHG_contexts],
                (coverage_file$mC[CHG_contexts]/(coverage_file$mC[CHG_contexts]+coverage_file$uC[CHG_contexts]))*100,
                coverage_file$mC[CHG_contexts],
                coverage_file$uC[CHG_contexts])

cov_CHH <- cbind(coverage_file$scaffold[CHH_contexts],
                 coverage_file$position[CHH_contexts],
                 coverage_file$position[CHH_contexts],
                 (coverage_file$mC[CHH_contexts]/(coverage_file$mC[CHH_contexts]+coverage_file$uC[CHH_contexts]))*100,
                 coverage_file$mC[CHH_contexts],
                 coverage_file$uC[CHH_contexts])

write.table(cov_CG, file = paste0(output,"_CG.cov"), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cov_CHG, file = paste0(output,"_CHG.cov"), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cov_CHH, file = paste0(output,"_CHH.cov"), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
