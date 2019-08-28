## This script takes as input candidate DM regions from dmrseq
## The output is a BED file with significant DM regions only (q<0.05)

library(data.table)

comm_args <- commandArgs(trailingOnly = TRUE)

# First argument: dmrseq output
dmrseq_output <- comm_args[1]

# Second argument: output name
output_name <- comm_args[2]


# Read file

candidate_regions <- fread(dmrseq_output)

# Select only significant regions

sig_regions <- candidate_regions[candidate_regions$qval < .05]

# Select range of significant regions

sig_regions_bed <- sig_regions[,c("seqnames", "start", "end")]

# Some regions might have start > end, so we need to fix this

# Get index of wrong positions

wrong_pos_index <- which(sig_regions_bed$start > sig_regions_bed$end)

#select rows with wrong index

wrong_pos <- sig_regions_bed[wrong_pos_index,]

#switch starting and ending position in the original file

sig_regions_bed[wrong_pos_index,2] <- wrong_pos$end
sig_regions_bed[wrong_pos_index,3] <- wrong_pos$start

#To check whether anything is wrong use following command:
#wrong_pos_index <- which(sig_regions_bed$start>sig_regions_bed$end)

# output bed file

write.table(sig_regions_bed, file=paste0(output_name, ".bed"), quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)

command <- paste("sort -k1,1 -k2,2n ", output_name, ".bed > ", output_name, "_sorted.bed", 
                 sep = "")

system(command = command)
