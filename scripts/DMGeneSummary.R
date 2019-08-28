## This script takes as input 1) the file with the intersection 
## between annotation and DM regions,
## 2) the DM regions file (from dmrseq) and 3) the geneID keyword
## in the 9th column of the annotation file.

library(data.table)
library(stringr)

comm_args <- commandArgs(trailingOnly = TRUE)

# First argument: intersection file

intersection <- comm_args[1]

# Second argument: DM regions (dmrseq)

regions <- comm_args[2]

# Third argument: geneID name

geneID <- comm_args[3]

# Fifth argument: output name

output_name <- comm_args[4]

# Read, clean and name file columns

intersection_file <- fread(intersection)
dm_regions <- fread(regions)

# Pick columns of interest
intersection_file <- intersection_file[, c(1, 3, 4, 5, 9, 11, 12, 13)]

# Rename columns
colnames(intersection_file) <- c("seqname", "feature", "start", "end", 
                                 "attribute", "overlap_start", "overlap_end",
                                 "length")

# Select only genes
intersection_file <- intersection_file[intersection_file$feature=="gene",]

# Pick geneID from attribute column
geneID_name <- paste0(geneID, "=") 

# Separate subcolumns in attribute column
attribute_col <- as.data.frame(str_split_fixed(intersection_file$attribute, ";", n=Inf))

# Look for column with geneID keyword
geneID_column <- which(grepl(geneID_name, unlist(attribute_col[1,])))

# Remove geneID keyword from column
intersection_file$attribute <- gsub(geneID_name, "", attribute_col[,geneID_column])

# Create final file (missing 1 column)
DM_genes_summary <- as.data.frame(cbind(geneID=intersection_file$attribute, 
                                        seqname=intersection_file$seqname,
                                        start=intersection_file$start, 
                                        end=intersection_file$end,
                                        region_start=intersection_file$overlap_start, 
                                        region_end=intersection_file$overlap_end, 
                                        overlap_length=intersection_file$length))

# Match DM regions coordinates

dm_coordinates <- paste(DM_genes_summary$seqname, ":", DM_genes_summary$region_start,
                        "-", DM_genes_summary$region_end, sep = "")
dm_regions <- cbind(dm_regions, coordinates=paste(dm_regions$seqname, ":", dm_regions$start,
                                      "-", dm_regions$end, sep = ""))
# Find which line of the summary corresponds to which line of the dmrseq file
corresponding_match <- match(dm_coordinates, dm_regions$coordinates)

# Add column to dmrseq with methylation status

m_status <- ifelse(dm_regions$stat>0, "decrease", "increase")

# Add final column to summary file

DM_genes_summary <- cbind(DM_genes_summary, methylation_status=m_status[corresponding_match])

# write file to output

write.table(DM_genes_summary, file=paste0(output_name, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)
