args <- commandArgs(trailingOnly = TRUE)
# Check if a file name is provided
if (length(args) == 0) {
  stop("No file name provided. Usage: Rscript script.R <filename>")
}

bim <- args[1]

bim2 <- args[2]

file <- read.csv(bim,header = FALSE, sep = "\t")

file2 <- read.csv(bim2,header = FALSE, sep = "\t")

common_snps <- intersect(file$V2, file2$V2)

write.csv(common_snps, file = "common_snps.txt",quote=FALSE, col.names = FALSE, row.names = FALSE)
