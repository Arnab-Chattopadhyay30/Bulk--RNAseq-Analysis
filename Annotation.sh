#!/usr/bin/env Rscript

.libPaths(c("~/R/library", .libPaths()))

#----------------------------
# Input files
#----------------------------
gtf_file <- "annotation.gtf"
count_file <- "/transcriptomics/counts_matrix/count_matrix.txt"

#----------------------------
# Read GTF
#----------------------------
gtf <- read.delim(
  gtf_file,
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(gtf) <- c(
  "seqname","source","feature","start","end",
  "score","strand","frame","attribute"
)

# Keep only gene entries
genes <- gtf[gtf$feature == "gene", ]

#----------------------------
# Function to extract fields
#----------------------------
#"x" denotes the column of interest, as attribute column
#"field" denotes the section present in attribute column, eg: gene_name
#"sapply" does iteration of each rows and store the row in function(z)
#"strsplit" is splitting each rows stored in z by detecting ";"
#"trimws: does the trimming of the unwanted informations
 
extract_field <- function(x, field){

  sapply(x, function(z){

    pieces <- strsplit(z, ";")[[1]]
    pieces <- trimws(pieces)

    hit <- pieces[grep(paste0("^", field, " "), pieces)]

    if(length(hit)==0) return(NA)

    value <- sub(paste0(field, ' "'), "", hit)
    value <- sub('"$', "", value)

    value
  })

}

gene_map <- data.frame(
  gene_id   = extract_field(genes$attribute, "gene_id"),
  gene_name = extract_field(genes$attribute, "gene_name"),
  stringsAsFactors = FALSE
)

# Remove version numbers
gene_map$gene_id <- sub("\\..*$","",gene_map$gene_id)

# Remove duplicates
gene_map <- unique(gene_map)

cat("Genes in GTF:", nrow(gene_map), "\n")

#----------------------------
# Read count matrix
#----------------------------
counts <- read.delim(
  count_file,
  sep="\t",
  header=TRUE,
  skip=1,
  check.names=FALSE,
  stringsAsFactors=FALSE
)

counts$Geneid <- sub("\\..*$","",counts$Geneid)

cat("Genes in counts:", nrow(counts), "\n")

#----------------------------
# Merge
#----------------------------
final <- merge(
  counts,
  gene_map,
  by.x="Geneid",
  by.y="gene_id",
  all.x=TRUE,
  sort=FALSE
)

#----------------------------
# Reorder columns
#----------------------------
final <- final[, c(
  "Geneid",
  "gene_name",
  setdiff(names(final), c("Geneid","gene_name"))
)]

cat("Annotated genes:", sum(!is.na(final$gene_name)), "\n")
cat("Missing genes :", sum(is.na(final$gene_name)), "\n")

#----------------------------
# Save
#----------------------------
write.csv(
  final,
  "final_counts_annotated.csv",
  row.names=FALSE,
  quote=FALSE
)

cat("\nDone!\n")
cat("Output: final_counts_annotated.csv\n")