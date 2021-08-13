## This file has been adapted from the file with the same name at
## https://github.com/markrobinsonuzh/diff_splice_paper

## ----- generate_flattened_gtf
## <<generate_flattened_gtf.R>>

## Generate manually flattened gtf files to use with featureCounts

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_gtf)
print(output_gtf)

library(GenomicRanges)
library(rtracklayer)

## Import gtf file
gtf0 <- import(input_gtf)

## Keep only exons
idx <- mcols(gtf0)$type == "exon" 
gtf0 <- gtf0[idx]
  
# Split by gene
gtf.g <- split(gtf0, mcols(gtf0)$gene_id)

# Flatten per gene
gtf.g.flat <- disjoin(gtf.g)

gtf.g.flat <- unlist(gtf.g.flat)
mcols(gtf.g.flat)$gene_id <- names(gtf.g.flat)
mcols(gtf.g.flat)$transcript_id <- mcols(gtf.g.flat)$gene_id
mcols(gtf.g.flat)$type <- "exon"

gtf.g.flat.sort <- gtf.g.flat[order(mcols(gtf.g.flat)$gene_id, decreasing = FALSE)]
exon.id <- split(mcols(gtf.g.flat.sort)$gene_id, mcols(gtf.g.flat.sort)$gene_id)

exon.id.new <- unlist(lapply(exon.id, function(g){ 
  seq(1, length(g)) 
}))

mcols(gtf.g.flat.sort)$exon_number <- exon.id.new
mcols(gtf.g.flat.sort)$exon_id <- paste0(mcols(gtf.g.flat.sort)$gene_id,
                                         ":", sprintf( "%03d", exon.id.new))
export(gtf.g.flat, output_gtf, format = "gtf")

