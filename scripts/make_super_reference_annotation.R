########################################################
## This script projects the coordinates of a genome gtf
## on to superTranscripts. It used the exon block positions
## in the superTranscript fasta file to get the coordinate
## transformation. If a transcript end up with adjacent 
## blocks on the superTranscripts, the blocks will be merged.
########################################################

library(GenomicRanges)

args <- (commandArgs(trailingOnly = TRUE))

#input and output file names
input_gtf<-args[1]
input_ST_fasta<-args[2] 
output_gtf<-args[3] 

# read the genome annotation gtf and get all the exon start and end positions
ann=read.delim(input_gtf,stringsAsFactors=F,header=F,skip=2)
pos1=ann$V4
pos2=ann$V5
chrom=ann$V1
extra_info=strsplit(ann$V9,";")
trans=gsub(" ","",gsub("transcript_id","",sapply(extra_info,function(x){x[2]})))
gene=gsub(" ","",gsub("gene_id","",sapply(extra_info,function(x){x[1]})))
genome_junctions=data.frame(chrom,pos1,pos2,trans,gene,stringsAsFactors=F)

## get the coordinate transform for each superTranscript from it's fasta ID.
command=paste("grep \">\" ",input_ST_fasta," | sed -e 's/:/\t/g' | sed -e 's/|/\t/g' | sed 's/>//g' > temp_list")
system(command)
blocks=read.delim("temp_list",stringsAsFactors=F,header=F)

blocks$V1<-gsub(" loc","",blocks$V1)
blocks$V5<-gsub(" segs","",blocks$V5)
blocks$V4<-gsub(" exons","",blocks$V4)

chrom_intervals=strsplit(blocks$V5,",")
block_intervals=strsplit(blocks$V6,",")
## reverse the block order if the strand is negative
block_intervals[blocks$V4=="-"]<-lapply(block_intervals[blocks$V4=="-"],rev)
lengths=sapply(chrom_intervals,length)
genes=rep(blocks$V1,lengths)
chroms=rep(blocks$V2,lengths)
strand=rep(blocks$V4,lengths)
chrom_intervals=unlist(chrom_intervals)
block_intervals=unlist(block_intervals)
get_bit<-function(interval,index){ as.numeric(sapply(strsplit(interval,"-"),function(x){x[index]})) } 
chrom_pos1=get_bit(chrom_intervals,1)
chrom_pos2=get_bit(chrom_intervals,2)
block_pos1=get_bit(block_intervals,1)

super_ref_intervals=data.frame(genes,chroms,chrom_pos1,chrom_pos2-chrom_pos1,block_pos1,stringsAsFactors=F)
colnames(super_ref_intervals)<-c("gene","chrom","chrom_start","width","block_start")
super_ref_intervals[strand=="-",]$block_start = super_ref_intervals[strand=="-",]$block_start + super_ref_intervals[strand=="-",]$width
super_ref_intervals$strand<-1
super_ref_intervals$strand[strand=="-"]<--1
gene_sri=split(super_ref_intervals,super_ref_intervals$gene)

## now work out the junction positions in the superTranscripts
get_ST_trans<-function(genome_junctions_list){
   i=1
   st_coords<-apply(genome_junctions_list,1,function(x){
      if(i %% 1000 == 0){ show(i) }
      csri=gene_sri[[x[5]]] #x[1] is chrom, x[5] is gene
      if(is.null(csri)){ return() }
      get_pos<-function(position){
         pos_diff = (position - csri$chrom_start)
         pos_ind=which(pos_diff <= csri$width & pos_diff >= 0)
         csri$block_start[pos_ind] + sign(csri$strand[pos_ind])*pos_diff[pos_ind]
      }
      start=get_pos(as.numeric(x[2])) #x[2] is pos1
      end=get_pos(as.numeric(x[3])) #x[3] is pos2
      if(start>end){ temp_start=start ; start=end ; end=temp_start }
      i<<-i+1
      IRanges(start,end)
   })
   message("making ranges object")
   st_coord_ranges=st_coords
   if(is.null(st_coord_ranges)){ return() }
   trans_ranges=split(unlist(IRangesList(st_coord_ranges)),genome_junctions_list$trans)
   message("reducing")
   unlist(IRangesList(sapply(trans_ranges,reduce)))
}
## convert the position chromosome-by-chromosome
## trying to do all at once takes too long
chrom_gj=split(genome_junctions,genome_junctions$chrom)
res=sapply(chrom_gj,get_ST_trans)
res=res[!sapply(res,is.null)]

trans_gene_map=unique(genome_junctions[,c("trans","gene")])
genes<-trans_gene_map$gene
names(genes)<-trans_gene_map$trans

trans=unlist(sapply(res,names))
contig=genes[trans]
start=unlist(sapply(res,start))
end=unlist(sapply(res,end))
n=length(trans)

gtf=data.frame(contig,rep("superTranscript",n),rep("exon",n),start,end,
  rep(".",n),rep("+",n),rep(".",n),
  paste("gene_id \"",contig,"\"; transcript_id \"",trans,"\"",sep=""))

## output the gtf of trancript for the superTranscriptome
write.table(gtf,output_gtf,quote=F,row.names=F,col.names=F,sep="\t")
