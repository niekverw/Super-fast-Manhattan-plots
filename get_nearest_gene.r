# install.
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# library(GenomicFeatures)
# library(IRanges)
# BiocManager::install("biomaRt")
library(biomaRt)
library(data.table)
library(matrixStats)
loadannotations <- function(){

  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  #listMarts()    # to see which database options are present
  #listDatasets(mart)     # function to see which datasets are present in ensembl

  mapping <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol",'band','chromosome_name',"start_position","end_position","gene_biotype","strand"), mart = mart)
  mapping.proteincoding <- mapping[mapping$gene_biotype=="protein_coding" & mapping$hgnc_symbol !="",]
  return(mapping.proteincoding)
}

findgenes <- function(mapping.proteincoding,chr,position,small_window=10000){
  # chr=12
  # position=2178165
  # chr=1;position=16339772
  chr=chr
  window=1000000
  genes.in.locus = data.frame()

  for (i in 1:10){
    startpos=position-window
    endpos=position+window

    startpos_smallwindow=position-small_window
    endpos_smallwindow=position+small_window
    genes.in.locus <- subset( mapping.proteincoding, mapping.proteincoding$chromosome_name==as.character(chr) & (
      ( mapping.proteincoding$start_position > startpos & mapping.proteincoding$start_position < endpos ) |
        ( mapping.proteincoding$end_position > startpos & mapping.proteincoding$end_position < endpos) ))

    window = window+1000000
    if (nrow(genes.in.locus) >=1){
      break
    }
  }

  genes.in.locus$DistanceSS <- genes.in.locus$start_position - position
  genes.in.locus$DistanceES <- genes.in.locus$end_position - position
  genes.in.locus$DistanceEE <- genes.in.locus$start_position - position
  genes.in.locus$DistanceSE<- genes.in.locus$end_position - position
  genes.in.locus$Distance <- as.numeric(rowMins(as.matrix(abs(genes.in.locus[,c("DistanceSS","DistanceES","DistanceEE","DistanceSE")]))))
  genes.in.locus$nearest.gene <- 0
  if(sum(genes.in.locus$DistanceES>0 & genes.in.locus$DistanceEE<0)>0){
    genes.in.locus[genes.in.locus$DistanceES>0 & genes.in.locus$DistanceEE<0,]$nearest.gene <- 1
  }
  if (max(genes.in.locus$nearest.gene)!=1){
    genes.in.locus[order(genes.in.locus$Distance),][1,]$nearest.gene <- 1
  }

  genes.in.locus$genes.in.locus <-0
  TF_genes.in.locus <- ( genes.in.locus$start_position > startpos_smallwindow & genes.in.locus$start_position < endpos_smallwindow ) |
    ( genes.in.locus$end_position > startpos_smallwindow & genes.in.locus$end_position < endpos_smallwindow)
  if (any(TF_genes.in.locus) ==TRUE){
    genes.in.locus[ TF_genes.in.locus ,]$genes.in.locus <-1
  }
  genes.in.locus$num_genes_1mb <- nrow(genes.in.locus)
  genes.in.locus = genes.in.locus[(genes.in.locus$nearest.gene ==1 | genes.in.locus$genes.in.locus ==1),]

  genes.in.locus$chr <- chr
  genes.in.locus$position <- position
  return(genes.in.locus[c('chr','position', 'ensembl_gene_id','hgnc_symbol','band', 'Distance','genes.in.locus','nearest.gene','num_genes_1mb')])

}
findgenes_df <- function(row,chr="CHR",position="hg19",snp="SNP"){
  #print(row["CHR"])
  #print(numeric(row$CHR))
  chr = as.numeric(unname(row[chr]))
  pos = as.numeric(unname(row[position]))
  snp = unname(row[snp])
  df.genes <- findgenes(mapping.proteincoding,chr=chr,position=pos,small_window=10000)
  df.genes$snp <- unname(snp)

  df.genes
}

process_snpstats_file <- function(snpstats){
  dfsnpstats <- as.data.frame(fread(file))
  df.anno <- do.call("rbind",apply(dfsnpstats,1, findgenes_df,chr="chromosome",position="position",snp="SNPID" ))
  #df.anno.agg <- aggregate(hgnc_symbol ~ snp + band, df.anno, paste, collapse = ",")

  df.anno.band <- aggregate(band ~ snp , unique(df.anno[,c("band","snp") ]), paste, collapse = ",")
  df.anno.gene <- aggregate(hgnc_symbol ~ snp , df.anno, paste, collapse = ",")
  df.anno.agg <- merge(df.anno.band,df.anno.gene,by=snp)
  names(df.anno.agg) <- c("SNPID","band","Gene")

  dfsnpstats <- merge(dfsnpstats,df.anno.agg,by="SNPID")
  return(dfsnpstats)
}

process_bolt_file <- function(file,chr="CHR",position="BP",snp="uniqid"){
  dfsnpstats <- as.data.frame(fread(file))
  df.anno <- do.call("rbind",apply(dfsnpstats,1, findgenes_df,chr=chr,position=position,snp=snp ))
  df.anno.band <- aggregate(band ~ snp , unique(df.anno[,c("band","snp") ]), paste, collapse = ",")
  df.anno.gene <- aggregate(hgnc_symbol ~ snp , unique(df.anno[,c("hgnc_symbol","snp") ]), paste, collapse = ",")
  df.anno.agg <- merge(df.anno.band,df.anno.gene,by="snp")
  names(df.anno.agg) <- c(snp,"band","Gene")

  dfsnpstats <- merge(dfsnpstats[,unique(colnames(dfsnpstats))],df.anno.agg,by=snp ) # deduplicate, sometimes uniqid is twice
  return(dfsnpstats)
}

get_genes_df <- function(df,chr="CHR",position="BP",snp="uniqid"){
  df<-data.frame(df)
  df.anno <- do.call("rbind",apply(df,1, findgenes_df,chr=chr,position=position,snp=snp ))
  df.anno.band <- aggregate(band ~ snp , unique(df.anno[,c("band","snp") ]), paste, collapse = ",")
  df.anno.gene <- aggregate(hgnc_symbol ~ snp , df.anno, paste, collapse = ",")
  df.anno.agg <- merge(df.anno.band,df.anno.gene,by="snp")
  names(df.anno.agg) <- c(snp,"band","Gene")
  df.genes <- merge(df[,unique(colnames(df))],df.anno.agg,by=snp ) # deduplicate, sometimes uniqid is twice
  return(df.genes)
}
#######################
#
# mapping.proteincoding <- loadannotations()
# df <- data.frame(fread("/Users/niekverw/Downloads/untitled.txt"))
#
# # 1 snp:
# findgenes(mapping.proteincoding,chr=1,position=16299312,small_window=10000)
#
# # dataframe
# df.anno <- do.call("rbind",apply(df,1, findgenes_df ))
#
# aggregate(hgnc_symbol ~ snp + band, df.anno, paste, collapse = ",")
#
#
# dfsnpstats <- process_snpstats_file("/Users/niekverw/Dropbox/Gwasshared_2/MR_iron/beverborg_grs/chrALL.snpstats")

args = commandArgs(trailingOnly=TRUE)

tsv <- args[1]
chrcol="CHR"
bpcol="BP"
uniqidcol="uniqid"
if(length(args)>1){
  chrcol=args[2]
  bpcol=args[3]
  uniqidcol=args[4]
}
if (!exists("mapping.proteincoding")){
  print("loading gene list from biomart (data.frame)")

  f.mapping.proteincoding = paste0(Sys.getenv("bitbuckethome"),"/Gwas-annotateStuff/nearbygenes/mapping.proteincoding.tsv")
  # mapping.proteincoding <- loadannotations()
  # fwrite(mapping.proteincoding,paste(f.mapping.proteincoding,row.names = FALSE,sep="\t",quote=FALSE  )
  try(mapping.proteincoding <- data.frame(fread(f.mapping.proteincoding)))

}
if (!exists("mapping.proteincoding")){
  print("loading gene list from biomart (online)")
  mapping.proteincoding <- loadannotations()
}


df <- process_bolt_file(tsv,chr=chrcol,position=bpcol,snp=uniqidcol)
print(paste0("writing: ", tsv,".anno"))
fwrite(df,paste0(tsv,".anno"),quote = FALSE,row.names = FALSE,sep = "\t")
