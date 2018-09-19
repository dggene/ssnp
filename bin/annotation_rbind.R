#coding=utf-8
#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
make_option(c("-r", "--rawbed"), type="character", default="input.txt",
              help="Input table file to read [default %default]"),
make_option(c("-w", "--gwasdb"), type="character", default="input.txt",
              help="gwasdb file to read [default %default]"),
make_option(c("-v", "--gwava"), type="character", default="input.txt",
              help="gwava file to read [default %default]"),
make_option(c("-b", "--sce"), type="character", default="input.txt",
              help="sce file to read [default %default]"), 
make_option(c("-e", "--eigen"), type="character", default="input.txt",
              help="Eigen file to read [default %default]"),
make_option(c("-s", "--spidex"), type="character", default="input.txt",
              help="spidex file to read [default %default]"),
make_option(c("-a", "--annovar"), type="character", default="input.txt",
              help="annovar to read [default %default]"),
make_option(c("-c", "--cadd"), type="character", default="input.txt",
              help="cadd file to read [default %default]"),
make_option(c("-i", "--silva_res"), type="character", default="input.txt",
              help="silva file to read [default %default]"),
make_option(c("-m", "--silva_mat"), type="character", default="input.txt",
              help="silva file to read [default %default]"), 
make_option(c("-p", "--primateai"), type="character", default="input.txt",
              help="primateai file to read [default %default]"),                            
make_option(c("-t", "--rnascore"), type="character", default="input.txt",
              help="rnascore file to read [default %default]")              
                            
)

opts = parse_args(OptionParser(option_list=option_list))
save.image('test.Rdata')

raw_res=read.table(opts$rawbed,sep='\t')[,c(1:4,6,7)]
if(length(grep('chr',raw_res$V1))>0){
    raw_res$V1=gsub('chr','',raw_res$V1)
    }
raw_res_by = transform(raw_res,merge=paste(raw_res$V1,raw_res$V3,raw_res$V6,raw_res$V7,sep='_'),merge1=paste(raw_res$V1,raw_res$V3,sep='_'))

gwasdb_res <- read.table(opts$gwasdb,head=T,check.names=F,sep='\t')
gwasdb_res_by<- transform(gwasdb_res,merge1=paste(gwasdb_res$Chr,gwasdb_res$Pos,sep='_'))

gwava_res <- read.table(opts$gwava,head=T,check.names=F,sep='\t')
gwava_res_by<- transform(gwava_res,merge1=paste(gwava_res$chr,gwava_res$end,sep='_'))

sce_res <- read.table(opts$sce,head=T,check.names=F,sep='\t')
sce_res_by<- transform(sce_res,merge1=paste(sce_res$chr,sce_res$end,sep='_'))

region_item=list(raw_res_by,gwasdb_res_by,gwava_res_by,sce_res_by)

raw_res_by <- Reduce(function(x,y) {merge(x,y,by ='merge1',all.x=TRUE)},region_item)
raw_res_by <- raw_res_by[!duplicated(raw_res_by["merge"]),]

eigen_res <- read.table(opts$eigen,head=T,check.names=F)
eigen_res_by <- transform(eigen_res,merge=paste(eigen_res$chr,eigen_res$position,eigen_res$ref,eigen_res$alt,sep='_'))
eigen_res_by<- eigen_res_by[!duplicated(eigen_res_by["merge"]),]

primateai_res <- read.table(opts$primateai,head=T,check.names=F)
primateai_res_by <- transform(primateai_res,merge=paste(primateai_res$chr,primateai_res$end,primateai_res$ref,primateai_res$alt,sep='_'))
primateai_res_by<- primateai_res_by[!duplicated(primateai_res_by["merge"]),]

spidex_res <- read.table(opts$spidex,head=T,check.names=F,sep='\t')
if(length(grep('chr',spidex_res$chromosome))>0){
  spidex_res$chromosome=gsub('chr','',spidex_res$chromosome)
 }
spidex_res_by<- transform(spidex_res,merge=paste(spidex_res$chromosome,spidex_res$position,spidex_res$ref_allele,spidex_res$mut_allele,sep='_'))
spidex_res_by <- spidex_res_by[!duplicated(spidex_res_by["merge"]),]

annovar_res <- read.csv(opts$annovar,head=T,check.names=F,sep=',')
annovar_res <- annovar_res[,-ncol(annovar_res)]
if(length(grep('chr',annovar_res$Chr))>0){
  annovar_res$Chr=gsub('chr','',annovar_res$Chr)
 }
annovar_res_by<- transform(annovar_res,merge=paste(annovar_res$Chr,annovar_res$End,annovar_res$Ref,annovar_res$Alt,sep='_'),check.names=F)
annovar_res_by <- annovar_res_by[!duplicated(annovar_res_by["merge"]),]

cadd_res <- read.table(opts$cadd,head=T,check.names=F,sep='\t')
cadd_res_sub <- subset(cadd_res,(AnnoType=='CodingTranscript' ))
cadd_res_sub <-cadd_res_sub[!duplicated(cadd_res_sub),]
cadd_res_by<- transform(cadd_res_sub,merge=paste(cadd_res_sub$'Chrom',cadd_res_sub$Pos,cadd_res_sub$Ref,cadd_res_sub$Alt,sep='_'),check.names=F)


silva_res <-read.table(opts$silva_res,head=T,sep='\t')
silva_res_by<- transform(silva_res,merge=paste(silva_res$'chrom',silva_res$pos,silva_res$ref,silva_res$alt,sep='_'))
silva_res_by <- silva_res_by[!duplicated(silva_res_by["merge"]),]

silva_mat <-read.table(opts$silva_mat,head=T,sep='\t')
silva_mat_by<- transform(silva_mat,merge=paste(silva_mat$'chrom',silva_mat$pos,silva_mat$ref,silva_mat$alt,sep='_'))
silva_mat_by <- silva_mat_by[!duplicated(silva_mat_by["merge"]),]

rna_res <-read.table(opts$rnascore,head=T)
rna_res_by<- transform(rna_res ,merge=paste(rna_res$'chr',rna_res$'pos',rna_res$'ref',rna_res$'alt',sep='_'))

dna_item <- list(raw_res_by,eigen_res_by,spidex_res_by,annovar_res_by,cadd_res_by,silva_res_by,silva_mat_by,primateai_res_by)
dna_res <- Reduce(function(x,y) {merge(x,y,by ='merge',all.x=TRUE)},dna_item)

write.table(dna_res,'dna_res.txt',row.names=F,quote=F,sep='\t')

all_res <- merge(dna_res,rna_res_by,by.x =c('merge','FeatureID'),by.y=c('merge','transcript_id'),all=TRUE)
all_res <- all_res[!duplicated(all_res),]

del_col=c(1,3,13,7,10,11,12,14,18,23,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,55,56,57,61,62,
63,64,65,66,73,74,75,94,95,96,97,100,102,107,108,109,110,111,112,113,114,115,124,125,126,127,128,129,130,131,132,133,134,
135,136,137,138,139,140,141,142,143,144,145,146,147,148,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,
168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,188,189,190,191,197,198,199,200,201,302,315,316,317,
318,319,320,321,322,323,324,325,326,327,328,329,332,355:365,367:371,373:375,382,387)

final_res <-all_res[-del_col]
final_res_colname <- names(final_res)
final_res_colname <-gsub('\\.y','',final_res_colname)
final_res_colname <-gsub('\\.x','',final_res_colname)
final_res_colname[179:181]=c('CADD_RawScore','CADD_PHRED','silva_rank')
final_res_colname[c(2:6,216:219)] =c('chr','start','end','ref','alt','remuRNA_MFE_wt','remuRNA_MFE_mu','remuRNA_dMFE','remuRNA_H')
final_res_colname[210:215] = c('Rnasnp_max_k','Rnasnp_max_d','Rnasnp_max_p','Rnasnp_interval','Rnasnp_d','Rnasnp_p')

colnames(final_res) = final_res_colname
write.table(final_res,'all_res.txt',row.names=F,quote=F,sep='\t')