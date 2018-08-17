#coding=utf-8
#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
make_option(c("-r", "--rawbed"), type="character", default="input.txt",
              help="Input table file to read [default %default]"),
make_option(c("-w", "--gwasdb"), type="character", default="input.txt",
              help="gwasdb file to read [default %default]"),              
make_option(c("-e", "--eigen"), type="character", default="input.txt",
              help="Eigen file to read [default %default]"),
make_option(c("-s", "--spidex"), type="character", default="input.txt",
              help="spidex file to read [default %default]"),
make_option(c("-a", "--annovar"), type="character", default="input.txt",
              help="annovar to read [default %default]"),
make_option(c("-c", "--cadd"), type="character", default="input.txt",
              help="cadd file to read [default %default]"),
make_option(c("-i", "--silva"), type="character", default="input.txt",
              help="silva file to read [default %default]"),
make_option(c("-t", "--rnascore"), type="character", default="input.txt",
              help="rnascore file to read [default %default]")              
                            
)

opts = parse_args(OptionParser(option_list=option_list))

raw_res=read.table(opts$rawbed,sep='\t')[,c(1:4,6,7)]
if(length(grep('chr',raw_res$V1))>0){
    raw_res$V1=gsub('chr','',raw_res$V1)
    }
raw_res_by = transform(raw_res,merge=paste(raw_res$V1,raw_res$V3,raw_res$V6,raw_res$V7,sep='_'),merge1=paste(raw_res$V1,raw_res$V3,sep='_'))

save.image('test.Rdata')
gwasdb_res <- read.table(opts$gwasdb,head=T,check.names=F,sep='\t')
gwasdb_res_by<- transform(gwasdb_res,merge1=paste(gwasdb_res$Chr,gwasdb_res$Pos,sep='_'))

raw_res_by = merge(raw_res_by,gwasdb_res_by,by='merge1',all.x=TRUE)


eigen_res <- read.table(opts$eigen,head=T,check.names=F)
eigen_res_by <- transform(eigen_res,merge=paste(eigen_res$chr,eigen_res$position,eigen_res$ref,eigen_res$alt,sep='_'))

spidex_res <- read.table(opts$spidex,head=T,check.names=F,sep='\t')
if(length(grep('chr',spidex_res$chromosome))>0){
  spidex_res$chromosome=gsub('chr','',spidex_res$chromosome)
 }
spidex_res_by<- transform(spidex_res,merge=paste(spidex_res$chromosome,spidex_res$position,spidex_res$ref_allele,spidex_res$mut_allele,sep='_'))


annovar_res <- read.csv(opts$annovar,head=T,check.names=F,sep=',')
annovar_res <- annovar_res[,-ncol(annovar_res)]
if(length(grep('chr',annovar_res$Chr))>0){
  annovar_res$Chr=gsub('chr','',annovar_res$Chr)
 }
annovar_res_by<- transform(annovar_res,merge=paste(annovar_res$Chr,annovar_res$End,annovar_res$Ref,annovar_res$Alt,sep='_'))

cadd_res <- read.table(opts$cadd,head=T,check.names=F,sep='\t')
cadd_res_sub <- subset(cadd_res,(AnnoType=='CodingTranscript' ))
cadd_res_sub <-cadd_res_sub[!duplicated(cadd_res_sub),]
cadd_res_by<- transform(cadd_res_sub,merge=paste(cadd_res_sub$'Chrom',cadd_res_sub$Pos,cadd_res_sub$Ref,cadd_res_sub$Alt,sep='_'),check.names=F)

silva_res <-read.table(opts$silva,head=T,sep='\t')
silva_res_by<- transform(silva_res,merge=paste(silva_res$'chrom',silva_res$pos,silva_res$ref,silva_res$alt,sep='_'))

rna_res <-read.table(opts$rnascore,head=T)
rna_res_by<- transform(rna_res ,merge=paste(rna_res$'chr',rna_res$'pos',rna_res$'ref',rna_res$'alt',sep='_'))

dna_item <- list(raw_res_by,eigen_res_by,spidex_res_by,annovar_res_by,cadd_res_by,silva_res_by)
dna_res <- Reduce(function(x,y) {merge(x,y,by ='merge',all.x=TRUE)},dna_item)

write.table(dna_res,'dna_res.txt',row.names=F,quote=F,sep='\t')
rna_res <-read.table(opts$rnascore,head=T)

all_res <- merge(dna_res,rna_res_by,by =c('merge','CDSpos'),all=TRUE)
all_res <- dna_res[!duplicated(all_res),]

write.table(all_res,'all_res.txt',row.names=F,quote=F,sep='\t')