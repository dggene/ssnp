#coding=utf-8
#!/usr/bin/env Rscript

options(warn = -1)
library("optparse")

option_list <- list(
make_option(c("-i", "--input"), type="character", default="input.txt",
              help="Input table file to read [default %default]"),
make_option(c("-d", "--database"), type="character", default="database",
              help="database directory  [default %default]") 
)            


opts = parse_args(OptionParser(option_list=option_list))

print(paste("The input file is ", opts$input,  sep = ""))
print(paste("The database dicectory is ", opts$database,  sep = ""))

input_bed <- read.delim(opts$input,head=F,check.names=F)[,1:3]
save.image('test.Rdata')
if( length(grep('chr',input_bed$V1)) > 0) input_bed$V1 =sub('chr','',input_bed$V1)
if( length(grep('X',input_bed$V1)) > 0) input_bed$V1 =sub('X','23',input_bed$V1)
if( length(grep('Y',input_bed$V1)) > 0) input_bed$V1 =sub('Y','24',input_bed$V1)

if(length(unique(input_bed$V1)) > 0) {
    run_temp <- lapply( as.numeric(as.character(unique(input_bed$V1))),function(nchr){
        print(nchr)
        if(nchr>=1 && nchr <=22) {
            input_bed_nchr <- subset(input_bed,subset = (V1==nchr))
            write.table(input_bed_nchr,paste('eigen_input_',nchr,'.bed',sep=''),quote=F,sep='\t',row.names=F,col.names=F)
            sh <- paste('tabix ',opts$database,'/Eigen_hg19_noncoding_annot_chr',nchr,'.tab.bgz',' -B ',
                    'eigen_input_',nchr,'.bed',' > ','eigen_score_',nchr,'.txt',sep='')
            print(sh)
            system(sh) }
        
        })    
  
    }else {
    print("The input file is NULL")
    }