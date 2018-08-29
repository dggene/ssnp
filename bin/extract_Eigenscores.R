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
if( length(grep('MT',input_bed$V1)) > 0) input_bed$V1 =sub('MT','25',input_bed$V1)

if(length(unique(input_bed$V1)) > 0) {
    run_temp <- lapply( as.numeric(as.character(unique(input_bed$V1))),function(nchr){
        print(nchr)
        if(nchr>=1 && nchr <=22) {
            input_bed_nchr <- subset(input_bed,subset = (V1==nchr))
            print(nrow(input_bed_nchr))
            input_region_nchr <- apply (input_bed_nchr,1,function(x){
                            
                             trim_blank <- function(x){
                                 gsub(' +','',x)
                                 }                            
                             res <- paste(trim_blank(x[1]),':',trim_blank(x[2]),'-',trim_blank(x[3]),sep='')
                             return(res)                                                                                            
                            }) 

            region_all <- paste(input_region_nchr,collapse=' ')  
            sh <- paste('tabix ',opts$database,'/Eigen_hg19_noncoding_annot_chr',nchr,'.tab.bgz ',
                    region_all,' > ','eigen_score_',nchr,'.txt',sep='')     
            #print(sh)
            write.table(sh,paste(nchr,'_sh.txt',sep=''),quote=F,row.names=F,col.names=F)
            system(paste('sh ',nchr,'_sh.txt',sep=''),intern =TRUE) }
        
        })    
  
    }else {
    print("The input file is NULL")
    }