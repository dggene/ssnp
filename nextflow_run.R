#coding=utf-8
#!/usr/bin/env Rscript

library("optparse")
library('vcfR')

option_list <- list(

make_option(c("-i", "--input"), type="character", default="input.txt",
              help="input.vcf file to read [default %default]"),
make_option(c("-r", "--resultdir"), type="character", default="input.txt",
              help="result directory to read [default %default]"),
make_option(c("-c", "--count"),type="integer", default=1,
              help="component [default %default]"),
make_option(c("-o", "--output"), type="character", default="input.txt",
              help="nextflowfile to output [default %default]")
)

opts = parse_args(OptionParser(option_list=option_list))

vcfSubmit <- function(vcf.head,vcf.fix,ind){
    
    n = nrow(vcf.fix) %/% ind
    print(paste("all set is split into ",n+1," group",sep=''))
    if(n == 0) {
        print(paste("run ",1," group",sep=''))
        file.name <- paste(inputdir,'/input_vcf/input_vcf.vcf',sep='')
        c <- file( file.name, "w" )
        writeLines( vcf.head, c )
        writeLines( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", c )
        vcf.fix.sub <- apply(vcf.fix,1,function(x){
                    paste(x,collapse='\t')
                })
        cat(vcf.fix.sub, c ,sep='\n')
        close( c )
        
        nextflow.sh <- paste('nextflow run main.nf --input=',file.name,' --output=',output,' -with-report',sep='')
        print(nextflow.sh)
        system(nextflow.sh)
        file.copy(paste('./',output,'/all_res.txt',sep=''),paste(inputdir,'/input_vcf/',output,'_res.txt',sep=''))
        
       }
 
    if(n > 0) {
        
        temp <- lapply(20:(n+1),function(x){
            print(x)
            print(paste("run ",x," group",sep=''))
            if( x <= n) {
                file.name <- paste(inputdir,'/input_vcf/input_vcf_',x,'.vcf',sep='')
                fileConn <- file( file.name, "wa" )
                writeLines( vcf.head, fileConn )
                writeLines( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", fileConn )
                print(dim(vcf.fix[(ind*(x-1)+1):(ind*x),]))   
                print((ind*(x-1)+1))
                print((ind*x))
                vcf.fix.sub = apply(vcf.fix[(ind*(x-1)+1):(ind*x),],1,function(x){
                        paste(x,collapse='\t')
                    })               
                writeLines(vcf.fix.sub,fileConn )
                close( fileConn )
                nextflow.sh <- paste('nextflow run main.nf --input=',file.name,' --output=',output,' -with-report',sep='')
                print(nextflow.sh)
                system(nextflow.sh)
                print(paste('./',output,'/all_res.txt',sep=''))
                print(paste(inputdir,'/out_result/',output,'_',x,'_res.txt',sep=''))
                file.copy(paste('./',output,'/all_res.txt',sep=''),paste(inputdir,'/out_result/',output,'_',x,'_res.txt',sep=''),overwrite = TRUE)
                
            }
            
            if( x > n) {
                file.name <- paste(inputdir,'/input_vcf/input_vcf_',x,'.vcf',sep='')
                c <- file( file.name, "w" )
                writeLines( vcf.head, c )
                writeLines( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", c )
                vcf.fix.sub <- apply(vcf.fix[(ind*(x-1)+1):nrow(vcf.fix),],1,function(x){
                        paste(x,collapse='\t')
                    })
                writeLines(vcf.fix.sub, c )
                close( c )
                
                nextflow.sh <- paste('nextflow run main.nf --input=',file.name,' --output=',output,' -with-report',sep='')
                print(nextflow.sh)
                system(nextflow.sh)
                print(paste('./',output,'/all_res.txt',sep=''))
                print(paste(inputdir,'/out_result/',output,'_',x,'_res.txt',sep=''))
                file.copy(paste('./',output,'/all_res.txt',sep=''),paste(inputdir,'/out_result/',output,'_',x,'_res.txt',sep=''),overwrite = TRUE)
            }                  
       })   

    }

  }

raw.vcf <- read.vcfR(opts$input)
#raw.vcf <- read.vcfR("/DG/home/zhum/sysdatabase/All_syn_1_snp150.vcf")

vcf.head <- raw.vcf@meta
vcf.fix <- raw.vcf@fix
args <- commandArgs(trailingOnly = F)
scriptPath <- dirname(sub("--file=","",args[grep("--file",args)]))
inputdir <- paste(scriptPath,opts$resultdir,sep='/')
output <- opts$output
ind <- opts$count

save.image('test.Rdata')
dir.create(paste(inputdir,'/input_vcf/',sep=''),recursive = TRUE,showWarnings=FALSE)
dir.create(paste(inputdir,'/out_result/',sep=''),recursive = TRUE,showWarnings=FALSE)  

print("input file is ok, the analysis is running")
res <- vcfSubmit(vcf.head,vcf.fix,ind)