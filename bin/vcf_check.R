library("optparse")

option_list <- list(
make_option(c("-i", "--input"), type="character", default="input.vcf",
              help="Input vcf file to read [default %default]"),
make_option(c("-l", "--log"), type="character", default="check.log",
              help="check log file to read [default %default]"),      
make_option(c("-o", "--output"), type="character", default="output.vcf",
              help="Output vcf file to read [default %default]")                                  
)

opts = parse_args(OptionParser(option_list=option_list))


vcf_error<-readLines(opts$log)
log <- strsplit(vcf_error[length(vcf_error)],'\t')[[1]][2]
modified <- strsplit(log,'/')[[1]][2]

if(as.numeric(modified)>0) {

library(Biostrings)

save.image('test.Rdata')
input_vcf<-readLines(opts$input)
vcf_head <- input_vcf[1:grep('#CHROM',input_vcf)]
vcf_info <- input_vcf[-(1:grep('#CHROM',input_vcf))]
vcf_tab = lapply(vcf_info,function(x){
    
        vcf_line <- strsplit(x,'\t')[[1]]
        ref <- vcf_line[4]
        alt <- vcf_line[5]
        if(nchar(alt)>1){
              allele <- strsplit(alt,',')[[1]]
              if (as.character(complement(DNAStringSet(allele[2]))) == ref) {
                  vcf_line[5] = as.character(complement(DNAStringSet(allele[1])))
                   }                     
             }
        if(alt=='.')print('Alt allele is error')        
        return(vcf_line)
        }) 
        
vcf_tab<-data.frame(do.call(rbind,vcf_tab))

vcf_tab_m=apply(vcf_tab,1,function(x){
                  paste(x,collapse='\t')})
                  
vcf_all <- data.frame(V1=c(vcf_head,vcf_tab_m))  
write.table(vcf_all[,1],'temp.vcf',quote=F,row.names=F,col.names=F)

sh=paste('bcftools norm -m- temp.vcf  > ',opts$output)
system(sh)
} else {
sh=paste('bcftools norm -m- ',opts$input, ' > ',opts$output)
system(sh)
}
