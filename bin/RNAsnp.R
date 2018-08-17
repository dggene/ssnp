#coding=utf-8
#!/usr/bin/env Rscript

library("optparse")

option_list <- list(

make_option(c("-w", "--wtfq"), type="character", default="input.txt",
              help="wt fasta to read [default %default]"),
make_option(c("-s", "--seqss"), type="character", default="input.txt",
              help="mt fasta to read [default %default]")                             
                                        
)

opts = parse_args(OptionParser(option_list=option_list))

sh=paste('RNAsnp -f ',opts$wtfq,' -s ',opts$seqss,' -m 2 >rnasnp.res',sep='')
system(sh)