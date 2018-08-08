#coding=utf-8
#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
make_option(c("-t", "--transcriptid"), type="character", default="input.txt",
              help="transcriptid file to read [default %default]"),
make_option(c("-d", "--database"), type="character", default="input.txt",
              help="database path to read [default %default]"),
make_option(c("-b", "--binpath"), type="character", default="input.txt",
              help="bin path to read [default %default]")                             
                                        
)


opts = parse_args(OptionParser(option_list=option_list))

transcript_id <-read.table(opts$transcriptid,head=TRUE)
transcript_length <- read.table(paste(opts$database,'/rfm_database.txt',sep=''),head=TRUE)

save.image('text.Rdata')
if(transcript_id[1,1] %in% transcript_length[,1]) {
       rfm_res=transcript_length[which(transcript_length[,1] %in% transcript_id[1,1]),2] 
       final_res <- data.frame (
                    rfm_wt = rfm_res,
                    rfm_mt = rfm_res,
                    rfm_rate = 1.000000,
                    rfm_rate_lg = 1.000000
                    )
        write.table(final_res,'rfm.res',quote=F,sep='\t',row.names=F)
    } else {

    sh=paste('java -jar ',opts$binpath,'/bin/rfm/RFMapp.jar ',opts$database,'/codon/huCodonFile.txt join.seq 25 initRateFile . 0 0',sep='')
    system(sh)
    sh2= paste('python ',opts$binpath,'/bin/collect_rfm.py -i RFM_Result.txt -o rfm.res',sep='')
    system(sh2)
    }
