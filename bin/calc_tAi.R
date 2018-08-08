#coding=utf-8
#!/usr/bin/env Rscript

library("optparse")

option_list <- list(

make_option(c("-d", "--database"), type="character", default="input.txt",
              help="database path to read [default %default]"),
make_option(c("-b", "--binpath"), type="character", default="input.txt",
              help="bin path to read [default %default]")                             
                                        
)

opts = parse_args(OptionParser(option_list=option_list))

source(paste(opts$binpath,'/bin/tAI.R',sep=''))
trna=scan(paste(opts$database,'/human.trna',sep=''))

ws=getws(tRNA = trna, sking=0)
w_m=matrix(scan("wt.m"), ncol=61,byrow=T)
w_tai_list=gettai(w_m, ws)

m_m=matrix(scan("mt.m"), ncol=61,byrow=T)
m_tai_list=gettai(m_m, ws)
save.image('text.Rdata')

len=nchar(readLines('mt.seq',skip=1)[2])


res=data.frame(
    WT_tAI=w_tai_list[2],
    MT_tAI=m_tai_list[2],
    tAI_rate=m_tai_list[2]/w_tai_list[2],
    tAI_rate_lg= (m_tai_list[2]/w_tai_list[2])**len 
)

write.table(res, "tai_res.txt", row.names = F,sep='\t', quote = F)






