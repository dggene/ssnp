source("/DG/home/gwj/personal-software/sSNP/pipeline/codonR/tAI.R")
trna=scan("/DG/home/gwj/personal-software/sSNP/pipeline/codonR/human.trna")
ws=getws(tRNA = trna, sking=0)
m=matrix(scan("/tmp/sSNP/tAI/tmp.m"), ncol=61,byrow=T)
tai_list=gettai(m, ws)
write.table(tai_list[2], file = "/tmp/sSNP/tAI/tai.tmp", row.names = F, col.names = F, quote = F)
q()





