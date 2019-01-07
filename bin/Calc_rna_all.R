#coding=utf-8
#!/usr/bin/env Rscript

library("optparse")

option_list <- list(

make_option(c("-d", "--database"), type="character", default="input.txt",
              help="database path [default %default]"),
make_option(c("-b", "--binpath"), type="character", default="input.txt",
              help="bin path [default %default]") )       


opts = parse_args(OptionParser(option_list=option_list))
save.image('test.Rdata')

if(all(file.exists(c('wt.fasta','transcript.id','seqss.txt','remurna.seq','rnafold_wt.seq')))){

    print('all files is ready')
    
    RNAsnp_sh=paste("Rscript ",opts$binpath,"/bin/RNAsnp.R -w wt.fasta -s seqss.txt",sep="")
    system(RNAsnp_sh)

    remuRNA_sh=paste(opts$binpath,"/bin/remuRNA remurna.seq >remurna.res",sep="")
    system(remuRNA_sh)

    system("RNAfold --noPS <rnafold_wt.seq >rnafold_wt.res")
    system("RNAfold --noPS <rnafold_mt.seq >rnafold_mt.res")
    rnafold_sh = paste("python ",opts$binpath,"/bin/collect_rnafold.py -w rnafold_wt.res -m rnafold_mt.res -o .",sep="")
    system(rnafold_sh)

    hcu_sh = paste ("python ",opts$binpath,"/bin/calc_hcu.py -w wt.fasta -m mt.fasta -f ",opts$database,"/codon/codon_frequency.txt -o hcu.res",sep="")
    system(hcu_sh)

    rscu_sh=paste("python ",opts$binpath,"/bin/calc_rscu.py -w wt.fasta -m mt.fasta -o rscu.res",sep="")
    system(rscu_sh)

    system(paste("perl ",opts$binpath,"/bin/codonM wt.fasta wt.m",sep=""))
    system(paste("perl ",opts$binpath,"/bin/codonM mt.fasta mt.m",sep=""))
    tai_sh=paste("Rscript ",opts$binpath,"/bin/calc_tAi.R -d ",opts$database," -b ",opts$binpath,sep="")
    system(tai_sh)

    system("printf 'GLOBAL_RATE\t0.06' > initRateFile")
    system("cat wt.fasta mt.fasta >join.seq")
    rfm_sh=paste("Rscript ",opts$binpath,"/bin/Calc_rfm.R -t transcript.id -d ",opts$database," -b ",opts$binpath,sep="")
    system(rfm_sh)

    system("sed -e '/Warnings/d' -e '/^[[:space:]]*$/d'  rnasnp.res > rnasnp_tmp.res")
    system("paste persnp.txt transcript.id rnasnp_tmp.res remurna.res rnafold.res tai_res.txt hcu.res rscu.res rfm.res > paste.res")

    quit(status=0)

}else{
    print("files is incomplete")
    quit(status=0)
}



