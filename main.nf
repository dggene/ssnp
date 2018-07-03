#!/usr/bin/env nextflow
params.input="$baseDir/test/test_input.vcf"
database="/DG/database/pub/ssnp"
params.silva_path="$database/silva"

params.gwasdb2_file="$database/gwasdb_20150819_annotation.gz"
params.spidex_file="$database/spidex_public_noncommercial_v1_0.tab.gz"

Channel.fromPath(params.input).into{input_vcf0;input_vcf1;input_vcf2}

/*
process silva{
    container="ssnp"

    input:
        file 'input.vcf' from input_vcf0
    
    output:
        file 'silva_result.tsv'
        file 'tmp/*.mat'
    """
    export PATH=${params.silva_path}:\$PATH
    silva-preprocess ./tmp input.vcf
    silva-run ./tmp >silva_result.tsv
    """
}
*/

input_vcf1.splitCsv(header:['chrom','pos','id','ref','alt','qual','filter','info'],skip:1,sep:'\t')
    .into{vcfs0;vcfs1}

/*
process gwasdb2{
    conda="/home/pyl/.conda/envs/ssnp_env"
    input:
        file 'input.vcf' from input_vcf1

    """
    echo \$PWD
    convert2bed --input=VCF < input.vcf > input.bed
    tabix ${params.gwasdb2_file} -R input.bed>gwasdb_res.tsv
    """
}

process spidex{
    conda="/home/pyl/.conda/envs/ssnp_env"
    input:
        file 'input.vcf' from input_vcf2

    """
        echo \$PWD
        convert2bed --input=VCF <input.vcf>input.bed
        tabix ${params.spidex_file} -R input.bed>spidex_res.tsv
    """
}
*/

process cadd{
    conda="pysam"
    input:
        file('input.vcf') from input_vcf2
    output:
        file('res.score') into cadd_score
    script:
   
    """
    echo \$PWD
    sed 's/chr//g' input.vcf > trim_chr.vcf
    cat trim_chr.vcf|python $baseDir/bin/extractCADDscores.py -p $database/whole_genome_SNVs_inclAnno.tsv.gz>res.score
    
    """
}
cadd_score.splitCsv(header:true,sep:'\t').set{transcripts}
process getSeq{
    conda="biopython"
    input:
        val row from transcripts

    when:
        row.ConsDetail=~'synonymous'

    output:
        set  file('wt.fasta'),file('mt.fasta') into seq_files
        file('transcript.id') into transcriptid_file
        file('seqss.txt') into rnasnp_input
        file('remurna.seq') into remurna_input
        set file('rnafold_wt.seq'),file('rnafold_mt.seq') into rnafold_input

    script:
    def transcriptid=row.FeatureID
    """
    echo \$PWD
    python $baseDir/bin/Transcript.py --ref ${row.Ref} --alt ${row.Alt} --transcriptid ${row.FeatureID} --cdsloc ${row.CDSpos} -f $database/Homo_sapiens.GRCh37.75.cds.all.fa -o .

    """
}

seq_files.into{seq_files0;seq_files1;seq_files2;seq_files3}
process rnasnp{
    conda="rnasnp"
    validExitStatus 0,160
    input:
        set file('wt.fasta'),file('mt.fasta') from seq_files0
        file('seqss.txt') from rnasnp_input
    output:
        file('rnasnp.res') into rnasnp_res
    script:
    """
    echo \$PWD
    RNAsnp -f wt.fasta -s seqss.txt >rnasnp.res
    """
}

process remuRNA{
    
    input:
       file('remurna.txt') from remurna_input
    output:
        file('remurna.res') into remurna_res
    script:
    """
    echo \$PWD
    $baseDir/bin/remuRNA remurna.txt >remurna.res
    """
}

process rnafold{
    conda='viennarna'
    input:
       set file('rnafold_wt.seq'),file('rnafold_mt.seq') from rnafold_input
    output:
        file('rnafold.res') into rnafold_res
    script:
    """
    echo \$PWD
    RNAfold --noPS < rnafold_wt.seq >rnafold_wt.res
    RNAfold --noPS <rnafold_mt.seq >rnafold_mt.res
    python $baseDir/bin/collect_rnafold.py -w rnafold_wt.res -m rnafold_mt.res -o .
    """        
}

process hcu{
    conda='biopython'
    input:
        set file('wt.seq'),file('mt.seq') from seq_files1
    output:
        file('hcu.res') into hcu_res
    script:
    """
    echo \$PWD
    python $baseDir/bin/calc_hcu.py -w wt.seq -m mt.seq -f $database/codon/codon_frequency.txt -o hcu.res
    """
}

process rscu{
    conda='biopython'
    input:
        set file('wt.seq'),file('mt.seq') from seq_files2
    output:
        file('rscu.res') into rscu_res
    
    script:
    """
    echo \$PWD
    python $baseDir/bin/calc_rscu.py -w wt.seq -m mt.seq -o rscu.res
    """
}

process rfm{
    input:
        set file('wt.seq'),file('mt.seq') from seq_files3
    output:
        file('rfm.res') into rfm_res
    script:
    """
    echo \$PWD
    printf "GLOBAL_RATE\t0.06" > initRateFile
    cat wt.seq mt.seq >join.seq
    java -jar $baseDir/bin/rfm/RFMapp.jar $database/codon/huCodonFile.txt join.seq 25 initRateFile . 0 0
    python $baseDir/bin/collect_rfm.py -i RFM_Result.txt -o rfm.res
    """
}

process paste_res{
    input:
        file('transcript.id') from transcriptid_file
        file('rnasnp.res') from rnasnp_res
        file('remurna.res') from remurna_res
        file('rnafold.res') from rnafold_res
        file('hcu.res') from hcu_res
        file('rscu.res') from rscu_res
        file('rfm.res') from rfm_res
    output:
        file('paste.res') into paste_res
    script:
    """
    echo \$PWD
    paste transcript.id rnasnp.res remurna.res rnafold.res hcu.res rscu.res rfm.res >paste.res
    """
}
paste_res.collectFile(name:'tran_score.txt',keepHeader:true).subscribe{println it}
