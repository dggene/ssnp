#!/usr/bin/env nextflow
//params.input="$baseDir/test/randomall.vcf"
params.input="$baseDir/test/pathogenic.vcf"
//params.input="$baseDir/test/test2.vcf"
params.output='pathogenic_output'
database="/DG/database/pub/ssnp"
params.silva_path="$database/silva"

params.gwasdb2_file="$database/gwasdb_20150819_annotation.gz"
params.spidex_file="$database/spidex_public_noncommercial_v1_0.tab.gz"
params.hg19genome ="$database/uchr1-x.fasta"

params.eigen_database="$database/eigen/eigen_noncoding"
params.annovar_database="/DG/database/genomes/Homo_sapiens/hg19/annovar"


Channel.fromPath(params.input).set { input_vcf0 }


process vcftocheck{
    conda="r-optparse bcftools bioconductor-biostrings"

    input:
        file 'input.vcf' from input_vcf0
    
    output:
        file 'adjust.vcf' into adjust_vcf
    
    script:
    """   
    echo \$PWD
    sed  "s/^chr//g" input.vcf > input_temp.vcf
    bcftools norm -c s input_temp.vcf -f ${params.hg19genome} 1> input_check.vcf 2> vcf_check.log
    Rscript $baseDir/bin/vcf_check.R -i input_check.vcf -l vcf_check.log -o adjust.vcf
    """
}

adjust_vcf.into{adjust_vcf0;adjust_vcf1;adjust_vcf2;adjust_vcf3}

process vcftobed{
    conda='bedops'

    input:
        file 'input.vcf' from adjust_vcf0
    
    output:
        file 'input.bed' into output_bed
    
    script:
    """   
    echo \$PWD
    convert2bed --input=VCF < input.vcf > input.bed
    """
}

output_bed.into{input0_bed;input1_bed;input2_bed;input3_bed}

process gwasdb2{
    conda="tabix"
    input:
        file 'input.bed' from input0_bed
    
    output:
        file 'gwasdb_res.tsv' into gwasdb_res

    script:
    """
    echo \$PWD
    zcat ${params.gwasdb2_file} |head -1 >head.txt
    tabix ${params.gwasdb2_file} -B input.bed > gwasdb_res_sub.tsv
    cat head.txt gwasdb_res_sub.tsv > gwasdb_res.tsv
    """
}

process spidex{
    conda="tabix"
    input:
        file 'input.bed' from input1_bed

    output:
        file 'spidex_res.tsv' into spidex_res

    script:
    """
        echo \$PWD
        zcat ${params.spidex_file} |head -1 > head.txt
        tabix ${params.spidex_file} -B input.bed > spidex_res_sub.tsv
        cat head.txt spidex_res_sub.tsv > spidex_res.tsv
    """
}


process eigen{
    conda="r-optparse tabix"

    input:
        file 'input.bed' from input2_bed
    
    output:
        file 'eigen_res.tsv' into eigen_res

    script:
    """
        echo \$PWD  
        Rscript $baseDir/bin/extract_Eigenscores.R -i input.bed -d ${params.eigen_database}
        cat ${params.eigen_database}/eigen_head.txt  eigen_score_* > eigen_res.tsv 
    """
}


process annovar{

    input:
        file 'input.vcf' from adjust_vcf1
    
    output:
        file 'out.hg19_multianno.csv' into annovar_res

    script:
    """
    echo \$PWD 
    perl $baseDir/bin/convert2annovar.pl -format vcf4 input.vcf -allsample -withfreq > input.avinput
    perl $baseDir/bin/table_annovar.pl --remove --otherinfo --csvout --buildver hg19 \
        -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,avsnp147,avsift,dbnsfp33a,caddgt20,clinvar_20170130 \
        -operation g,r,r,f,f,f,f,f,f,f,f \
        input.avinput \
        ${params.annovar_database} \
        -nastring NA \
        --outfile out
    """
}


process silva{
    conda='r-randomforest'
    input:
        file 'input.vcf' from adjust_vcf2
    
    output:
        file 'silva_result_final.tsv' into silva_res
    
    script:
    """
    export PATH=${params.silva_path}:\$PATH
    silva-preprocess ./tmp input.vcf
    silva-run ./tmp >silva_temp_result.tsv
    sed  -e '1d' -e 's/#//' -e 's/score/silva_score/' silva_temp_result.tsv |awk -F '\\t' '{for(i=1;i<11;i++)printf \$i"\\t";printf("\\n")}' > silva_result.tsv
    paste silva_result.tsv ./tmp/input.mat   >  silva_temp_result_final.tsv
    sed  -e 's/?//g'  -e 's/#//g' silva_temp_result_final.tsv > silva_result_final.tsv

    """

}


process cadd{
    conda="pysam"

    input:
        file('input.vcf') from adjust_vcf3
    output:
        file('res.score') into cadd_res
    script:
   
    """
    echo \$PWD
    sed 's/chr//g' input.vcf > trim_chr.vcf
    cat trim_chr.vcf|python $baseDir/bin/extractCADDscores.py -p $database/whole_genome_SNVs_inclAnno.tsv.gz>res_temp.score
    sed 's/#//g' res_temp.score > res.score
    
    """
}

cadd_res.into{cadd_res0;cadd_res1}
cadd_res0.splitCsv(header:true,sep:'\t').set{transcripts}


process getSeq{
    conda="biopython rnasnp viennarna r-optparse"
    input:
        val row from transcripts
        
    when:
        row.AnnoType=~'^CodingTranscript'

    output:
        file('paste.res') into paste_res

    script:
    def transcriptid=row.FeatureID
    
    """
    echo \$PWD  
   
    python $baseDir/bin/Transcript.py --ref ${row.Ref} --alt ${row.Alt} --transcriptid ${row.FeatureID} --cdsloc ${row.CDSpos} -f $database/Homo_sapiens.GRCh37.75.cds.all.fa -o .       
    if [ -f "wt.fasta" ]
    then
        echo "chr pos ref alt CDSpos" "\n"${row.Chrom} ${row.Pos} ${row.Ref} ${row.Alt} ${row.CDSpos} > persnp.txt
         
        RNAsnp -f wt.fasta -s seqss.txt -m 2 >rnasnp.res

        $baseDir/bin/remuRNA remurna.seq >remurna.res

        RNAfold --noPS <rnafold_wt.seq >rnafold_wt.res
        RNAfold --noPS <rnafold_mt.seq >rnafold_mt.res
        python $baseDir/bin/collect_rnafold.py -w rnafold_wt.res -m rnafold_mt.res -o .

        python $baseDir/bin/calc_hcu.py -w wt.fasta -m mt.fasta -f $database/codon/codon_frequency.txt -o hcu.res

        python $baseDir/bin/calc_rscu.py -w wt.fasta -m mt.fasta -o rscu.res

        perl $baseDir/bin/codonM wt.fasta wt.m
        perl $baseDir/bin/codonM mt.fasta mt.m
        Rscript $baseDir/bin/calc_tAi.R  -d $database -b $baseDir

        printf "GLOBAL_RATE\t0.06" > initRateFile
        cat wt.fasta mt.fasta >join.seq
        Rscript $baseDir/bin/Calc_rfm.R -t transcriptid.id -d $database -b $baseDir

        sed -e "/Warnings/d" -e "/^[[:space:]]*\$/d"  rnasnp.res > rnasnp_tmp.res
        paste persnp.res transcript.id rnasnp_tmp.res remurna.res rnafold.res tai.res hcu.res rscu.res rfm.res > paste.res

    fi    
    """

}


paste_res.collectFile(name:"${params.output}/trans_score.txt",keepHeader:true, sort: true, newLine: true).set{trans_res}


process merge_res{
    conda="r-optparse"

    publishDir {params.output}
    
    input:
        file('input.bed') from input3_bed
        file('trans.res') from trans_res
        file('gwasdb.res') from gwasdb_res
        file('spidex.res') from spidex_res
        file('eigen.res') from eigen_res
        file('annovar.res') from annovar_res
        file('cadd.score') from cadd_res1
        file('silva.score') from silva_res


    output:
        file('all_res.txt') into final_res

    """
    echo \$PWD
    Rscript $baseDir/bin/annotation_rbind.R -r input.bed -w gwasdb.res -e eigen.res -s spidex.res -a annovar.res -c cadd.score -i silva.score -t trans.res
    """
}


workflow.onComplete {
    def msg="""
Pipeline execution summary
---------------------------
ScriptId        :   ${workflow.scriptId}
ScriptName      :   ${workflow.scriptName}
scriptFile      :   ${workflow.scriptFile}
Repository      :   ${workflow.repository?:'-'}
Revision        :   ${workflow.revision?:'-'}
ProjectDir      :   ${workflow.projectDir}
LaunchDir       :   ${workflow.launchDir}
ConfigFiles     :   ${workflow.configFiles}
Container       :   ${workflow.container}
CommandLine     :   ${workflow.commandLine}
Profile         :   ${workflow.profile}
RunName         :   ${workflow.runName}
SessionId       :   ${workflow.sessionId}
Resume          :   ${workflow.resume}
Start           :   ${workflow.start}

Completed at    :   ${workflow.complete}
Duration        :   ${workflow.duration}
Success         :   ${workflow.success}
Exit status     :   ${workflow.exitStatus}
ErrorMessage    :   -
Error report    :   -
"""
    log.info(msg)
}