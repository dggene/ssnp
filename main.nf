#!/usr/bin/env nextflow
params.input="$baseDir/test/test2.vcf"
params.output='output'
database="/DG/database/pub/ssnp"
params.silva_path="$database/silva"

params.gwasdb2_file="$database/gwasdb_20150819_annotation.gz"
params.spidex_file="$database/spidex_public_noncommercial_v1_0.tab.gz"
params.gwavas_file="$database/gwava_scores.bed.gz"
params.hg19genome ="$database/uchr1-x.fasta"
params.sce_database ="$database/SCE"

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

output_bed.into{input0_bed;input1_bed;input2_bed;input3_bed;input4_bed;input5_bed}

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
        awk '{ print "chr"\$1 "\\t" \$2 "\\t" \$3}' input.bed > chr_input.bed
        tabix ${params.spidex_file} -B chr_input.bed | sed  "s/^chr//g" > spidex_res_sub.tsv
        cat head.txt spidex_res_sub.tsv  > spidex_res.tsv
    """
}

process gwava{
    conda="tabix"
    input:
        file 'input.bed' from input4_bed

    output:
        file 'gwava_res.tsv' into gwava_res

    script:
    """
        echo \$PWD
        echo -e "chr\tstr\tend\trs\tgwava_s1\tgwava_s2\tgwava_s3" > head.txt
        awk '{ print "chr"\$1 "\t" \$2 "\t" \$3}' input.bed > chr_input.bed
        tabix ${params.gwavas_file} -B chr_input.bed |sed  "s/^chr//g"> gwava_res_sub.tsv
        cat head.txt gwava_res_sub.tsv  > gwava_res.tsv
    """
}


process sce{
    conda="bedtools"
    input:
        file 'input.bed' from input5_bed

    output:
        file 'sce_final_res.tsv' into sce_res

    script:
    """
        echo \$PWD
        echo -e "chr\tstr\tend\tsce9\tsce15\tsce30" > head.txt
        awk '{ print "chr"\$1 "\t" \$2 "\t" \$3  }' input.bed > chr_input.bed
        bedtools intersect -c -a chr_input.bed -b ${params.sce_database}9_hg19.bed > sce9.txt
        bedtools intersect -c -a chr_input.bed -b ${params.sce_database}15_hg19.bed > sce15.txt
        bedtools intersect -c -a chr_input.bed -b ${params.sce_database}30_hg19.bed > sce30.txt
        paste sce9.txt sce15.txt sce30.txt > sce_temp.txt
        awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$4 "\t" \$8 "\t" \$12}' sce_temp.txt | sed  "s/^chr//g" > sce_res.tsv
        cat head.txt sce_res.tsv >sce_final_res.tsv
   
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
    conda='r-randomforest numpy'
    input:
        file 'input.vcf' from adjust_vcf2
    
    output:
        file 'silva_result.tsv' into silva_res0
        file 'silva_mat_result.tsv' into silva_res1
    
    script:
    """
    export PATH=${params.silva_path}:\$PATH
    silva-preprocess ./tmp input.vcf
    silva-run ./tmp >silva_temp_result.tsv
    sed  -e '1d' -e 's/#//' -e 's/score/silva_score/' silva_temp_result.tsv |awk -F '\\t' '{for(i=1;i<11;i++)printf \$i"\\t";printf("\\n")}' > silva_result.tsv
    sed  's/#//' ./tmp/input.syn |awk -F '\\t' '{for(i=1;i<6;i++)printf \$i"\\t";printf("\\n")}' > input.syn
    paste input.syn ./tmp/input.mat   >  silva_temp_mat_result.tsv
    sed  -e 's/?//g'  -e 's/#//g' silva_temp_mat_result.tsv > silva_mat_result.tsv

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
    conda="biopython rnasnp viennarna r-optparse bioconductor-biocparallel"
    validExitStatus 0,1,2
    input:
        val row from transcripts
        
    when:
        row.AnnoType=~'^CodingTranscript'

    output:
        file('paste.res') optional true into paste_res

    script:
    def transcriptid=row.FeatureID
    
    """
    echo \$PWD  
   
    python $baseDir/bin/Transcript.py --ref ${row.Ref} --alt ${row.Alt} --transcriptid ${row.FeatureID} --cdsloc ${row.CDSpos} -f $database/Homo_sapiens.GRCh37.75.cds.all.fa -o .       
    if [ -f "wt.fasta" ]
    then
        echo -e "chr\tpos\tref\talt\tCDSpos\n"${row.Chrom}"\t"${row.Pos}"\t"${row.Ref}"\t"${row.Alt}"\t"${row.CDSpos} > persnp.txt       
    fi

    Rscript $baseDir/bin/Calc_RNAscore.R -d $database -b $baseDir
    """
}

paste_res.collectFile(name:"${params.output}/trans_score.txt",keepHeader:true).set{trans_res}

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
        file('silva.score') from silva_res0
        file('silva.mat') from silva_res1
        file('gwava.res') from gwava_res
        file('sce.res') from sce_res


    output:
        file('all_res.txt') into final_res

    """
    echo \$PWD
    Rscript $baseDir/bin/annotation_rbind.R -r input.bed -w gwasdb.res -v gwava.res -b sce.res -e eigen.res -s spidex.res -a annovar.res -c cadd.score -i silva.score -m silva.mat -t trans.res
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