#!/usr/bin/env nextflow
params.input="$baseDir/test/test_input.vcf"
database="/DG/database/pub/ssnp"
params.silva_path="$database/silva"

params.gwasdb2_file="$database/gwasdb_20150819_annotation.gz"
params.spidex_file="$database/spidex_public_noncommercial_v1_0.tab.gz"

Channel.fromPath(params.input).into{input_vcf0;input_vcf1;input_vcf2}

process silva{
    container="ssnp"

    input:
        file 'input.vcf' from input_vcf0
    
    output:
        file 'silva_result.tsv'
        file 'tmp/*.mat'
    """
    which python
    export PATH=${params.silva_path}:\$PATH
    silva-preprocess ./tmp input.vcf
    silva-run ./tmp >silva_result.tsv
    """

}
process gwasdb2{
    conda="/home/pyl/.conda/envs/ssnp_env"
    input:
        file 'input.vcf' from input_vcf1

    """
    convert2bed --input=VCF < input.vcf > input.bed
    tabix ${params.gwasdb2_file} -R input.bed>result.tsv
    """
}

process spidex{
    conda="/home/pyl/.conda/envs/ssnp_env"
    input:
        file 'input.vcf' from input_vcf2

    """
        convert2bed --input=VCF <input.vcf>input.bed
        tabix ${params.spidex_file} -R input.bed>result.tsv
    """
}