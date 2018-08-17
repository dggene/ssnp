# -*- coding: utf-8 -*-
#!/usr/bin/env python

import csv
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

def get_antisense_genotype(genotype):
    '''
    获取反义基因型
    '''
    codes={'A':'T','C':'G','T':'A','G':'C'}
    return codes[genotype]

def check_valid(row):
    '''
    校验cadd的结果是否规范
    '''
    valid=True
    if row["ConsDetail"] != "synonymous":
        valid=False
    if row['FeatureID'] == 'NA':
        valid=False
    if row['CDSpos'] == 'NA':
        valid=False
    return valid

def cut_seq(seq,cdsLoc,length):
    '''
    获取指定序列按照长度切割的结果
    '''
    seq_list=list(seq)
    seq=''
    new_loc=cdsLoc
    count=len(seq_list)
    if count>length:
        if cdsLoc<length//2:
            seq=seq_list[0:length]
            new_loc=cdsLoc
        elif (count-cdsLoc)<length//2:
            seq=seq_list[0-length:]
            new_loc=length-(count-cdsLoc)
        else:
            seq=seq_list[cdsLoc-length//2:cdsLoc+length//2]
            new_loc=length//2
    else:
        seq=seq_list
        new_loc=cdsLoc
    return (''.join(seq),new_loc)

def get_mt_seq(seq,cdsLoc,wtype,mtype):
    """
    获取突变的序列
    """
    seq_list=list(seq)
    seq_list[cdsLoc-1]=mtype
    mt_seq_seq=''.join(seq_list)
    return mt_seq_seq

def get_seq_record(transcript_id,fasta_db):
    """
    根据fasta数据库获取bioseq的序列对象
    """
    for record in SeqIO.parse(fasta_db,'fasta'):
        if transcript_id in record.id:
            return record
    raise Exception('non transcript_id%s: found in fasta_db:%s'%(transcript_id,fasta_db))

def get_record_strand(record):
    """
    根据cds的fasta序列，判断正负链
    """
    ss=record.description.split(' ')[2]
    res=ss.split(':')[-1]
    return res

def get_wt_mt(seq_record,chrid,chrloc,ref,alt,transcript_id,cdsLoc):
    """
    获取对应链上的基因型
    """
    strand=get_record_strand(record=seq_record)
    
    if strand=='1':
        wtype=ref
        mtype=alt
    else:
        wtype=get_antisense_genotype(ref)
        mtype=get_antisense_genotype(alt)
    seq_list=list(seq_record.seq)
    if seq_list[cdsLoc-1]!=wtype:
        raise Exception('{chrid}-{chrloc}:{ref}:{alt} with {transcript_id} not equal to cds seq:{seq}'.format(
            chrid=chrid,
            chrloc=chrloc,
            ref=ref,
            alt=alt,
            transcript_id=transcript_id,
            seq=seq_record.format('fasta')
        ))
    return (wtype,mtype)

def generate_seq_file(seq_record,chrid,chrloc,wtype,mtype,transcript_id,cdsLoc,gene_id,gene_name):
    """
    生成原始型和突变型两个序列文件
    """
    cdsLen=len(seq_record.seq)
    wt_seq_id="{chrid}|{chrloc}|WT|{wtype}|{mtype}|{transcript_id}|{cdsLoc}|{cdsLen}|{gene_id}|{gene_name}".format(
        chrid=chrid,
        chrloc=chrloc,
        wtype=wtype,
        mtype=mtype,
        transcript_id=transcript_id,
        cdsLoc=cdsLoc,
        cdsLen=cdsLen,
        gene_id=gene_id,
        gene_name=gene_name)
    wt_seq_obj=SeqRecord(seq_record.seq,wt_seq_id,'','')
    SeqIO.write(wt_seq_obj,'wt.fasta','fasta')

    mt_seq_id="{chrid}|{chrloc}|MT|{wtype}|{mtype}|{transcript_id}|{cdsLoc}|{cdsLen}|{gene_id}|{gene_name}".format(
        chrid=chrid,
        chrloc=chrloc,
        wtype=wtype,
        mtype=mtype,
        transcript_id=transcript_id,
        cdsLoc=cdsLoc,
        cdsLen=cdsLen,
        gene_id=gene_id,
        gene_name=gene_name)
    
    mt_seq_seq=get_mt_seq(seq_record.seq,cdsLoc,wtype,mtype)
    mt_seq_obj=SeqRecord(Seq(mt_seq_seq),mt_seq_id,'','')
    SeqIO.write(mt_seq_obj,'mt.fasta','fasta')

def generate_rnasnp_file(wtype,mtype,cdsLoc):
    with open('seqss.txt','w') as snp_file:
        snp_file.write('%s%d%s'%(wtype,cdsLoc,mtype))

def generate_menurna_file(seq_record,wtype,mtype,cdsLoc,transcript_id):
    new_seq,new_loc=cut_seq(seq_record.seq,cdsLoc,150)
    with open('remurna.seq','w') as remu_file:
        remu_file.write('>{transcript_id}|{wtype}|{mtype}|{cdsLoc}|{cdsLen}\n'.format(
            transcript_id=transcript_id,
            wtype=wtype,
            mtype=mtype,
            cdsLoc=cdsLoc,
            cdsLen=len(seq_record.seq)
        ))
        remu_file.write(new_seq+'\n')
        remu_file.write('*{wtype}{new_pos}{mtype}'.format(wtype=wtype,new_pos=new_loc,mtype=mtype))

def generate_rnafold_file(seq_record,wtype,mtype,cdsLoc,transcript_id):
    wt_new_seq,new_loc=cut_seq(seq_record.seq,cdsLoc,175)
    mt_seq=get_mt_seq(seq_record,cdsLoc,wtype,mtype)
    mt_new_seq=cut_seq(mt_seq,cdsLoc,175)
    with open('rnafold_wt.seq','w') as wt_seq_file:
        wt_seq_file.write('>{transcript_id}|wt|{cdsLoc}|{wtype}|{mtype}\n'.format(
            transcript_id=transcript_id,
            cdsLoc=new_loc,
            wtype=wtype,
            mtype=mtype))
        wt_seq_file.write(wt_new_seq)
    
    with open('rnafold_mt.seq','w') as mt_seq_file:
        mt_seq_file.write('>{transcript_id}|mt|{cdsLoc}|{wtype}|{mtype}\n'.format(
            transcript_id=transcript_id,
            cdsLoc=new_loc,
            wtype=wtype,
            mtype=mtype))
        mt_seq_file.write(mt_new_seq)
    


def run(args):
    with open(args.score_file) as  csvfile:
        reader=csv.DictReader(csvfile,dialect='excel-tab')
        for row in reader:
            print(row)
            if check_valid(row):
                chrid=row['#Chrom']
                chrloc=row['Pos']
                ref=row['Ref']
                alt=row['Alt']
                transcript_id=row['FeatureID']
                cdsLoc=int(row['CDSpos'])
                geneid=row['GeneID']
                genename=row['GeneName']
                seq_record=get_seq_record(transcript_id,args.fasta)
                wtype,mtype=get_wt_mt(seq_record,chrid,chrloc,ref,alt,transcript_id,cdsLoc)
                generate_seq_file(seq_record,chrid,chrloc,wtype,mtype,transcript_id,cdsLoc,geneid,genename)
                generate_rnasnp_file(wtype,mtype,cdsLoc)
                generate_menurna_file(seq_record,wtype,mtype,cdsLoc,transcript_id)
                generate_rnafold_file(seq_record,wtype,mtype,cdsLoc,transcript_id)
            else:
                #todo invalid
                pass

def main():
    ''' Main entry '''

    # _check_seafile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--score_file',help='the cadd score file',required=True)
    parser.add_argument('-f','--fasta',help='the fasta file',required=True)
    parser.add_argument('-o','--output',help='the output file',required=True)

    if len(sys.argv) == 1:
        print(parser.format_help())
        return

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()