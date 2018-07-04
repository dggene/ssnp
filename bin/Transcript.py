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
    codes = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    return codes[genotype]


class Transcript(object):
    transcriptId = None
    chRom = None
    chLoc = None
    ref = None
    alt = None
    cdsLoc = None
    geneId = None
    geneName = None
    seqRecord = None
    wtype = None
    mtype = None
    wtSeq = None
    mtSeq = None
    strand = None
    fastaDb = None

    def __init__(self, transcriptId, cdsLoc, ref, alt, fastaDb):
        self.transcriptId = transcriptId
        self.cdsLoc = int(cdsLoc)
        self.ref = ref
        self.alt = alt
        self.fastaDb = fastaDb
        self._init_data()

    def _init_data(self):
        self._init_seq_record()
        self._init_record_strand()
        self._init_wt_mt()

    def _init_seq_record(self):
        """
        根据fasta数据库获取bioseq的序列对象
        """
        for record in SeqIO.parse(self.fastaDb, 'fasta'):
            if self.transcriptId in record.id:
                self.seqRecord = record
                if 'N' in str(record.seq):
                    print('transcript_id:%s has invalid seq:%s'%(self.transcriptId,record.seq))
                    exit(0)
                return
        print('non transcript_id%s: found in fasta_db:%s' %
                        (self.transcriptId, self.fastaDb))
        exit(0)

    def _init_record_strand(self):
        """
        根据cds的fasta序列，判断正负链
        """
        ss = self.seqRecord.description.split(' ')[2]
        res = ss.split(':')[-1]
        self.strand = res

    def _init_wt_mt(self):
        """
        获取对应链上的基因型
        """
        ref = self.ref
        alt = self.alt
        if self.strand == '1':
            wtype = ref
            mtype = alt
        else:
            wtype = get_antisense_genotype(ref)
            mtype = get_antisense_genotype(alt)
        seq_list = list(self.seqRecord.seq)
        if seq_list[self.cdsLoc-1] != wtype:
            raise Exception('{transcript_id} with cdsLoc={cdsLoc} wtype={wtype} not equal to cds seq:{seq}'.format(
                transcript_id=self.transcriptId,
                cdsLoc=self.cdsLoc,
                wtype=wtype,
                seq=self.seqRecord.format('fasta')
            ))
        self.wtype = wtype
        self.mtype = mtype
        self.wtSeq = str(self.seqRecord.seq)
        seq_list = list(self.wtSeq)
        seq_list[self.cdsLoc-1] = mtype
        self.mtSeq = ''.join(seq_list)

    def cut_seq(self, seq, cdsLoc, length):
        '''
        获取指定序列按照长度切割的结果
        '''
        seq_list = list(seq)
        seq = ''
        new_loc = cdsLoc
        count = len(seq_list)
        if count > length:
            if cdsLoc < length//2:
                seq = seq_list[0:length]
                new_loc = cdsLoc
            elif (count-cdsLoc) < length//2:
                seq = seq_list[0-length:]
                new_loc = length-(count-cdsLoc)
            else:
                seq = seq_list[cdsLoc-length//2:cdsLoc+length//2]
                new_loc = length//2
        else:
            seq = seq_list
            new_loc = cdsLoc
        return (''.join(seq), new_loc)

    def render_seq_to_file(self, seq_id, seq, filename, footer=''):
        #seq=seq.replace('N','')
        with open(filename, 'w') as f:
            f.write(seq_id+'\n')
            f.write(seq+'\n')
            f.write(footer+'\n')

    def get_default_seq_id(self):
        wt_seq_id = ">{wtype}|{mtype}|{transcript_id}|{cdsLoc}|{cdsLen}".format(
            chrid=self.chRom,
            chrloc=self.chLoc,
            wtype=self.wtype,
            mtype=self.mtype,
            transcript_id=self.transcriptId,
            cdsLoc=self.cdsLoc,
            cdsLen=len(self.wtSeq))
        return wt_seq_id

    def generateNormalSeqFile(self, wt_filename='wt.fasta', mt_filename='mt.fasta', length=0):
        wt_seq_id = self.get_default_seq_id()+"|WT"
        mt_seq_id = self.get_default_seq_id()+"|MT"
        wt_seq = self.wtSeq
        mt_seq = self.mtSeq
        if length > 0:
            wt_seq = self.cut_seq(wt_seq, self.cdsLoc, length)
            mt_seq = self.cut_seq(mt_seq, self.cdsLoc, length)
        self.render_seq_to_file(wt_seq_id, wt_seq, wt_filename)
        self.render_seq_to_file(mt_seq_id, mt_seq, mt_filename)
        with open('transcript.id','w') as f:
            f.write('transcript_id\n')
            f.write(self.transcriptId+'\n')
        

    def generateRNASnpSeqFile(self, filename='seqss.txt'):
        """
        生成RNASnp分析用的文件
        """
        with open(filename, 'w') as snp_file:
            snp_file.write('%s%d%s' % (self.wtype, self.cdsLoc, self.mtype))

    def generateMenuRnaSeqFile(self, filename='remurna.seq'):
        seq_id = self.get_default_seq_id()
        new_seq, new_loc = self.cut_seq(self.wtSeq, self.cdsLoc, 150)
        footer = '*{wtype}{new_pos}{mtype}'.format(
            wtype=self.wtype, new_pos=new_loc, mtype=self.mtype)
        self.render_seq_to_file(seq_id, new_seq, filename, footer)

    def generateRnaFoldSeqFile(self, wt_filename='rnafold_wt.seq', mt_filename='rnafold_mt.seq'):
        cut_len = 175
        wt_new_seq, wt_new_loc = self.cut_seq(self.wtSeq, self.cdsLoc, cut_len)
        wt_seq_id = self.get_default_seq_id()+"|WT|"+str(wt_new_loc)
        self.render_seq_to_file(wt_seq_id, wt_new_seq, wt_filename)

        mt_new_seq, mt_new_loc = self.cut_seq(self.mtSeq, self.cdsLoc, cut_len)
        mt_seq_id = self.get_default_seq_id()+"|MT|"+str(mt_new_loc)
        self.render_seq_to_file(mt_seq_id, mt_new_seq, mt_filename)

def check_valid(transcript_id,cdsLoc):
    '''
    校验cadd的结果是否规范
    '''
    valid=True
    if transcript_id == 'NA':
        valid=False
    if cdsLoc == 'NA':
        valid=False
    return valid

def run(args):
    transcript_id=args.transcriptid
    cdsLoc=args.cdsloc
    ref=args.ref
    alt=args.alt
    if check_valid(transcript_id,cdsLoc):
        tran=Transcript(transcript_id,cdsLoc,ref,alt,args.fasta)
        tran.generateNormalSeqFile()
        tran.generateMenuRnaSeqFile()
        tran.generateRNASnpSeqFile()
        tran.generateRnaFoldSeqFile()
        return 
    else:
        #todo invalid
        pass

def main():
    ''' Main entry '''

    # _check_seafile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--ref',help='the cadd score file',required=True)
    parser.add_argument('-a','--alt',help='the cadd score file',required=True)
    parser.add_argument('-t','--transcriptid',help='the cadd score file',required=True)
    parser.add_argument('-p','--cdsloc',help='the cadd score file',required=True)
    parser.add_argument('-f','--fasta',help='the fasta file',required=True)
    parser.add_argument('-o','--output',help='the output file',required=True)

    if len(sys.argv) == 1:
        print(parser.format_help())
        return

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()