#coding=utf-8
import argparse
import sys
import re
from Bio import SeqIO



STANDARD_CODON = {
    "TTT": 'F',
    "TTC": 'F',
    "TTA": 'L',
    "TTG": 'L',
    "TCT": 'S',
    "TCC": 'S',
    "TCA": 'S',
    "TCG": 'S',
    "TAT": 'Y',
    "TAC": 'Y',
    "TAA": 'B',
    "TAG": 'B',
    "TGT": 'C',
    "TGC": 'C',
    "TGA": 'B',
    "TGG": 'W',
    "CTT": 'L',
    "CTC": 'L',
    "CTA": 'L',
    "CTG": 'L',
    "CCT": 'P',
    "CCC": 'P',
    "CCA": 'P',
    "CCG": 'P',
    "CAT": 'H',
    "CAC": 'H',
    "CAA": 'Q',
    "CAG": 'Q',
    "CGT": 'R',
    "CGC": 'R',
    "CGA": 'R',
    "CGG": 'R',
    "ATT": 'I',
    "ATC": 'I',
    "ATA": 'I',
    "ATG": 'M',
    "ACT": 'T',
    "ACC": 'T',
    "ACA": 'T',
    "ACG": 'T',
    "AAT": 'N',
    "AAC": 'N',
    "AAA": 'K',
    "AAG": 'K',
    "AGT": 'S',
    "AGC": 'S',
    "AGA": 'R',
    "AGG": 'R',
    "GTT": 'V',
    "GTC": 'V',
    "GTA": 'V',
    "GTG": 'V',
    "GCT": 'A',
    "GCC": 'A',
    "GCA": 'A',
    "GCG": 'A',
    "GAT": 'D',
    "GAC": 'D',
    "GAA": 'E',
    "GAG": 'E',
    "GGT": 'G',
    "GGC": 'G',
    "GGA": 'G',
    "GGG": 'G'
}
CODON_NUM = {
    "TTT": 2,
    "TTC": 2,
    "TTA": 6,
    "TTG": 6,
    "TCT": 6,
    "TCC": 6,
    "TCA": 6,
    "TCG": 6,
    "TAT": 2,
    "TAC": 2,
    "TAA": 3,
    "TAG": 3,
    "TGT": 2,
    "TGC": 2,
    "TGA": 3,
    "TGG": 1,
    "CTT": 6,
    "CTC": 6,
    "CTA": 6,
    "CTG": 6,
    "CCT": 4,
    "CCC": 4,
    "CCA": 4,
    "CCG": 4,
    "CAT": 2,
    "CAC": 2,
    "CAA": 2,
    "CAG": 2,
    "CGT": 6,
    "CGC": 6,
    "CGA": 6,
    "CGG": 6,
    "ATT": 3,
    "ATC": 3,
    "ATA": 3,
    "ATG": 1,
    "ACT": 4,
    "ACC": 4,
    "ACA": 4,
    "ACG": 4,
    "AAT": 2,
    "AAC": 2,
    "AAA": 2,
    "AAG": 2,
    "AGT": 6,
    "AGC": 6,
    "AGA": 6,
    "AGG": 6,
    "GTT": 4,
    "GTC": 4,
    "GTA": 4,
    "GTG": 4,
    "GCT": 4,
    "GCC": 4,
    "GCA": 4,
    "GCG": 4,
    "GAT": 2,
    "GAC": 2,
    "GAA": 2,
    "GAG": 2,
    "GGT": 4,
    "GGC": 4,
    "GGA": 4,
    "GGG": 4
}


def calc_rscu(seq_record):
    (wtype, mtype, transcript_id, cdsLoc,
     cdsLen, ttype) = seq_record.description.split('|')
    codons = re.findall(r'\w{3}', str(seq_record.seq))
    count={}
    count_the=0
    count_syn=0

    #统计序列中每个密码子出现的个数
    for codon in codons:
        if codon in count.keys():
            count[codon]+=1
        else:
            count[codon]=1
    cdsLoc=int(cdsLoc)
    #获取当前突变位置的密码子
    if cdsLoc % 3 == 0:
        the_codon=codons[cdsLoc//3-1]
    else:
        the_codon=codons[cdsLoc//3]
    count_the=count[the_codon]

    for codon in codons:
        if STANDARD_CODON[codon]==STANDARD_CODON[the_codon]:
            count_syn+=1
    
    n=CODON_NUM[the_codon]
    print(the_codon)
    print(count)
    print('count_the:%d,count_syn:%d,n:%d'%(count_the,count_syn,n))
    rscu=n*count_the/count_syn
    
    return rscu

def run(args):
    wt_seq_obj = SeqIO.read(args.wt_seq, 'fasta')
    mt_seq_obj = SeqIO.read(args.mt_seq, 'fasta')
    wt_rscu = calc_rscu(wt_seq_obj)
    mt_rscu = calc_rscu(mt_seq_obj)
    d_hcu = mt_rscu-wt_rscu

    with open(args.output, 'w') as f:
        f.write('wt_rscu\tmt_rscu\td_rscu\n')
        f.write('%f\t%f\t%f' % (wt_rscu, mt_rscu, d_hcu))


def main():
    ''' Main entry '''

    # _check_seafile()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-w', '--wt_seq', help='the wtype seq file', required=True)
    parser.add_argument(
        '-m', '--mt_seq', help='the mtype seq file', required=True)
    parser.add_argument(
        '-o', '--output', help='the output file', required=True)

    if len(sys.argv) == 1:
        print(parser.format_help())
        return

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
