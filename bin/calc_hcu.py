import argparse
import sys
import re
from Bio import SeqIO


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


def get_frq(the_codon, cf_file):
    with open(cf_file) as f:
        lines = f.readlines()
        for line in lines:
            items = line.split('\t')
            if the_codon == items[0]:
                return float(items[2])
    raise Exception('the codon:%s not found in cf_file:%s' %
                    (the_codon, cf_file))


def calc_hcu(seq_record, cf_file):
    (wtype, mtype, transcript_id, cdsLoc,
     cdsLen,ttype) = seq_record.description.split('|')
    codons = re.findall(r'\w{3}', str(seq_record.seq))
    cdsLoc=int(cdsLoc)
    the_codon = codons[cdsLoc//3]
    n = CODON_NUM[the_codon]
    frq = get_frq(the_codon, cf_file)
    hcu = n*frq
    return hcu


def run(args):
    cf_file = args.codon_frequency
    wt_seq_obj = SeqIO.read(args.wt_seq, 'fasta')
    mt_seq_obj = SeqIO.read(args.mt_seq, 'fasta')
    wt_hcu = calc_hcu(wt_seq_obj, cf_file)
    mt_hcu = calc_hcu(mt_seq_obj, cf_file)
    d_hcu = mt_hcu-wt_hcu
    
    with open(args.output,'w') as f:
        f.write('wt_hcu\tmt_hcu\td_hcu\n')
        f.write('%f\t%f\t%f'%(wt_hcu,mt_hcu,d_hcu))

def main():
    ''' Main entry '''

    # _check_seafile()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-w', '--wt_seq', help='the wtype seq file', required=True)
    parser.add_argument('-m', '--mt_seq', help='the mtype seq file', required=True)
    parser.add_argument('-f', '--codon_frequency', help='codon_frequency file')
    parser.add_argument(
        '-o', '--output', help='the output file', required=True)

    if len(sys.argv) == 1:
        print(parser.format_help())
        return

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
