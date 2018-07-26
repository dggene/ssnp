# coding=utf-8
import argparse
import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def cut_seq(seq, cdsLoc, length):
    '''
    获取指定序列按照长度切割的结果
    '''
    seq_list = list(seq)
    seq = ''
    new_loc = cdsLoc
    count = len(seq_list)
    if count > length:
        if cdsLoc <= length//2:
            seq = seq_list[0:length]
            new_loc = cdsLoc
        elif (count-cdsLoc) <= length//2:
            seq = seq_list[0-length:]
            new_loc = length-(count-cdsLoc)
        else:
            left = cdsLoc-length//2-1
            seq = seq_list[left:left+length]
            new_loc = length//2+1
    else:
        seq = seq_list
        new_loc = cdsLoc
    return (''.join(seq), new_loc)


def run(args):
    if not os.path.exists(args.input):
        sys.stderr.write('miss input file:%s\n' % (args.input))
        exit(127)
    record = SeqIO.read(args.input, 'fasta')
    items = record.description.split('|')
    if args.cdsPos:
        cdsLoc=args.cdsPos
    else:
        if len(items) != 5:
            sys.stderr.write(
                'input cdsPos or fasta id must as Wtype|Mtype|TranscriptId|cdsLoc|cdsLen')
            exit(128)
    cds_loc = items[3]
    (new_seq, new_loc) = cut_seq(record.seq, cds_loc, args.length)
    


def main():
    ''' Main entry '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', help='the input seq file', required=True)
    parser.add_argument(
        '-l', '--length', help='the length to cut', required=True)
    parser.add_argument('-p', '--cdsPos', help='cds position', required=False)
    parser.add_argument(
        '-o', '--output', help='the output file', required=True)
    parser.add_argument('-h', '--help', help='the help', required=False)

    if len(sys.argv) == 1:
        print(parser.format_help())
        return

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
