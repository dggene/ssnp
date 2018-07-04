#coding=utf-8
import argparse
import sys
import re


def get_score(score_file):
    with open(score_file) as f:
        lines=f.readlines()
        last_line=lines[-1].strip().replace(' ','')
        res=re.findall(r"\([-+]?[0-9]*\.?[0-9]+\)",last_line)[0].replace('(','').replace(')','')
        return float(res)

def run(args):
    wt_score_file=args.wt_score
    mt_score_file=args.mt_score
    wt_score=get_score(wt_score_file)
    mt_score=get_score(mt_score_file)
    diff=mt_score-wt_score
    with open('rnafold.res','w') as f:
        f.write('rf_wt\trf_mt\trf_diff\n')
        f.write('%f\t%f\t%f'%(wt_score,mt_score,diff))
    

def main():
    ''' Main entry '''

    # _check_seafile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-w','--wt_score',help='the cadd score file',required=True)
    parser.add_argument('-m','--mt_score',help='the fasta file',required=True)
    parser.add_argument('-o','--output',help='the output file',required=True)

    if len(sys.argv) == 1:
        print(parser.format_help())
        return
    
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()