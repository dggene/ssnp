import argparse
import sys
import re


def run(args):
    wt_rfm=None
    mt_rfm=None
    cdsLen=None
    with open(args.input) as f:
        lines=f.readlines()
        for i in range(len(lines)):
            line=lines[i]
            if '>' not in line:
                continue
            rfm=float(lines[i+1].replace('Translation Rate = ',''))
            cdsLen=int(line.split('|')[4])
            if 'MT' in line:
                mt_rfm=rfm
            else:
                wt_rfm=rfm

    rfm_rate=mt_rfm/wt_rfm
    rfm_rate_lg=rfm_rate**cdsLen

    with open(args.output,'w') as f:
        f.write('rfm_wt\trfm_mt\trfm_rate\trfm_rate_lg\n')
        f.write('%f\t%f\t%f\t%f'%(wt_rfm,mt_rfm,rfm_rate,rfm_rate_lg))

def main():
    ''' Main entry '''

    # _check_seafile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help='the RFM_Result.txt',required=True)
    parser.add_argument('-o','--output',help='the output file',required=True)

    if len(sys.argv) == 1:
        print(parser.format_help())
        return
    
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()