#!/usr/bin/env python3
import argparse
import pandas as pd


def main(inf, outf, topnumber, p_value=None, readnumber=None, fdr_value=None, inclevel=None):
    df = pd.read_csv(inf, header=0, delimiter="\t")
    if p_value:
        df = df[df['PValue'] < p_value]
    if fdr_value:
        df = df[df['FDR'] < fdr_value]
    if inclevel:
        df = df[abs(df['IncLevelDifference']) > inclevel]
    if readnumber:
        df['ReadNumber'] = df['IJC_SAMPLE_1'] + df['SJC_SAMPLE_1'] + df['IJC_SAMPLE_2'] + df['SJC_SAMPLE_2']
        df = df[df['ReadNumber'] > readnumber]
        df.drop('ReadNumber', axis=1, inplace=True)

    df['sort_key'] = df.IncLevelDifference.abs()

    df.sort_values(by=['sort_key', 'FDR'], ascending=[False, True], inplace=True)
    df.drop('sort_key', axis=1, inplace=True)

    if topnumber:
        df = df[:topnumber]
    if outf:
        df.to_csv(outf, sep="\t", index=False)

    print(f"{outf}: {df.shape[0]} genes.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input_file", type=str,
                        required=True, help="Input file gnerated by rmats")
    parser.add_argument("-o", "--output", dest="output_file",
                        type=str, default=None, help="Output file")
    parser.add_argument("-p", "--pvalue", dest="pvalue",
                        type=float, default=None, help="Pvalue cutoff")
    parser.add_argument("-f", "--fdr", dest="fdr",
                        type=float, default=None, help="FDR value cutoff")
    parser.add_argument("-l", "--inclevel", dest="inclevel",
                        type=float, default=0.01, help="IncLevel Difference cutoff")
    parser.add_argument("-r", "--readnumber", dest="readnumber",
                        type=int, default=None, help="Reads number cutoff")
    parser.add_argument("-n", "--number", dest="topnumber",
                        type=int, default=None, help="How many records to save")

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.topnumber, args.pvalue,
         args.readnumber, args.fdr, args.inclevel)
