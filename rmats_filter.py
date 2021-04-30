#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd


def main(inf, outf, topnumber, p_value=0.05, fdr_value=None, inclevel=None):
    df = pd.read_csv(inf, header=0, delimiter="\t")
    df = df[df['PValue'] < p_value]
    if fdr_value:
        df = df[df['FDR'] < fdr_value]
    if inclevel:
        df = df[abs(df['IncLevelDifference']) > inclevel]

    df['sort_key'] = df.IncLevelDifference.abs()

    df.sort_values(by=['sort_key', 'FDR'], ascending=[False, True], inplace=True)
    df.drop('sort_key', axis=1, inplace=True)

    df = df[:topnumber]
    df.to_csv(outf, sep="\t", index=False)

    print(f"{df.shape[0]} lines.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input_file", type=str,
                        required=True, help="Input file gnerated by rmats")
    parser.add_argument("-o", "--output", dest="output_file",
                        type=str, help="Output file")
    parser.add_argument("-p", "--pvalue", dest="pvalue",
                        type=float, default=0.05, help="Pvalue cutoff")
    parser.add_argument("-f", "--fdr", dest="fdr",
                        type=float, help="FDR value cutoff")
    parser.add_argument("-l", "--inclevel", dest="inclevel",
                        type=float, help="IncLevel Difference cutoff")
    parser.add_argument("-n", "--number", dest="topnumber",
                        type=int, default=50, help="How many records to save")

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.topnumber,
         args.pvalue, args.fdr, args.inclevel)
