#!/usr/bin/env python3
"""
Remove duplicated reads for paired-end sequencing.
Author: liqiming@whu.edu.cn
"""
import os
import gzip
import argparse


def main(read1, read2, out_read1, out_read2, dup_read1, dup_read2):
    unique = set()
    input_reads = 0
    output_reads = 0
    id1 = read1.readline()
    while id1:
        id2 = read2.readline()
        name1 = id1.decode().split()[0]
        name2 = id2.decode().split()[0]
        assert name1 == name2
        seq1 = read1.readline()
        seq2 = read2.readline()
        symbol1 = read1.readline()
        symbol2 = read2.readline()
        qual1 = read1.readline()
        qual2 = read2.readline()
        input_reads += 1

        whole_reads = seq1 + seq2
        if whole_reads in unique:
            dup_read1.write(id1+seq1+symbol1+qual1)
            dup_read2.write(id2+seq2+symbol2+qual2)
        else:
            out_read1.write(id1+seq1+symbol1+qual1)
            out_read2.write(id2+seq2+symbol2+qual2)
            unique.add(whole_reads)
            output_reads += 1
        id1 = read1.readline()
    out_read1.close()
    out_read2.close()
    dup_read1.close()
    dup_read2.close()

    duplicate_ratio = (input_reads - output_reads) * 100 / input_reads
    print(
        f"input read pairs: {input_reads}\noutput read pairs: {output_reads}\nduplicate ratio: {duplicate_ratio:.2f}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r1", dest="input1",
                        type=str, required=True, help="input read1.fq.gz")
    parser.add_argument("-r2", dest="input2",
                        type=str, required=True, help="input read2.fq.gz")
    parser.add_argument("-o1", dest="output1",
                        type=str, required=True, help="output read1.fq.gz")
    parser.add_argument("-o2", dest="output2",
                        type=str, required=True, help="output read2.fq.gz")
    parser.add_argument("-d1", dest="duplicate1", type=str,
                        default=os.devnull, help="output read1.dup.fq.gz")
    parser.add_argument("-d2", dest="duplicate2", type=str,
                        default=os.devnull, help="output read2.dup.fq.gz")
    args = parser.parse_args()

    print(f"processing {args.input1} {args.input2} ...")

    main(gzip.open(args.input1, "rb"), gzip.open(args.input2, "rb"),
         gzip.open(args.output1, "wb"), gzip.open(args.output2, "wb"),
         gzip.open(args.duplicate1, "wb"), gzip.open(args.duplicate2, "wb"))
