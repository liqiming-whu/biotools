#!/usr/bin/env python3
import re
import sys
import argparse
from itertools import groupby


def reader(fname, header=None, sep="\t", skip_while=None):
    sep = re.compile(sep)
    for line in open(fname):
        toks = sep.split(line.rstrip("\r\n"))
        if skip_while:
            if skip_while(toks):
                continue
        if header:
            yield header(toks)
        else:
            yield toks

class Bed6:
    """Convert GTF to Bed6."""
    __slots__ = ["chrom", "start", "end", "name", "score", "strand"]

    def __init__(self, gtf, attr):
        self.chrom = gtf.chrom
        self.start = gtf.start - 1
        self.end = gtf.end
        self.name = attr
        self.score = 0 if gtf.score == '.' else int(gtf.score)
        self.strand = gtf.strand

    def __repr__(self):
        return "Bed6({chrom}:{start}-{end})".format(chrom=self.chrom,
                start=self.start, end=self.end)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)


class Bed12:
    """Convert GTF to Bed12."""
    __slots__ = ["chrom", "start", "end", "name", "score", "strand",
                    "thickStart", "thickEnd", "itemRgb", "blockCount",
                    "blockSizes", "blockStarts"]
    def __init__(self, gtf, attr, thickStart, thickEnd):
        self.chrom = gtf.chrom
        self.start = gtf.start - 1
        self.end = gtf.end
        self.name = attr
        self.score = 0 if gtf.score == '.' else int(gtf.score)
        self.strand = gtf.strand
        self.thickStart = thickStart if thickStart else self.start
        self.thickEnd = thickEnd if thickEnd else self.start
        self.itemRgb = 0
        self.blockCount = 1
        self.blockSizes = "{size},".format(size=gtf.end - gtf.start + 1)
        self.blockStarts = "{start},".format(start=0)

    def __repr__(self):
        return "Bed12({chrom}:{start}-{end})".format(chrom=self.chrom,
                start=self.start, end=self.end)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)

    def add_exon(self, gtf):
        self.end = gtf.end
        self.blockSizes += "{size},".format(size=gtf.end - gtf.start + 1)
        self.blockStarts += "{start},".format(start=gtf.start - 1 - self.start)
        self.blockCount += 1


class Intron:
    """Convert GTF to Intron Bed6."""
    __slots__ = ["chrom", "start", "end", "name", "score", "strand"]

    def __init__(self, exon_1, exon_2, attr):
        self.chrom = exon_1.chrom
        self.start = exon_1.end
        self.end = exon_2.start - 1
        assert self.end >= self.start
        self.name = attr
        self.score = 0
        self.strand = exon_1.strand

    def __repr__(self):
        return "Intron({chrom}:{start}-{end})".format(chrom=self.chrom,
                start=self.start, end=self.end)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)


class GTF(object):
    __slots__ = ['seqname', 'source', 'feature', 'start', 'end', 'score',
                    'strand', 'frame', 'attributes', 'chrom']

    def __init__(self, args):
        for s, v in zip(self.__slots__[:9], args):
            setattr(self, s, v)
        self.start = int(self.start)
        self.end = int(self.end)
        self.chrom = self.seqname

    def __repr__(self):
        return "GTF({seqname}:{start}-{end})".format(seqname=self.seqname,
                start=self.start, end=self.end)


def main(gtf, attribute_id, feature, outformat):
    parse_attr = re.compile('{attribute_id} "([^"]+)"'.format(attribute_id=attribute_id))
    if outformat == 'bed12':
        cds = dict()
        for k, group in groupby(reader(gtf, header=GTF, skip_while = lambda toks: toks[0].startswith("#") or not (toks[2] == 'CDS' or toks[2] == 'stop_codon')), lambda x: parse_attr.findall(x.attributes)[0]):
            gtfs = list(sorted(group, key=lambda gtf: int(gtf.start)))
            thickStart = gtfs[0].start - 1
            thickEnd = gtfs[-1].end
            cds[k] = (thickStart, thickEnd)

    for k, group in groupby(reader(gtf, header=GTF, skip_while = lambda toks: toks[0].startswith("#") or not toks[2] == feature), lambda x: parse_attr.findall(x.attributes)[0]):
        # loop over the group building the single bed12 line
        assert outformat in ("bed12", "bed6", "intron"), "output file format: {} not support.".format(outformat)
        if outformat == 'bed12':
            bed = None
            for gtf_entry in sorted(group, key=lambda gtf: int(gtf.start)):
                assert gtf_entry.feature == feature
                if not bed:
                    try:
                        thickStart, thickEnd = cds[k]
                    except Exception:
                        thickStart = thickEnd = None
                    bed = Bed12(gtf_entry, k, thickStart, thickEnd)
                else:
                    bed.add_exon(gtf_entry)
            print(bed)
        if outformat == 'bed6':
            for gtf_entry in group:
                assert gtf_entry.feature == feature
                bed = Bed6(gtf_entry, k)
                print(bed)
        if outformat == 'intron':
            exons = []
            for gtf_entry in sorted(group, key=lambda gtf: int(gtf.start)):
                assert gtf_entry.feature == feature
                exons.append(gtf_entry)
            assert len(exons) > 0
            if len(exons) < 2:
                continue
            for i in range(len(exons) - 1):
                intron = Intron(exons[i], exons[i+1], k)
                print(intron)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('gtf', help="gtf file sorted and grouped by `attr`.")
    p.add_argument('attr', help="attribute ID by which to group bed12 entries.\
            example: gene_id or gene_name for gene.bed6, transcript_id for exon.bed12")
    p.add_argument('--feature-type', dest='feature', default='exon', help="feature type to \
            join -- all others are filtered out")
    p.add_argument('--outfile-format', dest='outformat', default='bed12', help="choose output file format: \
            bed12, bed6, intron")
    args = p.parse_args()
    main(args.gtf, args.attr, args.feature, args.outformat)
