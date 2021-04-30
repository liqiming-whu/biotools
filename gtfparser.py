#!/usr/bin/env python3
import re
import argparse
from itertools import groupby


def reader(fname, header=None, sep="\t", skip_while=None):
    sep = re.compile(sep)
    with open(fname) as f:
        for line in f:
            toks = sep.split(line.rstrip("\r\n"))
            if skip_while:
                if skip_while(toks):
                    continue
            if header:
                yield header(toks)
            else:
                yield toks


def parse_attributes(term, attributes):
    parse_attr = re.compile('{term} "([^"]+)"'.format(term=term))
    patterns = parse_attr.findall(attributes)

    if patterns:
        return patterns[0]
    else:
        return 'NA'


class Alias:
    def __init__(self, data):
        self.data = data

    def update(self, data):
        self.data = data

    def __repr__(self):
        return str(self.data)

    __str__ = __repr__


class TransInfo:
    """Get transcript info"""
    __slots__ = ["transcript_id", "transcript_name", "transcript_type",
                 "gene_id", "gene_name", "gene_type", "strand"]

    def __init__(self, gtf, attr):
        self.strand = gtf.strand
        attributes = gtf.attributes
        self.transcript_id = attr
        self.transcript_name = parse_attributes("transcript_name", attributes)
        self.transcript_type = parse_attributes("transcript_type", attributes)
        self.gene_id = parse_attributes("gene_id", attributes)
        self.gene_name = parse_attributes("gene_name", attributes)
        self.gene_type = parse_attributes("gene_type", attributes)

    def __repr__(self):
        return "Transcript({transcript_id} {transcript_type} {strand})".format(
            transcript_id=self.transcript_id,
            transcript_type=self.transcript_type,
            strand=self.strand
            )

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)


class GeneInfo:
    """Get gene info"""
    __slots__ = ["gene_id", "gene_name", "gene_type", "strand"]

    def __init__(self, gtf, attr):
        self.strand = gtf.strand
        attributes = gtf.attributes
        self.gene_id = attr
        self.gene_name = parse_attributes("gene_name", attributes)
        self.gene_type = parse_attributes("gene_type", attributes)

    def __repr__(self):
        return "Gene({gene_id} {gene_type} {strand})".format(
            gene_id=self.gene_id,
            gene_type=self.gene_type,
            strand=self.strand
            )

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)


class Bed6:
    """Convert GTF to Bed6."""
    __slots__ = ["chrom", "start", "end", "name", "score", "strand"]

    def __init__(self, gtf, attr):
        self.chrom = gtf.chrom
        self.start = gtf.start - 1
        self.end = gtf.end
        self.name = attr
        self.score = gtf.score
        self.strand = gtf.strand

    def __repr__(self):
        return "Bed6({chrom}:{start}-{end})".format(
            chrom=self.chrom, start=self.start, end=self.end)

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
        self.end = Alias(gtf.end)
        self.name = attr
        self.score = 0 if gtf.score == '.' else int(gtf.score)
        self.strand = gtf.strand
        if thickStart and thickEnd:
            self.thickStart = thickStart
            self.thickEnd = thickEnd
        else:
            self.thickStart = self.end
            self.thickEnd = self.end
        self.itemRgb = 0
        self.blockCount = 1
        self.blockSizes = "{size},".format(size=gtf.end - gtf.start + 1)
        self.blockStarts = "{start},".format(start=0)

    def __repr__(self):
        return "Bed12({chrom}:{start}-{end})".format(
            chrom=self.chrom, start=self.start, end=self.end)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)

    def add_exon(self, gtf):
        self.end.update(gtf.end)
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
        self.score = '.'
        self.strand = exon_1.strand

    def __repr__(self):
        return "Intron({chrom}:{start}-{end})".format(
            chrom=self.chrom, start=self.start, end=self.end)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)


class GTF:
    __slots__ = ['seqname', 'source', 'feature', 'start', 'end', 'score',
                 'strand', 'frame', 'attributes', 'chrom']

    def __init__(self, args):
        for s, v in zip(self.__slots__[:9], args):
            setattr(self, s, v)
        self.start = int(self.start)
        self.end = int(self.end)
        self.chrom = self.seqname

    def __repr__(self):
        return "GTF({seqname}:{start}-{end})".format(
            seqname=self.seqname, start=self.start, end=self.end)


def main(gtf, attribute_id, feature, outformat):
    assert outformat in (
        "bed12", "bed6", "intron", "info"
        ), "output file format: {} not support.".format(outformat)
    parse_attr = re.compile(
        '{attribute_id} "([^"]+)"'.format(attribute_id=attribute_id))
    if outformat == 'bed12':
        cds = dict()
        for k, group in groupby(
            reader(
                gtf,
                header=GTF,
                skip_while=lambda toks: toks[0].startswith("#") or not (
                    toks[2] == 'CDS' or toks[2] == 'stop_codon')
                ),
            lambda x: parse_attr.findall(x.attributes)[0]
        ):
            gtfs = list(sorted(group, key=lambda gtf: gtf.start))
            thickStart = gtfs[0].start - 1
            thickEnd = gtfs[-1].end
            cds[k] = (thickStart, thickEnd)
    if outformat == 'info':
        if attribute_id == 'transcript_id':
            head = ["transcript_id", "transcript_name", "transcript_type",
                    "gene_id", "gene_name", "gene_type", "strand"]
            print("\t".join(head))
        elif attribute_id == 'gene_id':
            head = ["gene_id", "gene_name", "gene_type", "strand"]
            print("\t".join(head))
        else:
            raise Exception("GTF attribute and outformat not support.")

    for k, group in groupby(
        reader(
            gtf,
            header=GTF,
            skip_while=lambda toks: toks[0].startswith("#") or not (
                toks[2] == feature)
            ),
        lambda x: parse_attr.findall(x.attributes)[0]
    ):
        # loop over the group building the single bed12 line
        if outformat == 'bed12':
            bed = None
            for gtf_entry in sorted(group, key=lambda gtf: gtf.start):
                if not bed:
                    try:
                        thickStart, thickEnd = cds[k]
                    except Exception:
                        thickStart = thickEnd = None
                    bed = Bed12(gtf_entry, k, thickStart, thickEnd)
                else:
                    bed.add_exon(gtf_entry)
            print(bed)
        elif outformat == 'bed6':
            for gtf_entry in group:
                bed = Bed6(gtf_entry, k)
                print(bed)
        elif outformat == 'intron':
            exons = []
            for gtf_entry in sorted(group, key=lambda gtf: gtf.start):
                exons.append(gtf_entry)
            assert len(exons) > 0
            if len(exons) < 2:
                continue
            for i in range(len(exons) - 1):
                intron = Intron(exons[i], exons[i+1], k)
                print(intron)

        elif outformat == 'info' and attribute_id == 'transcript_id':
            for gtf_entry in group:
                transcript_info = TransInfo(gtf_entry, k)
                print(transcript_info)

        elif outformat == 'info' and attribute_id == 'gene_id':
            for gtf_entry in group:
                gene_info = GeneInfo(gtf_entry, k)
                print(gene_info)


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('gtf', help="gtf file download from Genecode.")
    p.add_argument('--attribute', dest='attr', default='transcript_id', help="""
    attribute ID by which to group bed entries.
    'gene_id' or 'gene_name' for gene.bed6, gene_info
    'transcript_id' for genome.bed12, exon, intron, transcrpt_info
    """)
    p.add_argument('--type', dest='feature', default='exon', help="""
    annotation type to join, all others are filtered out:
    'exon' genome.bed12, exon, intron, 'gene' for gene.bed6, gene_info,
    'transcript' for transcript_info.
    """)
    p.add_argument('--format', dest='outformat', default='bed12', help="""
    choose output file format:bed12, bed6, intron, info
    """)
    args = p.parse_args()
    main(args.gtf, args.attr, args.feature, args.outformat)
