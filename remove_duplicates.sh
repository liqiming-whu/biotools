#!/bin/sh
if [ $# != 5 ];then
    echo "Usage: ./remove_duplicates.sh threads read1.fq.gz read2.fq.gz out1.fq.gz out2.fq.gz"
    exit 1;
fi
# Accept parameters.
threads=$1
read1=$2
read2=$3
out1=$4
out2=$5

# Define temp filles.
prefix=$out1
tempdir=${prefix}.temp
decompressed1=${prefix}.temp/R1.decompressed.fastq
decompressed2=${prefix}.temp/R2.decompressed.fastq
uncompress1=${prefix}.temp/R1.undecompress.fastq
uncompress2=${prefix}.temp/R2.undecompress.fastq

# Run.
mkdir -p $tempdir
pigz -p $threads -d -c $read1 > $decompressed1
pigz -p $threads -d -c $read2 > $decompressed2
paste $decompressed1 $decompressed2 | \
    awk '{if(NR%4==0){printf $0"\n"}else{printf $0"\t"}}' | \
    sort -T $tempdir -t $'\t' -k3,3 -u | \
    awk -v FS='\t' '{print $1"\n"$3"\n"$5"\n"$7 > "'$uncompress1'"; print $2"\n"$4"\n"$6"\n"$8 > "'$uncompress2'"}'
pigz -p $threads -c $uncompress1 > $out1
pigz -p $threads -c $uncompress2 > $out2

# Remove temp files.
rm -r $tempdir
