# biotools

## Usage

    # Convert gtf to bed12, bed6, intron ...
    ./gtfparser.py gtf
    # Calculate fragment sizes for mapped bam file.
    ./fragment_sizes.py -i bam -r bed12 -g gtf -o output
    # merge rows
    ./merge_row.py -i input -o output -header 0 -by 0,1,2