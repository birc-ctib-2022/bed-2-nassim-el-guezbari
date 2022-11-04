"""Tool for cleaning up a BED file."""

import argparse  # we use this module for option parsing. See main for details.

import sys
from bed import (
    read_bed_file, print_line, BedLine
)


def extract_region(features: list[BedLine], start: int, end: int) -> list[BedLine]:
    """Extract region chrom[start:end] and write it to outfile."""
    up_bound = len(features)
    low_bound = 0
    while up_bound > low_bound:
        i = (up_bound + low_bound) // 2
        if features[i].chrom_start < start:
            low_bound = i + 1
        else:
            up_bound = i
    if low_bound >= len(features):
        return []
    feats=[]
    while features[low_bound].chrom_start < end:
        feats.append(features[low_bound])
        low_bound += 1
        if low_bound == len(features):
            break
    return feats


def main() -> None:
    """Run the program."""
    # Setting up the option parsing using the argparse module
    argparser = argparse.ArgumentParser(
        description="Extract regions from a BED file")
    argparser.add_argument('bed', type=argparse.FileType('r'))
    argparser.add_argument('query', type=argparse.FileType('r'))

    # 'outfile' is either provided as a file name or we use stdout
    argparser.add_argument('-o', '--outfile',  # use an option to specify this
                           metavar='output',  # name used in help text
                           type=argparse.FileType('w'),  # file for writing
                           default=sys.stdout)

    # Parse options and put them in the table args
    args = argparser.parse_args()

    # With all the options handled, we just need to do the real work
    features = read_bed_file(args.bed)
    for query in args.query:
        chrom, start, end = query.split()
        # Extract the region from the chromosome, using your extract_region()
        # function. If you did your job well, this should give us the features
        # that we want.
        region = extract_region(
            features.get_chrom(chrom), int(start), int(end))
        for line in region:
            print_line(line, args.outfile)


if __name__ == '__main__':
    main()
