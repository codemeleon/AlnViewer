#!/usr/bin/env python

import sys
from Bio import SeqIO, AlignIO


def aln_format_predict(infile):
    """
        This function tries to find the format based in the file extesion. If
        file extesion is not given, it will return an error.
    """
    extensions = {
        'fasta': 'fasta',
        'fa': 'fasta',
        'fas': 'fasta',
        'phy': 'phylip',
        'aln': 'clustal'
    }

    file_extenstion = infile.split('.')[-1]
    if file_extenstion in extensions:
        return extensions[file_extenstion]
    else:
        print("Supported file extension not found")
        sys.exit()


def alignment_dict(infile):
    return SeqIO.to_dict(AlignIO.parse(open(infile),
                                       aln_format_predict(infile)))


if __name__ == '__main__':
    alignment_dict('infile.txt')  # Write a nose test
