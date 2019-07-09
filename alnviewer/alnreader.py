#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : alnreader.py
# Date              : 07.01.2019
# Last Modified Date: 07.01.2019
# -*- coding: utf-8 -*-
# File              : alnviewer/alnreader.py
# Date              : 23.10.2018
# Last Modified Date: 23.10.2018

"""Function part of the program"""


from Bio import SeqIO # , AlignIO


def aln_format_predict(infile):
    """
        This function tries to find the format based in the file extesion. If
        file extesion is not given, it will return an error.
    """
    extensions = {
        'fasta': 'fasta',
        'fa': 'fasta',
        'faa': 'fasta',
        'fna': 'fasta',
        'fas': 'fasta',
        'phy': 'phylip',
        'aln': 'clustal'
    }

    file_extenstion = infile.split('.')[-1]
    if file_extenstion in extensions:
        return extensions[file_extenstion]
    return None


def alignment_dict(infile):
    """Generates alignment sequence index based on indentified file format"""
    aln_format = aln_format_predict(infile)
    if aln_format:
        return SeqIO.index(infile, aln_format)
    return None
    # return SeqIO.to_dict(AlignIO.parse(open(infile),
                                       # aln_format_predict(infile)))


if __name__ == '__main__':
    alignment_dict('infile.txt')  # Write a nose test
