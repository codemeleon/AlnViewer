#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : alnviewer/alterations.py
# Date              : 23.10.2018
# Last Modified Date: 23.10.2018
# -*- coding: utf-8 -*-
# File              : alnviewer/alterations.py
# Date              : 23.10.2018
# Last Modified Date: 23.10.2018
# -*- coding: utf-8 -*-



import pickle
import re
from os import path

import pandas as pd


def aln_to_df(original):
    "Coneverts alignment to DataFrame"
    sequences = dict()
    for key, val in original.items():
        sequences[key] = val
    # sequences = dict([(k, list(v)) for k, v in original.items()])
    sequences = pd.DataFrame.from_dict(sequences)
    return sequences


def df_to_aln(seq_df):
    """Convert DataFrame to Sequence dictionary"""
    seq_df = pd.DataFrame.to_dict(seq_df, orient='list')
    seq_dict = dict()
    for key, val in seq_df.items():
        seq_dict[key] = ''.join(val)
    return seq_dict



def trim_map(original, gpc, hpc):
    # NOTE: gpc is gap percent and hpc is hetero percent
    """Return data fram with location where there is more than one kind
    nucleotide or amino acid"""
    # You might consider to return indexes and sequences
    sequences = aln_to_df(original)
    # TODO: Might need some corrections here
    sequences = sequences[(sequences == '-'
                           ).astype(int
                                    ).sum(axis=1)/sequences.shape[1] < gpc]
    sequences = sequences[sequences.apply(lambda x: x.tolist().count(
        x.mode().values[0]), axis=1)/sequences.shape[1] > hpc]
    # sequences = sequences.mode(axis=1)
    sequences = sequences.loc[sequences.apply(pd.Series.nunique, axis=1) != 1]
    indexes = sequences.index + 1
    indexes = list(map(str, indexes))
    max_len = max(map(len, indexes))
    indexes = ["."*(max_len-len(n))+n for n in indexes]
    sequences = df_to_aln(sequences)
    return sequences, indexes, max_len


def ref_map(original, ref):
    """Replace base/amino with . if it matches reference"""
    sequence_ids = list(original.keys())
    sequences = aln_to_df(original)
    ref_seq = sequences[ref]
    sequence_ids.remove(ref)
    del sequences[ref]
    for col in sequences.columns:
        sequences.loc[(sequences[col] == ref_seq) & (ref_seq != '-'),
                      col] = '.'
    sequences = df_to_aln(sequences)
    sequence_ids.sort()
    sequence_ids.insert(0, ref)
    sequences[ref] = ''.join(ref_seq)
    return sequences, sequence_ids


def pattern_ranges(sequences, pattern):
    """Indentifying pattern range"""
    search_movement = []
    pattern_pos = {}
    sequence_ids = list(sequences.keys())
    sequence_ids.sort()
    pat_len = len(pattern)
    for i, k in enumerate(sequence_ids):
        pattern_pos[k] = [[m.start(), m.start() + pat_len]
                          for m in re.finditer(pattern,
                                               str(sequences[k]))]
        for p in pattern_pos[k]:
            search_movement.append((i, p))
    return pattern_pos, search_movement


def old_status(current_file, hiddenpath):
    """Back to old parameters"""
    try:
        hidden_path = path.split(current_file)
        hidden_path = "%s/%s" % (hiddenpath,
                                 hidden_path[1])
        return pickle.load(open(hidden_path, "rb"))

    except Exception as ex:
        return {
            'display_mode': 'o',
            'display_color': 'normal',
            'display_ref': None,
            'display_x': 0,
            'display_y': 0
        }
