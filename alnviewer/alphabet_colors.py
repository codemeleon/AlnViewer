#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : alnviewer/alphabet_colors.py
# Date              : 23.10.2018
# Last Modified Date: 23.10.2018
# -*- coding: utf-8 -*-

"""Set color or nucleotides and amino acids"""




import pandas as pd


def color_type(cl_type):
    """Assigns color"""
    amino_property_table = pd.DataFrame.from_dict(
        {
            'amino': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
                      'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                      'U', 'V', 'W', 'X', 'Y', 'Z', '-', '.', '*'],
            'phobicity': [20, 200, 80, 100, 100, 20, 20, 80, 20, 200, 100, 20,
                          80, 80, 20, 200, 80, 100, 80, 80,
                          200, 20, 80, 200, 80, 200, 231, 255, 255],
            # Include Polar, 1: 20
            # Hydrophobic 2: 80
            # and charged 3: 100, 4: 200
            'size': [100, 200, 100, 200, 200, 200, 100, 200, 200, 200,
                     200, 200, 200, 200, 200, 200, 200, 200, 100, 200,
                     200, 200, 200, 200, 200, 200, 231, 255, 255],
            'normal': [20, 3, 80, 11, 5, 6, 100, 110, 12, 13,
                       14, 21, 27, 28, 57, 58, 88, 94, 190, 200,
                       189, 201, 219, 226, 183, 160, 231, 255, 255]
        }
    )
    colors = {}
    # print(amino_property_table)
    for _, row in amino_property_table[['amino', cl_type]].iterrows():
        colors[row['amino']] = row[cl_type]
    return colors


def color_pairs(curses, mode):
    """Forground and background color paires"""
    # Modes are normal, seize and phobicity
    alphabet_colors = color_type(mode)
    color_pair = {}
    alternative_color_pair = {}
    for k in alphabet_colors:
        curses.init_pair(ord(k), -1, alphabet_colors[k])
        color_pair[k] = curses.color_pair(ord(k))
        # Check in 1000 is important
        curses.init_pair(ord(k)+100, alphabet_colors[k], -1)
        # curses.init_pair(ord(k), alphabet_colors[k], 80)
        alternative_color_pair[k] = curses.color_pair(ord(k)+100)
    return color_pair, alternative_color_pair
