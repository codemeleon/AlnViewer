# About

Alnview is a simple interactive command-line key binding based viewer/editor to display/edit pre-processed multiple sequence alignments. Alnviewer can be used locally or over remote login to visualize one or multiple alignments without leaving the command-line user interface.

A desire to be able to able to visualise large multiple sequence alignments, stored on a remote host, initiated the design of the program.

The program is being developed by the MiDEP group at the [Malawi Liverpool Wellcome Trust Clinical Research Programme](http://www.mlw.medcol.mw/index.php/microbial-ecology.html), as part of the activities of the [H3ABionet Consortium](http://www.h3abionet.org/).


# Dependecies

- Python 3.+ (Not tested on Python 2.+)
- Biopython
- Click
- Pandas

# Installation

python setup install

or

pip install git+https://github.com/codemeleon/AlnView.git

# Usage

alnview --alnpath [aln_file|aln_folder] --alntype [amn|nuc] --undocount num --modified path_for_altered_files


# Amino acid highligting

Amino acid color in **three** different modes. Only valid for normal mode where all amino acid sequences are being displayed.

1. All the amino acids will be displayed with distincs colors. (use **o** key)
2. Fobicity base grouping has beed done using : Polar, Hydrophobic, charged, uncharged. (use **b** key)
3. Sizes: large and small. (use **z** key)
4. *Your suggestions for groupping*


# Editing Options

## Removing bases
Press **d** and click on bases of any sequence to remove the base. From other sequences, bases from ends will be removed. When you are done, press **esc**.

## Removing columns
Press **D** and click on columns to remove columns. When you are done, press **esc**.

## Inserting gaps
Press **i** and click on a base of sequence to insert gap. For rest of the sequences, gaps will be appended at the end to keep the sequences of same length. When you are done, press **esc**.

## Deleting a sequence
Press **S** and click on a sequence to delete. When you are done, press **esc**.

## UnDo and ReDo
**ctrl+u** and **ctrl+r** for undoing and redoing above changes in normal mode.

# Key binding options

Please press **?** inside the program

![Normal View][i6]

# Screenshots

**Normal view**

![Normal View][i1]

**Reference view**

![Reference View][i2]

**Reference view**

![Mismatch View][i3]

**Genomic location search**

![Genomic Location Search][i4]

**Pattern search**

![Pattern search][i5]

# Licence

GPLv3

# TODOs

- [ ] Mostly speed improvements

[i1]: figures/NormalView.png "Normal View"
[i2]: figures/Reference_based.png "Reference View"
[i3]: figures/Mismaches.png "Mismatch View"
[i4]: figures/Genomic_Location.png "Genomic Location Search"
[i5]: figures/Pattern.png "Pattern Search"
[i6]: figures/Options.png "Options Page"
