#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : alnviewer.py
# Date              : 07.01.2019
# Last Modified Date: 07.01.2019


#pylint: disable-msg=too-many-arguments
#pylint: disable-msg=no-value-for-parameter


"""Main file for the task."""

from glob import glob
from os import getcwd, makedirs, path

import click

from alnviewer import viewer


@click.command()
@click.option("--alnpath", help="Alignment file/folder",
              type=str, default=None,
              show_default=True)
@click.option("--alntype", help="nucleotide/Amino Acid",
              type=click.Choice(["nuc", "amn"]),
              default="amn", show_default=True)
@click.option("--undocount", help="Number of UnDo/ReDo count",
              type=int,
              default=20, show_default=True)
@click.option("--modified", help="Folder for modified file",
              type=str,
              default="modified", show_default=True)
@click.option("-gpc", help="Gap fraction to suppress", type=float,
              default=0.5, show_default=True)
@click.option("-hpc", help="Homogygous fraction to suppress", type=float,
              default=0.5, show_default=True)
def run(alnpath, alntype, undocount, modified, gpc, hpc):
    """Commandline Sequence alignemnt viewer.
    Please prefer to use it in full screen mode.
    This script assume that the file contains alignment.
    If it not an alignment, program fill ends with gaps to make them
    of same sizes."""
    if not alnpath or not path.exists(alnpath):
        exit("Given path %s doesn't exist or not given\nExiting . . . . . .")
    alnpath = path.abspath(alnpath)
    file_list = []
    if path.isfile(alnpath):
        dirpath = path.split(alnpath)[0]
        if not dirpath:
            dirpath = getcwd()
        hiddenpath = dirpath + "/.alnview"
        file_list = [alnpath]
        # if modified == "modified":
        # modification_memory = dirpath + "/" + modified
    elif path.isdir(alnpath):
        hiddenpath = alnpath + "/.alnview"
        # if modified == "modified":
        # modified = path.split(alnpath)[0] + "/" + modified
        file_list = glob("%s/*" % alnpath)
    else:
        exit("Given path is neithe file nor folder. Exiting . . . . .")
    if not path.exists(hiddenpath):
        makedirs(hiddenpath)
    if file_list:
        viewer.show_me_the_alignment(file_list, hiddenpath,
                                     alntype, undocount, modified,
                                     gpc, hpc)
    else:
        exit("No file found. Exiting . . . . .")


if __name__ == '__main__':
    run()
