#!/usr/bin/env python
"""Main file for the task."""

import click
from alnviewer import viewer
from glob import glob
from os import path, makedirs


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
def run(alnpath, alntype, undocount, modified):
    """Commandline Sequence alignemnt viewer."""
    """Please prefer to use it in full screen mode."""
    if not alnpath or not path.exists(alnpath):
        click.echo("%s doesn't exist\nExiting ...")
        exit(1)
    file_list = []
    if path.isfile(alnpath):
        dirpath = path.split(alnpath)[0]
        if not len(dirpath):
            dirpath = '.'
        hiddenpath = dirpath + "/.alnview"
        file_list = [alnpath]
        # if modified == "modified":
        modified = dirpath + "/" + modified
    elif path.isdir(alnpath):
        hiddenpath = alnpath + "/.alnview"
        # if modified == "modified":
        modified = path.split(alnpath)[0] + "/" + modified
        file_list = glob("%s/*" % alnpath)
    else:
        click.echo("Stranger thing(s). Don't know what to do.")
        exit(1)
    if not path.exists(hiddenpath):
        makedirs(hiddenpath)
    if len(file_list):
        if modified.startswith("/"):
            modified = modified[1:]
        viewer.show_me_the_alignment(file_list, hiddenpath,
                                     alntype, undocount, modified)
    else:
        click.echo("No input file found in give folder")
        exit(1)


if __name__ == '__main__':
    run()
