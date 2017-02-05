#!/usr/bin/env python

import click
from alnview import viewer
from glob import glob
from os import path, makedirs


@click.command()
@click.option("--alnpath", help="Alignment file/folder",
              type=str, default="../AlnView_bkp/aln_fold",
              show_default=True)
@click.option("--alntype", help="nucleotide/Amino Acid",
              type=click.Choice(["nuc", "amn"]),
              default="amn", show_default=True)
def run(alnpath, alntype):
    """Commandline Sequence alignemnt viewer."""
    """Please prefer to use it in full screen mode."""
    if not path.exists(alnpath):
        click.echo("%s doesn't exist\nExiting ...")
        exit(1)
    file_list = []
    if path.isfile(alnpath):
        hidden_path = path.split(alnpath)[0] + "/.alnview"
        file_list = [alnpath]
    elif path.isdir(alnpath):
        hidden_path = alnpath + "/.alnview"
        file_list = glob("%s/*" % alnpath)
    else:
        click.echo("Stranger thing(s). Don't know what to do.")
        exit(1)
    if not path.exists(hidden_path):
        makedirs(hidden_path)
    if len(file_list):
        viewer.show_me_the_alignment(file_list, alntype)
    else:
        click.echo("No input file found in give folder")
        exit(1)

if __name__ == '__main__':
    run()
