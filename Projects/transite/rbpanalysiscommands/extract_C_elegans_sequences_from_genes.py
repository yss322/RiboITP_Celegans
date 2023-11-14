#!/usr/bin/env python3
"""Extracts sequences from a FASTA file based on gene name."""

import argparse
import os
import pysam
import sys

HEADER_DELIM = "|"
GENE_DELIM = ":"
GENE_FIELD_IDX = 6
GENE_IDX = 1
FOREGROUND_EXT = ".foreground.fa"
BACKGROUND_EXT = ".background.fa"


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    parser.add_argument("-i", "--gene_ids_file", type=str, required=True, help='Gene IDs to select for foreground.')

    parser.add_argument("-f", "--fasta", type=str, required=True, help='Curated transcriptome FASTA.')

    parser.add_argument("-o", "--output_dir", type=str, required=False, default=".",
                        help='Optional output directory. Default current working directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    outdir = parsed_args["output_dir"]
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    gene_file = parsed_args["gene_ids_file"]

    foreground_outfile = os.path.join(outdir, os.path.basename(os.path.splitext(gene_file)[0]) + FOREGROUND_EXT)
    background_outfile = os.path.join(outdir, os.path.basename(os.path.splitext(gene_file)[0]) + BACKGROUND_EXT)

    with open(gene_file, "r") as gene_fh:
        genes = {e.strip() for e in gene_fh}

    fasta_file = parsed_args["fasta"]

    with pysam.FastxFile(fasta_file, "r") as in_fh, \
            open(foreground_outfile, "w") as foreground_fh, \
            open(background_outfile, "w") as background_fh:

        for rec in in_fh:

            rec_split = rec.name.split(HEADER_DELIM)
            rec_gene_field = rec_split[GENE_FIELD_IDX]
            rec_gene = rec_gene_field.split(GENE_DELIM)[GENE_IDX]
            rna_seq = "".join(e if e != "T" else "U" for e in rec.sequence)
            rec.sequence = rna_seq

            if rec_gene in genes:
                foreground_fh.write(str(rec) + "\n")

            background_fh.write(str(rec) + "\n")


if __name__ == "__main__":
    main()
