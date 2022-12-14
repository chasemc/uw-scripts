#!/usr/bin/env python3

import argparse
import logging
import os
import pandas as pd
import subprocess
import tempfile

parser = argparse.ArgumentParser(
    description="Run MMseqs2 search on an input amino acid sequence"
)

parser.add_argument(
    "input_fasta_path",
    help="Path to fasta file",
)

parser.add_argument(
    "target_database_path",
    help="Path to fasta file",
)

MMSEQS_OUT_COLUMNS = {
    "query": str,
    "target": str,
    "pident": "Float32",
    "alnlen": "Int32",
    "mismatch": "Int32",
    "gapopen": "Int32",
    "qstart": "Int32",
    "qend": "Int32",
    "tstart": "Int32",
    "tend": "Int32",
    "evalue": "Float32",
    "bits": "Int32",
}


def search(
    fasta_path,
    target_database,
):
    """Search an input fasta file against a target database with external MMseqs2 program

    Args:
        fasta_path (str): Path to fasta file
        target_database (str): Path to MMseqs2 database

    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        outpath = os.path.join(tmpdirname, "result.m8")
        command_list = [
            "mmseqs",
            "easy-search",
            fasta_path,
            target_database,
            outpath,
            tmpdirname,
            "--format-mode",
            "0",
        ]
        cmd = " ".join(command_list)
        logging.info(cmd)
        mes = subprocess.run(
            args=command_list, check=False, shell=False, capture_output=True
        )
        logging.info(mes)
        if not os.path.exists(outpath):
            logging.warning("No output from MMseqs2 search")
        else:
            return pd.read_csv(
                outpath, sep="\t", names=MMSEQS_OUT_COLUMNS, dtype=MMSEQS_OUT_COLUMNS
            )


def main():

    args = parser.parse_args()

    fasta_path = args.input_fasta_path
    target_database = args.target_database_path

    if not os.path.exists(fasta_path):
        raise FileExistsError(fasta_path)
    if not os.path.exists(target_database):
        raise FileExistsError(target_database)

    print(search(fasta_path, target_database))


if __name__ == "__main__":
    main()
