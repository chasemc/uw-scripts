#!/usr/bin/python3

<!-- MIT License

Copyright (c) Chase M. Clark

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. -->

import gzip
import logging
import sqlite3
from pathlib import Path
from typing import IO
from collections.abc import Generator
import bz2
import gzip
import lzma
from contextlib import contextmanager
from enum import Enum, auto
from pathlib import Path
from itertools import islice

import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Create a database from a FASTA file")

parser.add_argument(
    "input",
    metavar="fasta_path",
    type=Path,
    help="Path to the FASTA file to read",
)
parser.add_argument(
    "output",
    metavar="output_path",
    type=Path,
    help="Path where SQLite file will be written",
)
parser.add_argument(
    "-b",
    "--batchsize",
    metavar="Batch size",
    type=int,
    help="Add x-sequences to SQLite at a time",
    required=False,
    default=50000,
)
parser.add_argument(
    "-m",
    "--inmem",
    default=False,
    help="Create the SQLite database in RAM first then copy the database to disk",
    required=False,
    action=argparse.BooleanOptionalAction,
)


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    uncompressed = auto()


def is_compressed(filepath: Path) -> Compression:
    """
    Determines the compression type of a file based on its signature.

    Args:
      filepath (Path): The `filepath` parameter is a `Path` object that represents the path to the file
    that you want to check for compression.

    Returns:
      The function `is_compressed` returns the type of compression used for the file specified by the
    `filepath` parameter. The possible return values are `Compression.gzip`, `Compression.bzip2`,
    `Compression.xz`, or `Compression.uncompressed`.
    """
    with open(filepath, "rb") as f:
        signature = f.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        else:
            return Compression.uncompressed


@contextmanager
def open_read(filepath: Path) -> IO:
    """
    The function `open_read` opens a file for reading, taking into account different compression
    formats.

    Args:
      filepath (Path): The `filepath` parameter is the path to the file that you want to open and read.
    It should be a string representing the file path.
    """
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        f = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        f = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        f = lzma.open(filepath, "rt")
    else:
        f = open(filepath, "rt")
    try:
        yield f
    finally:
        f.close()


def parse_fasta(filepath: str) -> Generator[str, str]:
    seq_id = ""
    with open_read(filepath) as h:
        seq_id = ""
        seq = ""
        for i in h:
            if i[0] == ">":
                if seq_id:
                    yield (seq_id, seq)
                seq_id = i[1:].strip()
                seq = ""
            else:
                seq += i.strip()
        yield (seq_id, seq)


def batched(iterable, n):
    "Batch data into tuples of length n. The last batch may be shorter."
    "https://docs.python.org/3/library/itertools.html#itertools-recipes"
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def create_db(batch_iterator, batch_size, outpath, inmem=True):
    sqliteout = Path(outpath)
    # files = list(Path("/home/chase/Downloads/a").glob("*.tsv"))
    if inmem:
        # build database in-memory then write out to disk at the end
        conn = sqlite3.connect(":memory:")
    else:
        # build database on disk
        conn = sqlite3.connect(sqliteout)
    c = conn.cursor()
    c.execute("""DROP TABLE IF EXISTS sequences""")
    c.execute(
        """
        CREATE TABLE sequences (
            hash CHAR(32) PRIMARY KEY,
            sequence VARCHAR(255) NOT NULL
        )""",
    )
    conn.commit()
    c.execute("""PRAGMA synchronous = 0""")
    c.execute("""PRAGMA journal_mode = OFF""")
    c.execute("""PRAGMA temp_store = MEMORY""")
    c.execute("""PRAGMA cache_size = 1000000""")
    c.execute("""PRAGMA locking_mode = EXCLUSIVE""")
    c.execute("""PRAGMA ignore_check_constraints = true""")

    cntr = 0
    for batch in batch_iterator:
        c.executemany("insert into sequences(hash, sequence) values (?,?)", batch)
        conn.commit()
        cntr += batch_size
        print(f"Added: {cntr} entries", end="\r")
    conn.close()
    if inmem:
        # save in-memory sqlite database to disk
        db_disk = sqlite3.connect(sqliteout)
        conn.backup(db_disk)
    logging.info("Processed {cntr} entries")


def main():
    args = parser.parse_args()

    batch_iterator = batched(parse_fasta(args.input), args.batchsize)

    create_db(
        batch_iterator=batch_iterator,
        batch_size=args.batchsize,
        outpath=args.output,
        inmem=args.inmem,
    )


if __name__ == "__main__":
    main()
