#!/usr/bin/env python

##########
# If you do this follow NCBI's guidelines of when to do large downloads. This script won't to check your available hard drive space and **RefSeq** alone will require > 1 Terabyte** 
##########

# Import Module
import argparse
import ftplib
from tqdm import tqdm
import requests
from pathlib import Path
import logging
import multiprocessing
import csv
import hashlib
from collections.abc import Generator
from typing import List
import numpy as np

logging.basicConfig(filename="info",
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.INFO)

# Fill Required Information
HOSTNAME = "130.14.250.12"
USERNAME = "anonymous"
PASSWORD = "@anonymous"

parser = argparse.ArgumentParser(
    description="Download lots of files from RefSeq/Genbank"
)

parser.add_argument(
    "--refseq",
    metavar="bool",
    help="Download from RefSeq",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction
)

parser.add_argument(
    "--genbank",
    metavar="bool",
    help="Download from Genbank",
    required=False,
    default=False,
    action=argparse.BooleanOptionalAction
)
parser.add_argument(
    "--cpus",
    metavar="int",
    help="number of processes to spawn, max 6",
    required=False,
    default=1
)

HEADERS = ['assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'isolate', 'version_status', 'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter', 'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path', 'excluded_from_refseq', 'relation_to_type_material', 'asm_not_live_date']


def assembly_report_url(x):
    match x:
        case "r":
             return "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
        case "g":
            return "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"



def chunk_a_list_with_numpy(input_list: List, n_chunks: int) -> Generator:
    """Chunk a list into n-lists

    Args:
        input_list (list): input list to be chunked
        n_chunks (int): integer of number chunks to split list into
    """
    if len(input_list) < n_chunks:
        logging.info(
            "chunk_a_list_with_numpy(): n_chunks is < len(input_list), proceeding with len(input_list)"
        )
        n_chunks = len(input_list)
    return (list(i) for i in np.array_split(input_list, n_chunks))


def check_md5(md5_filepath, filename_to_check, md5_to_check):    
    expected_md5 = None
    with open(md5_filepath.resolve(), "r") as f:
        for line in f:
            if Path(line.split()[1]).name == filename_to_check:
                expected_md5 = line.split()[0].strip()
        if not expected_md5:
            logging.warning(f"{filename_to_check}\tmd5 not found in md5checksums.txt")
    if expected_md5 == md5_to_check:
        pass
    else:
        logging.warning(f"{filename_to_check}\tmd5 mismatch")
        gbff_filepath.unlink()



def check_or_delete(ftp_host, gbff_filepath ,md5_filepath):
    if md5_filepath.is_file() and gbff_filepath.is_file():
        logging.info(f"Skipping {gbff_filepath.name}")
    else:
        Path(gbff_filepath).parent.mkdir(parents=True, exist_ok=True)
        with open(md5_filepath, 'wb') as f:
            _=ftp_host.retrbinary(f"RETR {md5_filepath}", f.write)
        m = hashlib.md5()
        with open(gbff_filepath, 'wb') as f:
            def callback(data,):
                m.update(data)
                f.write(data)
            ftp_host.retrbinary('RETR %s' % gbff_filepath, callback)
        check_md5(md5_filepath=md5_filepath, filename_to_check=gbff_filepath.name, md5_to_check=m.hexdigest())
        
        


def download_records_generator(url):
    r = requests.get(url, stream=True)
    headers=HEADERS
    for line in r.iter_lines():
        # filter out keep-alive new lines
        assembly_summary = None
        if line:
            dl = line.decode('utf-8')
            # if dl.startswith("# assembly_accession"):
            #     headers = dl.removeprefix("#").strip().split('\t')
            if not dl.startswith("#"):
                assembly_summary = {x:y for x,y in zip(headers, dl.split('\t'))}
        if assembly_summary:
            yield assembly_summary

def downloader(records):
    pass_fail ={}
    with ftplib.FTP(HOSTNAME, USERNAME, PASSWORD) as ftp_host:
        for record in tqdm(records):
            pass_fail[record['assembly_accession']] = {}
            dirpath = Path(record['ftp_path'].removeprefix("https://ftp.ncbi.nlm.nih.gov/"))
            gbff_filepath = Path(dirpath, f"{dirpath.name}_genomic.gbff.gz")
            md5_filepath = Path(dirpath, "md5checksums.txt")
            gbff_filepath.parent.mkdir(parents=True, exist_ok=True)
            check_or_delete(ftp_host=ftp_host,gbff_filepath=gbff_filepath,md5_filepath=md5_filepath)
    return pass_fail



def main():
    args = parser.parse_args()

    if not args.refseq or args.genbank:
        raise argparse.ArgumentTypeError('either--refseq OR --genbank must be used')

    if args.refseq and args.genbank:
        raise argparse.ArgumentTypeError('either--refseq OR --genbank must be used')
    
    if args.refseq:
        db_from = "r"
    elif args.genbank:
        db_from = "g"


    logging.info("starting getting assembly info")
    req = download_records_generator(assembly_report_url(db_from))
    logging.info("finished getting assembly info")

    if int(args.cpus) >= 6:
        connections_to_make=4
    elif int(args.cpus) > 0:
        connections_to_make = int(args.cpus)
    else:
        connections_to_make = 1

    if connections_to_make == 1:
        for record in req:
            downloader(req)
    else:
        record_list_of_lists_generator = chunk_a_list_with_numpy([i for i in req], connections_to_make)
        out_res = {}
        with multiprocessing.Pool(
                    processes=4) as pool:
            for zzz in tqdm(pool.imap_unordered(downloader, record_list_of_lists_generator), total=connections_to_make):
                out_res = out_res | zzz
        log.info(out_res)

if __name__ == "__main__":
    main()
