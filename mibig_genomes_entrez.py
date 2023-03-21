#!/usr/bin/env python3

# script to find and download genomes associated with mibig

# pip install biopython rich numpy


from pathlib import Path
import json
import math
import os
import urllib
import hashlib
import requests

from Bio import Entrez
import numpy as np
from rich.progress import Progress


Entrez.email = os.environ["ENTREZ_EMAIL"]
Entrez.api_key = os.environ["ENTREZ_API"]
json_dir = "/home/chase/Documents/data/mibig/3_1/mibig_json_3.1"
outdir = "/home/chase/Documents/data/mibig/3_1/mibig_genomes"


# get all json file paths
pathlist = Path(json_dir).glob("*.json")

# read json files and create dict  {mibig_id:locus_accession}
mibig_dict = {}
for path in pathlist:
    with open(path, "r") as f:
        data = json.load(f)
        mibig_dict[data["cluster"]["mibig_accession"]] = data["cluster"]["loci"][
            "accession"
        ]

chunked_list = np.array_split(
    list(mibig_dict.values()), math.ceil(len(list(mibig_dict.values())) / 10)
)

a = []

with Progress(transient=True) as progress:
    task = progress.add_task("Progress...", total=len(chunked_list))
    for id_list in chunked_list:
        for single_id in id_list:
            esearch_file = Entrez.esearch(db="nuccore", term=single_id, retmax=1)
            esearch_record = Entrez.read(esearch_file)
            temp = {}
            nucleo_ids = esearch_record["IdList"]
            results = Entrez.read(
                Entrez.elink(
                    dbfrom="nucleotide",
                    db="assembly",
                    LinkName="nuccore_assembly",
                    from_uid=nucleo_ids,
                )
            )
            for result in results:
                temp[single_id] = result["LinkSetDb"]
            a.append(temp)
        progress.update(task, advance=1)


temp = []
for i in a:
    for k, v in i.items():
        for ii in v:
            if ii["LinkName"] == "nuccore_assembly":
                for iii in ii["Link"]:
                    temp.append(iii["Id"])


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def extract_md5sum_from_filetext(url, filename):
    try:
        md5_req = requests.get(url, stream=True)
        md5 = [i.split("  ") for i in md5_req.text.split("\n")]
        md5 = [i for i in md5 if len(i) == 2]
        md5 = [i for i in md5 if i[1] == filename][0]
        return md5
    except:
        return [["this is will fail md5sum", "this is will fail md5sum"]]


def get_assemblies(ids, download=True, endswith="_genomic.gbff.gz"):
    print(f"Downloading {len(ids)} genomes")
    status = {"succeeded": [], "failed": []}
    for id in ids:
        # get summary
        summary = get_assembly_summary(id)
        # get ftp link
        url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
        if url == "":
            continue
        label = os.path.basename(url)
        # get the assembly link - change this to get other formats
        link = os.path.join(url, label + endswith)
        link = link.replace("ftp://", "https://", 1)
        md5_url = os.path.join(url, "md5checksums.txt")
        md5_url = md5_url.replace("ftp://", "https://", 1)
        if download == True:
            # download link
            outpath = Path(outdir, f"{label}.gbff.gz")
            genome_req = requests.get(link)
            expected_md5 = extract_md5sum_from_filetext(
                url=md5_url, filename=f"./{label + endswith}"
            )
            if hashlib.md5(genome_req.content).hexdigest() == expected_md5[0]:
                with open(outpath, "wb") as f:
                    f.write(genome_req.content)
                with open(Path(outdir, "md5sums"), "a") as f:
                    f.writelines(f"{expected_md5[0]}  {expected_md5[1]}\n")
                status["succeeded"].append(id)
            else:
                status["failed"].append(id)


get_assemblies(set(temp))
