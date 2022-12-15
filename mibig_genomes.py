# script to find and download genomes associated with mibig
from pathlib import Path
import json
from Bio import Entrez
import numpy as np
import math
from rich.progress import Progress
import os
import urllib

Entrez.email = os.environ['ENTREZ_EMAIL']
Entrez.api_key = os.environ['ENTREZ_API']
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


def get_assemblies(ids, download=True):
    print(f"found {len(ids)} ids")
    links = []
    for id in ids:
        # get summary
        summary = get_assembly_summary(id)
        # get ftp link
        url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
        if url == "":
            continue
        label = os.path.basename(url)
        # get the fasta link - change this to get other formats
        link = os.path.join(url, label + "_genomic.gbff.gz")
        print(link)
        links.append(link)
        if download == True:
            # download link
            urllib.request.urlretrieve(link, f"{Path(outdir)}/{label}.gbff.gz")
    return links

get_assemblies(set(temp))
