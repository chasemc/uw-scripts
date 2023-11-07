from collections import defaultdict
import gzip
from rich.progress import (
    Progress,
    MofNCompleteColumn,
    TextColumn,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
    TextColumn,
    TimeRemainingColumn,
)
from sortedcontainers import SortedSet

res = defaultdict(set)

proteins = {}
hmms = {}

with Progress(transient=True) as pg:
    task = pg.add_task("Progress...", total=270923787)
    with gzip.open(
        "/media/socialgene_nvme/v0.2.3/refseq/socialgene_neo4j/import/protein_info/d32e30570390a1dc9108993e505fe438.protein_ids.gz",
        "r",
    ) as h:
        for i, line in enumerate(h):
            proteins[line.decode().strip().split("\t")[0]] = i
            pg.update(task, advance=1)

with Progress(transient=True) as pg:
    task = pg.add_task("Progress...", total=25562)
    with gzip.open(
        "/media/socialgene_nvme/v0.2.3/refseq/socialgene_neo4j/import/hmm_info/5ba584e2a9efcd9319d6b9f252102475.sg_hmm_nodes"
    ) as h:
        for i, line in enumerate(h):
            hmms[line.decode().strip()] = i
            pg.update(task, advance=1)

progress_bar = Progress(
    TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
    BarColumn(),
    MofNCompleteColumn(),
    TextColumn("•"),
    TimeElapsedColumn(),
    TextColumn("•"),
    TimeRemainingColumn(),
)

with progress_bar as pg:
    task = pg.add_task("Progress...", total=918179086)
    with open("/media/socialgene_nvme/v0.2.3/refseq/sparse.txt", "w") as hout:
        with gzip.open(
            "/media/socialgene_nvme/v0.2.3/refseq/socialgene_neo4j/import/parsed_domtblout/3db4fb69b9a5fdfb58eefc84426afb85.parseddomtblout.gz",
            "r",
        ) as h:
            for i, line in enumerate(h):
                res[proteins[line.decode().split("\t")[0]]].add(
                    hmms[line.decode().split("\t")[1]]
                )
                pg.update(task, advance=1)


# docker run \
#     --user=$(id -u):$(id -g) \
#     -p7474:7474 -p7687:7687 \
#     -v $sg_neoloc/data:/data \
#     -v $sg_neoloc/logs:/logs \
#     -v $sg_neoloc/import:/var/lib/neo4j/import \
#     -v $sg_neoloc/plugins:/plugins \
#     -v $sg_neoloc/conf:/var/lib/neo4j/conf \
#         --env NEO4J_AUTH=neo4j/test12345 \
#         --env NEO4J_PLUGINS='["apoc"]' \
#         --env NEO4J_dbms_security_procedures_unrestricted=algo.*,apoc.*,n10s.*,gds.*, \
#         --env NEO4J_dbms_security_procedures_allowlist=algo.*,apoc.*,n10s.*,gds.* \
#         --env NEO4J_server_config_strict__validation_enabled=false \
#         --env NEO4J_server_memory_heap_initial__size='100G' \
#         --env NEO4J_server_memory_heap_max__size='100G' \
#         --env NEO4J_server_memory_pagecache_size='800G' \
#         --env NEO4J_server_jvm_additional='-XX:+ExitOnOutOfMemoryError' \
#         --env NEO4J_ACCEPT_LICENSE_AGREEMENT='yes' \
#     neo4j:5.7.0-enterprise
