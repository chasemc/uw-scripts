import base64
import hashlib
from Bio.SeqUtils import CheckSum

import sys

def sha512t24u(input):
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7714221/
    # To standardize to caps-only input, use hash_aminos()
    sha512_digest = hashlib.sha512(bytes(input,"utf8")).digest()[:24]
    sha512t24u = base64.urlsafe_b64encode(sha512_digest).decode("ascii")
    return sha512t24u


def hash_aminos(input, **kwargs):
    # make sure everything is uppercase before hashing
    return hasher(input=input.upper(), **kwargs)


def use_hashlib(input, algo):
    allow_algos = ("sha512", "sha256", "sha384", "md5", "sha224")
    if algo not in allow_algos:
        # TODO: this should also output the sha512t24u ans crc algos
        raise ValueError(f"algo must be one of: {allow_algos}")
    hasher = getattr(hashlib, algo)
    return hasher(bytes(input, "utf8")).hexdigest()


def hasher(input, algo):
    match algo:
        case "sha512t24u":
            return sha512t24u(input)
        case "crc64":
            return CheckSum.crc64(input).removeprefix("CRC-")
        case "crc32":
            return CheckSum.crc32(input)
        case "seguid":
            return CheckSum.seguid(input)
        case _:
            return use_hashlib(input=input, algo=algo)


def do_the_thing(id, input):
    return(
            id,
            str(hasher(input,"seguid")),
            str(hasher(input,"crc64")),
            str(hasher(input,"crc32")),
            str(hasher(input,"md5")),
            str(hasher(input,"sha512")),
            str(hasher(input,"sha512t24u")),
        )

import fileinput

with open(sys.argv[1], "w") as h:
    current_id=""
    current_aa=""
    for line in fileinput.input():
        if 'Exit' == line.rstrip():
            break
        if line.startswith(">"):            
            if current_id:
                tmp="\t".join(do_the_thing(current_id, current_aa))
                h.write(tmp + "\n")
            current_id=""
            current_aa=""
            current_id=line.split(" " ,1)[0].removeprefix(">")
        else:
            current_aa= current_aa + line.strip()
