
#!/bin/bash

# for i in `seq 128`
# do
#     curl "https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/fasta/active/uniparc_active_p{$i}.fasta.gz" |\
#         zcat |
#         python3 /media/socialgene_nvme/hashing_uniparc.py /media/socialgene_nvme/hash_uniprot/$i.tsv - 
# done


getit() {
  curl  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/fasta/active/uniparc_active_p{$1}.fasta.gz" > "uniparc_active_p$1.fasta.gz"
}
export -f getit


parallel -j 12 getit ::: $(seq 1 128)

#june28 13:06


