
#!/bin/bash


getit() {

  cut -f $1 /media/socialgene_nvme/hash_uniprot/uniparc_hash_test.tsv | sort | uniq | wc -l > /media/socialgene_nvme/hash_uniprot/results/size_${1}

}
export -f getit


parallel -j 7 getit ::: $(seq 1 7)

#june28 13:06

wc -l /media/socialgene_nvme/hash_uniprot/uniparc_hash_test.tsv > /media/socialgene_nvme/hash_uniprot/results/size_truth


