#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH -c 4 
#SBATCH --mail-type=ALL
#SBATCH --job-name=pepseq_design
#SBATCH --mail-user=ekelley@tgen.org
#SBATCH --nice=50
##SBATCH --begin=now+22hours
#SBATCH --partition=compute

# conda activate pepseq_encoding

/home/ekelley/bin/Library-Design/oligo_encoding/main \
    -r output_ratio \
    -s out_seqs \
    -n 300 \
    -c 2 \
    -p codon_weights_updated.csv \
    -i DP3_named_peptides.csv \
    -t 10000 \
    -g 0.55 \
