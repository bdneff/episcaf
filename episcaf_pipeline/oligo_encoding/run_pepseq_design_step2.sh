#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH -c 4 
#SBATCH --mail-type=ALL
#SBATCH --job-name=pepseqdesign2
#SBATCH --mail-user=ekelley@tgen.org
#SBATCH --nice=50
##SBATCH --begin=now+22hours
#SBATCH --partition=compute

# conda activate pepseq_encoding 

/home/ekelley/bin/Library-Design/oligo_encoding/encoding_with_nn.py \
        -m /home/ekelley/bin/Library-Design/oligo_encoding/DeepLearning_model_R_1539970074840_1_20181019 \
        -r output_ratio \
        -s out_seqs \
        -o DP3_best_encodings \
	--subsample 300 \
	--read_per_loop 10 \
	-n 1 \
 
