#!/bin/bash

#SBATCH -J rf2_test
#SBATCH -p kuhlab
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=40g
#SBATCH -t 09:00:00
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1

source ~/.bashrc
module load gcc
module load cuda
conda activate RF2
python /proj/kuhl_lab/RosettaFold2/RoseTTAFold2/run/get_msa.py --input_dir inputs --custom_msa_path /proj/kuhl_lab/RosettaFold2/RoseTTAFold2/examples_longleaf/complex/mmseqs2/query_0/
python /proj/kuhl_lab/RosettaFold2/RoseTTAFold2/run/run_rf2.py -input_dir mmseqs2/
