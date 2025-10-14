#!/bin/bash
  
#SBATCH --account=pi-cqcampos 
#SBATCH --partition=standard
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-12:00:00
#SBATCH --job-name=build_descentralized
#SBATCH --output=/project/lausd/decentralized_choice/code/0_build_data/logs/0_build_data_%j.log  # %j will be replaced by the job ID

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

module load stata/19.0

#srun stata-mp  -b  do /project/lausd/decentralized_choice/code/0_build_data/a_build.do
srun stata-mp  -b  do /project/lausd/decentralized_choice/code/0_build_data/b_build_mle_data.do
# srun stata-mp  -b  do /project/lausd/decentralized_choice/code/0_build_data/c_build_mle_data_2004_2013.do


