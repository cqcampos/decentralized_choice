#!/bin/bash
  
#SBATCH --account=pi-cqcampos
#SBATCH --partition=highmem
#SBATCH --mem=600GB       
#SBATCH --cpus-per-task=32
#SBATCH --time=4-00:00:00
#SBATCH --job-name=magnet_model
#SBATCH --output=/project/lausd/decentralized_choice/code/2_structural_estimation/logs/estimation_%j.log 
#SBATCH --mail-user=ryan.lee2@chicagobooth.edu  # Replace with your email
#SBATCH --nodelist=mcn62

module load python/booth/3.10

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"


# Run the python script
echo "Starting estimation at $(date)"
srun python3 estimate_mixture.py
echo "Job completed at $(date)"

