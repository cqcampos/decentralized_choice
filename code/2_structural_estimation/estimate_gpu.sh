#!/bin/bash
  
#SBATCH --account=pi-cqcampos
#SBATCH --partition=gpu_h100
#SBATCH --gres=gpu:3
#SBATCH --cpus-per-task=24
#SBATCH --mem=242GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=mix_model_eta_out
#SBATCH --output=/project/lausd/decentralized_choice/code/2_structural_estimation/logs/estimation_gpu_%j.log 

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

module load python/booth/3.12
module load cuda/12.2
source ~/venv/pytorch-gpu/bin/activate


echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"
echo "Node Name: $SLURM_NODELIST"


# Run the estimation
echo "Starting estimation at $(date)"
srun python3 estimate_mixture.py
echo "Job completed at $(date)"

