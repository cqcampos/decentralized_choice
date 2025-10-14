#!/bin/bash
  
#SBATCH --account=pi-cqcampos
#SBATCH --partition=gpu_h100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=counterfactual_gpu
#SBATCH --output=/project/lausd/decentralized_choice/code/logs/gpu_counterfactual_%j.log
#SBATCH --mail-user=ryan.lee2@chicagobooth.edu

nvidia-smi   # should show a GPU

source ~/venv/pytorch-gpu/bin/activate


# ---- # On host machine (mgpu) # ---- 
# set path to singularity container 
container=/apps/containers/pytorch-notebook-cuda12.sif 

# get path to R_LIBS_USER inside container 
rlibs_user_dir=$(singularity exec ${container} printenv R_LIBS_USER)
echo ${rlibs_user_dir} 

# create R_USER_LIBS path if not exists 
mkdir -p ${rlibs_user_dir} 

# prepare environment 
unset LD_PRELOAD 

# launch container with exec command apptainer 
exec --nv ${container} R 


# ---- # Inside container # ---- 
kind <- "cu124" 
version <- "0.15.1" 
url <- sprintf("https://torch-cdn.mlverse.org/packages/%s/%s/src/contrib/torch_%s_R_x86_64-pc-linux-gnu.tar.gz", kind, version, version) 
install.packages(url) 
library("torch") 
packageVersion("torch") 
torch::cuda_is_available()

#module avail

#module load R/4.3/4.3.2
#module  load cuda/12.4


echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"
echo "Node Name: $SLURM_NODELIST"

# Set up personal R library
#export R_LIBS_USER=/project/lausd/decentralized_choice/code/Rlib


# Install required packages if they don't exist
#Rscript -e "pkgs <- c('torch', 'haven', 'doParallel', 'foreach', 'dplyr', 'tidyr', 'matrixStats'); 
#            new_pkgs <- pkgs[!pkgs %in% installed.packages(lib.loc='$R_LIBS_USER')[,'Package']]; 
#            if(length(new_pkgs)>0) install.packages(new_pkgs, lib='$R_LIBS_USER', repos='https://cran.rstudio.com/')"


# Run the R script
#echo "Starting estimation at $(date)"
#srun Rscript --vanilla run_counterfactual.R

# echo "Job completed at $(date)"
