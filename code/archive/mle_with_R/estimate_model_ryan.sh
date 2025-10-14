#!/bin/bash
  
#SBATCH --account=pi-cqcampos
#SBATCH --partition=highmem
#SBATCH --mem=512GB       
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00:00
#SBATCH --job-name=magnet_model
#SBATCH --output=/project/lausd/decentralized_choice/code/logs/post_estimation_%j.log  # %j will be replaced by the job ID
#SBATCH --mail-user=ryan.lee2@chicagobooth.edu  # Replace with your email



module load R/4.3/4.3.2
module load python/booth/3.10
which python3

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"


# Set up personal R library
export R_LIBS_USER=$HOME/R/library
mkdir -p $R_LIBS_USER


# Install required packages if they don't exist
Rscript -e "pkgs <- c('haven', 'doParallel', 'foreach', 'numDeriv', 'matlib'); 
            new_pkgs <- pkgs[!pkgs %in% installed.packages(lib.loc='$R_LIBS_USER')[,'Package']]; 
            if(length(new_pkgs)>0) install.packages(new_pkgs, lib='$R_LIBS_USER', repos='https://cran.rstudio.com/')"


# Set up R environment variables to use the correct number of cores
export MC_CORES=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


# Run the R script
echo "Starting estimation at $(date)"

# srun Rscript estimate.R
# srun Rscript estimate_mixture.R
srun Rscript post_estimation.R

echo "Job completed at $(date)"

