#!/bin/bash
#SBATCH --account=pi-cqcampos
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=10
#SBATCH --mem=600GB
#SBATCH --time=4-00:00:00
#SBATCH --job-name=cfals
#SBATCH --output=logs/counterfactual_%j.log

module load R/4.3/4.3.2
module load gcc/9.2.0
module load python/booth/3.10
echo "Python binary: $(which python3)"


# Set up personal R library
export R_LIBS_USER=$HOME/R/library
mkdir -p $R_LIBS_USER

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

# Install required packages if they don't exist
Rscript -e "pkgs <- c('haven', 'doParallel', 'foreach', 'numDeriv', 'doRNG', 'doSNOW', 'Rmpi'); 
            new_pkgs <- pkgs[!pkgs %in% installed.packages(lib.loc='$R_LIBS_USER')[,'Package']]; 
            if(length(new_pkgs)>0) install.packages(new_pkgs, lib='$R_LIBS_USER', repos='https://cran.rstudio.com/')"

# Set up R environment variables to use the correct number of cores
export MC_CORES=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


# Run the R script for counterfactuals
echo "Starting estimation at $(date)"
srun Rscript --vanilla run_counterfactual.R
echo "Job completed at $(date)"