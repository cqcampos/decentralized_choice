#!/bin/bash
  
#SBATCH --account=faculty 
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --mem=640G         # 20GB per core for 64 cores
#SBATCH --time=6-12:00:00
#SBATCH --job-name=magnet_model_par2
#SBATCH --output=estimation_%j.log  # %j will be replaced by the job ID
#SBATCH --mail-type=BEGIN,END,FAIL  # Email notifications
#SBATCH --mail-user=Christopher.Campos@chicagobooth.edu  # Replace with your email

module load R/4.3/4.3.2

# Set up personal R library
export R_LIBS_USER=$HOME/R/library
mkdir -p $R_LIBS_USER

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

# Install required packages if they don't exist
Rscript -e "pkgs <- c('haven', 'doParallel', 'foreach', 'numDeriv', 'parallel'); 
            new_pkgs <- pkgs[!pkgs %in% installed.packages(lib.loc='$R_LIBS_USER')[,'Package']]; 
            if(length(new_pkgs)>0) install.packages(new_pkgs, lib='$R_LIBS_USER', repos='https://cran.rstudio.com/')"


# Set up R environment variables to use the correct number of cores
export MC_CORES=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


# Run the R script
echo "Starting estimation at $(date)"
Rscript --vanilla estimate.R

echo "Job completed at $(date)"

