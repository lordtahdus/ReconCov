#!/bin/bash
#SBATCH --job-name=sim_run
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Load the R module
module load r/4.4.1

export R_LIBS_USER=~/yi61/tsuu0007/R/library
# Run the script
Rscript ~/yi61/tsuu0007/ReconCov/sim/job/test_job.R

