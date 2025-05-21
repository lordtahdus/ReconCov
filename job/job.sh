#!/bin/bash
#SBATCH --job-name=sim_chunk
#SBATCH --output=job/logs/output_%A_%a.log
#SBATCH --error=job/logs/error_%A_%a.log
#SBATCH --array=1-2
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Load the R module
module load r/4.4.0-mkl

export R_LIBS_USER=~/yi61/tsuu0007/R/library

# Calculate simulation index
INDEX=$(( $SLURM_ARRAY_TASK_ID ))
# START=$(( ($SLURM_ARRAY_TASK_ID - 1) * 100 + 1 ))
# END=$(( $SLURM_ARRAY_TASK_ID * 100 ))

# Run the script
Rscript ~/yi61/tsuu0007/ReconCov/job/S16-4-1/run1/jobscript.R $INDEX
