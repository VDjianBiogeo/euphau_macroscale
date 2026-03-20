#!/bin/bash
#SBATCH --partition=std
#SBATCH -n 15
#SBATCH --mem-per-cpu 32000
#SBATCH --mail-user=valentin.djian@sorbonne-universite.fr
#SBATCH --mail-type=ALL
#SBATCH --time=03-00:00:00

# load R module
module load R/4.4.1
# launch parallel execution instances of R running script file test-dompi.R
mpirun Rscript 10_parsim_ensemble_mod.R

exit 0