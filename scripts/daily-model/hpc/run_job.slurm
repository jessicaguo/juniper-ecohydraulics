#!/bin/bash
#SBATCH --job-name=R-Plotting-Job
#SBATCH --ntasks=1
#SBATCH --nodes=1 
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=dlebauer

module load R/4.0.0
Rscript run_model.r