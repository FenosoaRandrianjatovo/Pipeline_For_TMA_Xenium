#!/bin/bash
#SBATCH --job-name=run_TMA
#SBATCH --account=def-amadou
#SBATCH --mem=150G
#SBATCH --time=0-2:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --output=/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/Pipeline_For_TMA_Xenium/logs/run_pipeline_TMA_%j.out
#SBATCH --error=/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/Pipeline_For_TMA_Xenium/logs/run_pipeline_TMA_%j.err


# Output message
echo "Job is starting..."

# Load the necessary environment and Python module

module load StdEnv/2023
module load  gcc/12.3  r/4.4.0

echo "Rscript  is starting..."

Rscript pipeline_first_draft.R

echo "Job finished."
