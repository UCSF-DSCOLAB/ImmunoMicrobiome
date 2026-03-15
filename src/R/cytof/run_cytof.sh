#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x-%j.out
##SBATCH --mail-type=END
##SBATCH --mail-user=abc@def.ghi

#set -e
#set -o nounset

singularity exec -B /scratch:/scratch singularity_images/RSingleCell/v3/RSingleCell.sif Rscript ImmunoMicrobiome/src/cytof/cyclone_script_cytof.R -c ${CONFIG}

