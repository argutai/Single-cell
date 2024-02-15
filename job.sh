#!/bin/bash -l
#SBATCH --job-name=BD-pipeline
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=160G
#SBATCH --time=0-72:00
#SBATCH --error=/scratch/prj/cb_hormad1_sc/BD-pipeline/BD-pipeline.err
#SBATCH --output=/scratch/prj/cb_hormad1_sc/BD-pipeline/BD-pipeline.out

cwl-runner \
--basedir ./ \
--cachedir ./cache_bd_c1/ \
--tmpdir-prefix ./tmp_bd_c1/ \
--leave-tmpdir \
--parallel \
--singularity \
--disable-pull \
--outdir results \
workflow.cwl \
input.yml
