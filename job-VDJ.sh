#!/bin/bash -l
#SBATCH --output=/scratch/prj/cb_hormad1_sc/BD-pipeline/BD-pipeline-VDJ.out
#SBATCH --job-name=BD-pipeline-VDJ
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=80G
#SBATCH --time=0-48:00
#SBATCH --error=/scratch/prj/cb_hormad1_sc/BD-pipeline/BD-pipeline-VDJ.err

cwl-runner \
--basedir ./ \
--cachedir ./cache_bd_vdj/ \
--tmpdir-prefix ./tmp_bd_vdj/ \
--leave-tmpdir \
--parallel \
--singularity \
--outdir results-VDJ \
--debug \
workflow.cwl \
input-VDJ.yml
