#!/bin/bash
#SBATCH --job-name=SIR_SA2_1.0_1.2_1.0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=1000
#SBATCH --time=20:00:00
#SBATCH --output=/home/jbender/Documents/SCI499/slurmout/%j.out

julia --project=/home/jbender/Documents/SCI499 src/DistributedSims.jl 1.0 1.2 1.0 SA2 10
