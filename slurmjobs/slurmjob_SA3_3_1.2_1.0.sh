#!/bin/bash
#SBATCH --job-name=SIR_SA3_3_1.2_1.0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=1000
#SBATCH --time=20:00:00
#SBATCH --output=/home/jbender/Documents/SCI499/slurmout/%j.out

julia --project=/home/jbender/Documents/SCI499 DistributedSims_moss.jl 3 1.2 1.0 SA3 100
