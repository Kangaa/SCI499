#!/bin/bash
#SBATCH --job-name=SIRsims

#SBATCH --nodes=7
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-task=1000
#SBATCH --time=00:30:00


module load julia/1.8.9
julia --project=/home/jbender/Documents/SCI499 -threads=1, Scripts/SIRsims.jl
rm julia-*.out 