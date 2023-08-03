#!/bin/bash
#SBATCH --job-name=SIRsims

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1000
#SBATCH --time=00:30:00


module load julia
julia --project=/home/jbender/Documents/SCI499 -threads=1, Scripts/WIP/slurmsims2.jl
rm julia.out 


# 1 Node, ~cpus per node. 