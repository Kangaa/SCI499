#!/bin/bash
    #SBATCH --job-name=SIR_HPMM_SA4_0.0_total_1.8_1.0
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=40
    #SBATCH --mem=1000
    #SBATCH --time=20:00:00
    #SBATCH --output=/home/jbender/Documents/SCI499_new/slurmout/%j.out
    
    julia --project=/home/jbender/Documents/SCI499_new src/DistributedSims.jl HPMM 0.0 1.8 1.0 SA4 total 50
    