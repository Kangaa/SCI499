## Activate project and get src
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/SCI499.jl")
using .SCI499

## load packages
using Distributed
using CSV
using DataFrames
using ProfileView
using BenchmarkTools

mm_parameter = 0.5
beta = 1.5
gamma = 1.0
SA_scale =  "SA2"
nsims = 1

VicPop = CSV.read("data/Gmelb$(SA_scale)Pop21.csv", DataFrame)
names = convert(Vector{String}, VicPop[:,2]) 
codes = VicPop[:, 1]|> x -> string.(x)
popns  = VicPop[:, 3] .+ 1

CodePops = DataFrame(
    Codes = codes,
    Pop = popns)

SA_scale == "SA2" ?  ξ = [1/4,1/4,1/4,1/4] : 
SA_scale == "SA3" ?  ξ = [1/2,1/4,1/4] :
SA_scale == "SA4" ?  ξ = [3/4,1/4] : error("SA_Scale must be SA2, SA3 or SA4")

mixmat = SCI499.MixingMatrices.SpatialMixingMatrix(CodePops, ξ, mm_parameter)

const  sim_params = (
        patch_names = names,
        population_per_patch = popns,
        mixing_matrix = mixmat,
        β = beta,
        γ = gamma
    )
using ProfileView
using Cthulhu
 #run sims
map(function batch_sim(sim) 
        tot_data, summary_stats, patch_sus, patch_inf = SCI499.simulate(sim_params, 10, 1,false, sim)

        CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
        CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchinf_$sim.csv" , patch_inf)
        CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchsus_$sim.csv" , patch_sus)
        CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_summary_$sim.csv" , summary_stats)
    end,
1:nsims)