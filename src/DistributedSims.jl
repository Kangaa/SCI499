## Activate project and get src
using Pkg
Pkg.activate("../SCI499")
Pkg.instantiate()
include("src/SCI499.jl")
using .SCI499

## load packages
using Distributed
using CSV
using DataFrames

## parse arguments

@everywhere begin 
    mm_parameter = parse(Float64, ARGS[1])
    beta = parse(Float64, ARGS[2])
    gamma = parse(Float64, ARGS[3])
    SA_scale =  ARGS[4]
    nsims = parse(Int64, ARGS[5])
end

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


addprocs(40,   exeflags="--project")

@everywhere begin
    include("src/SCI499.jl")
    using .SCI499
    using CSV
    using DataFrames
end

@everywhere begin
    popns = $popns
    codes = $codes 
    names = $names   
    mm_parameter = $mm_parameter
    mixmat = $mixmat
    beta = $beta
    gamma = $gamma
    SA_scale = $SA_scale
end

@everywhere begin
const  sim_params = (
        patch_names = names,
        population_per_patch = popns,
        mixing_matrix = mixmat,
        β = beta,
        γ = gamma
    )
nsims = $nsims

end

@everywhere begin
    #run sims
    pmap(function batch_sim(sim) 
        tot_data, summary_stats, patch_sus, patch_inf = CompartmentalModels.simulate(sim_params, 10, 1,false, sim)

        CSV.write("data/sims/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
        CSV.write("data/sims/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchinf_$sim.csv" , patch_inf)
        CSV.write("data/sims/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchsus_$sim.csv" , patch_sus)
        CSV.write("data/sims/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_summary_$sim.csv" , summary_stats)
    end,
    1:nsims)
end

rmprocs(workers())

