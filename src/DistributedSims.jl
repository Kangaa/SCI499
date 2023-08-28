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
if ARGS[1] == "HMM"
    @everywhere const CI_params = (
        MM_type = "HMM",
        mm_parameter = "",
        beta = parse(Float64, ARGS[2]),
        gamma = parse(Float64, ARGS[3]),
        SA_scale =  ARGS[4],
        nsims = parse(Int64, ARGS[5]))::Tuple{String, Float64, Float64, Float64, String, Int64}
end

if ARGS[1] == "HPMM"
    @everywhere const CI_params = ( 
        MM_type = "HPMM",
        mm_parameter = parse(Float64, ARGS[2]),
        beta = parse(Float64, ARGS[3]),
        gamma = parse(Float64, ARGS[4]),
        SA_scale =  ARGS[5],
        nsims = parse(Int64, ARGS[6]))::Tuple{String, Float64, Float64, Float64, String, Int64}
end

const VicPop = CSV.read("data/Gmelb$(CI_params.SA_scale)Pop21.csv", DataFrame)
const patchnames = convert(Vector{String}, VicPop[:,2]) 
const codes = VicPop[:, 1]|> x -> string.(x)
const popns  = VicPop[:, 3] .+ 1

SA_scale == "SA2" ?  ξ = [1/4,1/4,1/4,1/4] : 
SA_scale == "SA3" ?  ξ = [1/2,1/4,1/4] :
SA_scale == "SA4" ?  ξ = [3/4,1/4] : error("SA_Scale must be SA2, SA3 or SA4")

if ARGS[1] == "HHMM"
    const    mixmat = SCI499.MixingMatrices.HierarchicalPopulationMixingMatrix(
        DataFrame(
            Codes = codes,
            Pop = popns),
        ξ, mm_parameter)::Matrix{Float64}
elseif ARGS[1] == "HPMM"
    const    mixmat = SCI499.MixingMatrices.HMixingMatrix(CodePops, ξ)::Matrix{Float64}
end




addprocs(40,   exeflags="--project")

@everywhere begin
    include("src/SCI499.jl")
    using .SCI499
    using CSV
    using DataFrames
end

@everywhere const metaparams = (
    mixmat_type = $(CI_params.MM_type),
    popns = $popns,
    codes = $codes,
    patchnames = $patchnames,   
    mm_parameter = $(CI_params.mm_parameter),
    mixmat = $mixmat,
    beta = $(CI_params.beta),
    gamma = $(CI_params.gamma),
    SA_scale = $(CI_params.SA_scale))::Tuple{Vector{Int64}, Vector{String}, Vector{String}, Float64, Matrix{Float64}, Float64, Float64, String}


@everywhere begin

const  sim_params = (
        patch_names = metaparams.patchnames,
        population_per_patch = metaparams.popns,
        mixing_matrix = metaparams.mixmat,
        β = metaparams.metabeta,
        γ = metaparams.gamma
    )
end

@everywhere begin
    #run sims
    pmap(function batch_sim(sim) 
        tot_data, summary_stats, patch_sus, patch_inf = CompartmentalModels.simulate(sim_params, 10, 1, false, sim)

        CSV.write("data/sims/$(metaparams.SA_scale)/$(metaparams.SA_scale)_$(metaparams.beta)_$(metaparams.gamma)_$(metaparams.mixmat_type)$(metaparams.mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
        CSV.write("data/sims/$(metaparams.SA_scale)/$(metaparams.SA_scale)_$(metaparams.beta)_$(metaparams.gamma)_$(metaparams.mixmat_type)$(metaparams.mm_parameter)_patchinf_$sim.csv" , patch_inf)
        CSV.write("data/sims/$(metaparams.SA_scale)/$(metaparams.SA_scale)_$(metaparams.beta)_$(metaparams.gamma)_$(metaparams.mixmat_type)$(metaparams.mm_parameter)_patchsus_$sim.csv" , patch_sus)

    end,
    1:ci_params.nsims)
end

rmprocs(workers())

