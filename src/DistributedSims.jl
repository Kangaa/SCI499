## Activate project and get src
using Pkg
Pkg.activate("../SCI499_new")
Pkg.instantiate()
include("SCI499.jl")
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
        Intervention_type = ARGS[5],
        nsims = parse(Int64, ARGS[6]))
elseif ARGS[1] == "HPMM"
    @everywhere const CI_params = ( 
        MM_type = "HPMM",
        mm_parameter = parse(Float64, ARGS[2]),
        beta = parse(Float64, ARGS[3]),
        gamma = parse(Float64, ARGS[4]),
        SA_scale =  ARGS[5],
        Intervention_type = ARGS[6],
        nsims = parse(Int64, ARGS[7]))
elseif ARGS[1] == "OD"
    @everywhere const CI_params = ( 
        MM_type = "OD",
        deltaH = parse(Float64, ARGS[2]),
        beta = parse(Float64, ARGS[3]),
        gamma = parse(Float64, ARGS[4]),
        SA_scale =  ARGS[5],
        Intervention_type = ARGS[6],
        nsims = parse(Int64, ARGS[7]))
end

const VicPop = CSV.read("data/Gmelb$(CI_params.SA_scale)Pop21.csv", DataFrame)
const patchnames = convert(Vector{String}, VicPop[:,2]) 
const codes = VicPop[:, 1]|> x -> string.(x)
const popns  = VicPop[:, 3] .+ 1

CI_params.SA_scale == "SA2" ?  ξ = [1/4,1/4,1/4,1/4] : 
CI_params.SA_scale == "SA3" ?  ξ = [1/2,1/4,1/4] :
CI_params.SA_scale == "SA4" ?  ξ = [3/4,1/4] : error("SA_Scale must be SA2, SA3 or SA4")

if ARGS[1] == "HPMM"
    const    mixmat = SCI499.MixingMatrices.HPMixingMatrix(
        DataFrame(
            Codes = codes,
            Pop = popns),
        ξ, mm_parameter)::Matrix{Float64}
        "generated HMM mixmat" 
elseif ARGS[1] == "HMM"
    const    mixmat = SCI499.MixingMatrices.HMixingMatrix(codes, ξ)::Matrix{Float64}
    "Generated HPMM Mixmat"

elseif ARGS[1] == "OD"
    const    mixmat = SCI499.MixingMatrices.ODMixingMatrix(CI_params.SA_scale, CI_params.deltaH)::Matrix{Float64}
    "Generated OD Mixmat"
end

addprocs(40,   exeflags="--project")

@everywhere begin
    include("SCI499.jl")
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
    SA_scale = $(CI_params.SA_scale),
    Intervention_type = $(CI_params.Intervention_type))


@everywhere begin

const  sim_params = (
        patch_names = metaparams.patchnames,
        population_per_patch = metaparams.popns,
        mixing_matrix = metaparams.mixmat,
        β = metaparams.beta,
        γ = metaparams.gamma
    )
end

@everywhere begin
    #run sims
    pmap(function batch_sim(sim) 
        tot_data, patch_sus, patch_inf = CompartmentalModels.simulate(sim_params, 10, 1, metaparams.Intevention_type, false, sim)

        CSV.write("data/sims/$(metaparams.SA_scale)/$(metaparams.SA_scale)_$(metaparams.beta)_$(metaparams.gamma)_$(metaparams.mixmat_type)_$(metaparams.mm_parameter)_$(metaparams.Intevention_type)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
        CSV.write("data/sims/$(metaparams.SA_scale)/$(metaparams.SA_scale)_$(metaparams.beta)_$(metaparams.gamma)_$(metaparams.mixmat_type)_$(metaparams.mm_parameter)_$(metaparams.Intevention_type)_patchinf_$sim.csv" , patch_inf)
        CSV.write("data/sims/$(metaparams.SA_scale)/$(metaparams.SA_scale)_$(metaparams.beta)_$(metaparams.gamma)_$(metaparams.mixmat_type)_$(metaparams.mm_parameter)_$(metaparams.Intevention_type)_patchsus_$sim.csv" , patch_sus)

    end,
    1:$CI_params.nsims)
end

rmprocs(workers())

