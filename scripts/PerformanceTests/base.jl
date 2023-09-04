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

function Distrosim(mm_parameter, beta, gamma, SA_scale, Intervention_type)
    VicPop = CSV.read("data/Gmelb$(SA_scale)Pop21.csv", DataFrame)
    codes::Vector{String} = VicPop[:, 1]::Vector{Int64} |> x -> string.(x)
    names::Vector{String} = VicPop[:, 2] |> x -> string.(x)
    popns::Vector{Int64}  = (VicPop[:, 3] .÷ 100) .+ 1
    CodePops = DataFrame(
        Codes = codes,
        Pop = popns)
    SA_scale == "SA2" ?  ξ = [9/10, 1/30, 1/30, 1/30] : 
    SA_scale == "SA3" ?  ξ = [9/10, 1/20, 1/20] :
    SA_scale == "SA4" ?  ξ = [9/10,1/10] : error("SA_Scale must be SA2, SA3 or SA4")
    mixmat = SCI499.MixingMatrices.HPMixingMatrix(CodePops, ξ, mm_parameter)
      sim_params = (
            patch_names = names::Vector{String},
            population_per_patch = popns::Vector{Int64},
            mixing_matrix = mixmat::Matrix{Float64},
            β = beta::Float64,
            γ = gamma::Float64,
        )
    #run sims
SCI499.simulate(sim_params, 10, 1, Intervention_type, false, 1)

end


using BenchmarkTools
using Profile
using Cthulhu

Distrosim(0.5, 1.4, 1.0, "SA3", "total")