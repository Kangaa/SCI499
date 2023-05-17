using Distributed, ClusterManagers


@everywhere using Pkg
@everywhere Pkg.activate("../SCI499")
@everywhere Pkg.instantiate()
@everywhere include("src/SCI499.jl")
@everywhere using .SCI499
addprocs(SlurmManager(4), exeflags=["--project"])

@everywhere begin    
using CSV, DataFrames

VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
codes = VicPop[:,1]|> x -> string.(x)
mixmat = MixingMatrices.SpatialMixingMatrix(codes)

params = (
    patch_names =  string.(VicPop[:,1]),
    population_per_patch = VicPop[:,3] .+1,
    mixing_matrix = mixmat,
    β = 0.5,
    γ = 0.1
)
end

@everywhere function batch_sim(sim)
    CompartmentalModels.simulate(params, 10, 7,true, sim)
end

pmap(batch_sim, 1:5)