using Distributed, ClusterManagers

    using DrWatson
    @quickactivate :SCI499


using CSV, DataFrames
VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
codes = VicPop[:,1]|> x -> string.(x)
mixmat = SpatialMixingMatrix(codes)

params = (
    patch_names =  string.(VicPop[:,1]),
    population_per_patch = VicPop[:,3],
    mixing_matrix = mixmat,
    β = 0.5,
    γ = 0.1
)

CompartmentalModels.simulate(params, 10, 7,true)