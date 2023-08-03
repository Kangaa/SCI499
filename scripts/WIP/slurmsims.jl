#Distributed SIR sims
using DrWatson
@quickactivate :SCI499

#load packages
using Distributed, ClusterManagers
using CSV, DataFrames

#load data (in main process)
Vicpops = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)

codes = Vicpops[:,1]
pops = vec(Vicpops[:,3])
intra = 0.5

#calculate mixing matrix
mixmat = SpatialMixingMatrix(codes, intra)


# define simulation parameters as dictionary
params = Dict(
    patch_names => Vicpops[:,1],
    population_per_patch => Vicpops[:,3],
    mixing_matrix => mixmat,
    β => 0.5,
    γ => 0.1
)

#create some worker processes
addprocs(2)

r = DataFrame()
#Pass mixing matrix to all processes
@distributed for i = 1:5
    simulate!(params)
end

#run simulations

#fetch results from worker processes and store in data structure (?)

#save results in BSON format on disk


