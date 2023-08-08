## Next-Generation operator for spatial SIR model
## Mixing matrix

using CSV
using DataFrames

SA2_Mixmat = CSV.read("data/SA2_Mixmat_0.5_0.5.csv", DataFrame)|> x -> Matrix(x)
SA3_Mixmat = CSV.read("data/SA3_Mixmat_0.5_0.5.csv", DataFrame)|> x -> Matrix(x)


using LinearAlgebra
β = 1.9
γ = 0.5

## SA2
SA2_T = SA2_Mixmat * β
SA2_Σ = diagm(fill(-γ, size(SA2_Mixmat, 1)))
SA2_A = inv(SA2_Σ)
SA2_K = -(SA2_T*SA2_A)

maximum(abs.(eigvals(SA2_K)))

## SA3
SA3_T = SA3_Mixmat * β
SA3_Σ = diagm(fill(-γ, size(SA3_Mixmat, 1)))
SA3_A = inv(SA3_Σ)

SA3_K = -(SA3_T*SA3_A)

maximum(abs.(eigvals(SA3_K)))



##Plot
using StatsPlots
heatmap(SA2_Mixmat)
heatmap(SA3_Mixmat)


#Hss representation
using HssMatrices

SA2_hss_mixmat = hss(SA2_Mixmat)
SA3_hss_mixmat = hss(SA3_Mixmat)

norm(SA2_hss_mixmat) - norm(SA2_Mixmat)


using BenchmarkTools
SA2_ran_vec = rand(size(SA2_Mixmat, 1))
@benchmark SA2_Mixmat*SA2_ran_vec
@benchmark SA2_hss_mixmat*SA2_ran_vec


SA3_ran_vec = rand(size(SA3_Mixmat, 1))
@benchmark SA3_Mixmat*SA3_ran_vec
@benchmark SA3_hss_mixmat*SA3_ran_vec

plotranks(SA2_hss_mixmat)
plotranks(SA3_hss_mixmat)


test_1, test_2 = generators(SA2_hss_mixmat, (1,2))

cluster_tree = cluster(SA2_hss_mixmat)
cluster_tree

plot(cluster_tree)