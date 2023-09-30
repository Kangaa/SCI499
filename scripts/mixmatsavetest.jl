## saving mixing matrices
using CSV
using StatsPlots
using DataFrames
include("../src/SCI499.jl")
using .SCI499

# SAMats
## SA2
VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
patchnames = convert(Vector{String}, VicPop[:,2]) 
codes = VicPop[:, 1]|> x -> string.(x)
popns  = VicPop[:, 3] .+ 1

ξ = [1/4,1/4,1/4,1/4]

MixingMatrices.HPMixingMatrix(
        DataFrame(
            Codes = codes,
            Pop = popns),
        ξ, 1.0) |>
        x -> DataFrame([codes x], ["i"; codes]) |>
        x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/HPMM_SA2_1.0.csv", x)
## SA3

 VicPop = CSV.read("data/GmelbSA3Pop21.csv", DataFrame)
 patchnames = convert(Vector{String}, VicPop[:,2]) 
 codes = VicPop[:, 1]|> x -> string.(x)
 popns  = VicPop[:, 3] .+ 1

ξ = [1/2,1/4,1/4]
MixingMatrices.HPMixingMatrix(
        DataFrame(
            Codes = codes,
            Pop = popns),
        ξ, 1.0) |>
        x -> DataFrame([codes x], ["i"; codes]) |>
        x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/HPMM_SA3_1.0.csv", x)
### SA4
    VicPop = CSV.read("data/GmelbSA4Pop21.csv", DataFrame)
    patchnames = convert(Vector{String}, VicPop[:,2]) 
    codes = VicPop[:, 1]|> x -> string.(x)
    popns  = VicPop[:, 3] .+ 1

    ξ = [3/4,1/4]
    MixingMatrices.HPMixingMatrix(
            DataFrame(
                Codes = codes,
                Pop = popns),
            ξ, 1.0) |>
            x -> DataFrame([codes x], ["i"; codes]) |>
            x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/HPMM_SA4_1.0.csv", x)
