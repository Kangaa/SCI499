## saving mixing matrices
using CSV
using StatsPlots
using DataFrames
include("../src/SCI499.jl")
using .SCI499

## ODMMs

for δH in 0:0.05:1, SA in ["SA2", "SA3", "SA4"]
    OD = MixingMatrices.ODMixingMatrix(SA, δH)

    OD |>
     x -> Tables.table(x) |>
     x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/ODMM_$(SA)_$(δH).csv", x)
end


MixingMatrices.ODMixingMatrix("SA4", 0.1) |>
 x -> Tables.table(x) |>
 x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/ODMM_SA4_0.1.csv", x)


# SAMats
## SA2
 VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
 patchnames = convert(Vector{String}, VicPop[:,2]) 
 codes = VicPop[:, 1]|> x -> string.(x)
 popns  = VicPop[:, 3] .+ 1

ξ = [1/4,1/4,1/4,1/4]


MixingMatrices.HMixingMatrix(codes, ξ)  |>
x -> DataFrame([codes x], ["i"; codes])  |>
 x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/SAMM_SA2.csv", x)

MixingMatrices.HPMixingMatrix(
        DataFrame(
            Codes = codes,
            Pop = popns),
        ξ, 1.0) |>
         x -> DataFrame([codes x], ["i"; codes]) |>
         x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/HPMM_SA2_1.0.csv", x)
 
         MixingMatrices.HPMixingMatrix(
            DataFrame(
                Codes = codes,
                Pop = popns),
            ξ, 0.0) |>
            x -> DataFrame([codes x], ["i"; codes]) |>
             x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/PPMM_SA2.csv", x)

 ## SA3

 VicPop = CSV.read("data/GmelbSA3Pop21.csv", DataFrame)
 patchnames = convert(Vector{String}, VicPop[:,2]) 
 codes = VicPop[:, 1]|> x -> string.(x)
 popns  = VicPop[:, 3] .+ 1

ξ = [1/2,1/4,1/4]

MixingMatrices.HMixingMatrix(codes, ξ)  |>
x -> DataFrame([codes x], ["i"; codes]) |>
 x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/SAMM_SA3.csv", x)

MixingMatrices.HPMixingMatrix(
        DataFrame(
            Codes = codes,
            Pop = popns),
        ξ, 1.0) |>
        x -> DataFrame([codes x], ["i"; codes]) |>
        x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/HPMM_SA3_1.0.csv", x)

MixingMatrices.HPMixingMatrix(
        DataFrame(
            Codes = codes,
            Pop = popns),
        ξ, 0.0) |>
        x -> DataFrame([codes x], ["i"; codes]) |>
        x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/PPMM_SA3.csv", x)

### SA4
    VicPop = CSV.read("data/GmelbSA4Pop21.csv", DataFrame)
    patchnames = convert(Vector{String}, VicPop[:,2]) 
    codes = VicPop[:, 1]|> x -> string.(x)
    popns  = VicPop[:, 3] .+ 1

    ξ = [3/4,1/4]

    MixingMatrices.HMixingMatrix(codes, ξ)  |>
    x -> DataFrame([codes x], ["i"; codes]) |>
     x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/SAMM_SA4.csv", x)

    MixingMatrices.HPMixingMatrix(
            DataFrame(
                Codes = codes,
                Pop = popns),
            ξ, 1.0) |>
            x -> DataFrame([codes x], ["i"; codes]) |>
            x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/HPMM_SA4_1.0.csv", x)

    MixingMatrices.HPMixingMatrix(
            DataFrame(
                Codes = codes,
                Pop = popns),
            ξ, 0.0) |>
            x -> DataFrame([codes x], ["i"; codes]) |>
            x ->  CSV.write("Papers/Manuscript/data/MixingMatrices/PPMM_SA4.csv", x)









