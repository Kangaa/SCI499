include("scripts/OLD/MixingMatricesOLD2.jl")


using CSV
using DataFrames
SA2_codePop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
Codes = expandCodes(SA2_codePop[:, 1])
MM = SpatialMixingMatrix(SA2_codePop[:, 1], SA2_codePop[:,3], 0.5, 1 )

using StatsPlots
heatmap(log10.(MM))

SA3_codePop = CSV.read("data/GmelbSA3Pop21.csv", DataFrame)

transform!(SA3_codePop, :SA3_Code => ByRow(x -> string(x)) => :SA3_Code)
Codes = expandCodes(SA3_codePop[:, 1])
MM = SpatialMixingMatrix(SA3_codePop[:, 1], SA3_codePop[:,2], 0.5, 1 )


SA4_CodedPop = expandCodes(SA3_codePop[:, 1]) |> x ->
leftjoin(x, SA3_codePop, on = (:SA3 => :SA3_Code)) |> x ->
groupby(x, :SA4)|> x -> 
combine(x, :SA3_pop => sum => :Pop)

MM4 = SpatialMixingMatrix(SA4_CodedPop[:, 1], SA4_CodedPop[:,2], 0.5, 1 )

heatmap(log10.(MM4))