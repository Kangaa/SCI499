##Create mixing matrices
using DrWatson
@quickactivate :SCI499
using .SCI499
using CSV
using DataFrames
### SA2

SA2_VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
SA2_codes = SA2_VicPop[:, 1]|> x -> string.(x)
SA2_names = SA2_VicPop[:,2]
SA2_popns = SA2_VicPop[:, 3] .+ 1  

SA2_CodePop = DataFrame([:SA2_Code => SA2_codes, :SA2_pop => SA2_popns .+ 1])
SA2_mixmat = SCI499.MixingMatrices.SpatialMixingMatrix(SA2_CodePop, [1/4, 1/4, 1/4, 1/4], 1.0)

CSV.write("data/SA2_Mixmat_0.5_0.5.csv", DataFrame(SA2_mixmat, SA2_codes))

## SA3 


using Tables
using DataFramesMeta
using Shapefile

shps = Shapefile.Table("data\\ASGS_GDA2020\\SA3_2021_AUST_SHP_GDA2020\\SA3_2021_AUST_GDA2020.shp")
SA3_melb = Tables.subset(shps, shps.GCC_NAME21 .== "Greater Melbourne")

DataFrame([:SA3_Name =>SA3_melb.SA3_NAME21, :SA3_Code => SA3_melb.SA3_CODE21])



SA3_VicPop = @chain SA2_VicPop begin
  @transform @byrow codestring = string(:("SA2 code"))
  @transform @byrow :SA3_Code = first(:codestring, 5)
  groupby(:SA3_Code)
  @combine :SA3_pop = sum(:("Pop")) 
end

SA3_VicPop = innerjoin(SA3_VicPop,DataFrame([:SA3_Name =>SA3_melb.SA3_NAME21, :SA3_Code => SA3_melb.SA3_CODE21]), on = :SA3_Code)





SA3_codes = SA3_VicPop[:, 1]|> x -> string.(x)
SA3_names = SA3_VicPop[:,3]
SA3_popns = SA3_VicPop[:, 2] .+ 1  

SA3_mixmat = SCI499.MixingMatrices.SpatialMixingMatrix(SA3_codes, SA3_popns, 0.5, 1)

CSV.write("data/SA3_Mixmat_0.5_0.5.csv", DataFrame(SA3_mixmat, SA3_codes))