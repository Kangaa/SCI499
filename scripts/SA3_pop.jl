using DrWatson
@quickactivate :SCI499
using Shapefile
using CSV
using DataFrames
using Tables
shps = Shapefile.Table("data\\ASGS_GDA2020\\SA3_2021_AUST_SHP_GDA2020\\SA3_2021_AUST_GDA2020.shp")
SA3_melb = Tables.subset(shps, shps.GCC_NAME21 .== "Greater Melbourne")

DataFrame([:SA3_Name =>SA3_melb.SA3_NAME21, :SA3_Code => SA3_melb.SA3_CODE21])


using DataFramesMeta

VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)

VicPop = @chain VicPop begin
  @transform @byrow codestring = string(:("SA2 code"))
  @transform @byrow :SA3_Code = first(:codestring, 5)
  groupby(:SA3_Code)
  @combine :SA3_pop = sum(:("Pop")) 
end

SA3_VicPop = innerjoin(VicPop,DataFrame([:SA3_Name =>SA3_melb.SA3_NAME21, :SA3_Code => SA3_melb.SA3_CODE21]), on = :SA3_Code)

CSV.write("data/GmelbSA3Pop21.csv", SA3_VicPop) 

