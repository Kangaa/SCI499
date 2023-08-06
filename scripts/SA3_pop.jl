using Pkg
Pkg.activate(".")

using Shapefile
using CSV
using DataFrames
using Tables
using DataFramesMeta

shps = Shapefile.Table("data\\ASGS_GDA2020\\SA3_2021_AUST_SHP_GDA2020\\SA3_2021_AUST_GDA2020.shp")
Melb_shp = Tables.subset(shps, shps.GCC_NAME21 .== "Greater Melbourne")

VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)

SA3_VicPop = @chain VicPop begin
  @transform @byrow codestring = string(:("SA2_Code"))
  @transform @byrow :SA3_Code = first(:codestring, 5)
  groupby(:SA3_Code)
  @combine :SA3_pop = sum(:("Pop")) 
  leftjoin(DataFrame([:SA3_Name =>Melb_shp.SA3_NAME21, :SA3_Code => Melb_shp.SA3_CODE21]), on = :SA3_Code)
  select(:SA3_Code, :SA3_Name, :SA3_pop)
end

CSV.write("data/GmelbSA3Pop21.csv", SA3_VicPop) 

SA4_VicPop = @chain VicPop begin
  @transform @byrow codestring = string(:SA2_Code)
  @transform @byrow :SA4_Code = first(:codestring, 3)
  groupby(:SA4_Code)
  @combine(:SA4_pop = sum(:Pop))
  leftjoin(DataFrame([:SA4_Name =>unique(Melb_shp.SA4_NAME21), :SA4_Code => unique(Melb_shp.SA4_CODE21)]), on = :SA4_Code)
  select(:SA4_Code, :SA4_Name, :SA4_pop)
end
CSV.write("data/GmelbSA4Pop21.csv", SA4_VicPop) 