using DrWatson
@quickactivate 

using Shapefile
using Tables
using Plots

include(srcdir("SpatialSEIR_module.jl"))
using .SpatialSEIRutils.MixingMatrices

AUS_SA3_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA3_2021_AUST_SHP_GDA2020/SA3_2021_AUST_GDA2020.shp"))
GMelb_SA2_SHP = Tables.subset(AUS_SA3_SHP, AUS_SA3_SHP.GCC_NAME21 .== "Greater Melbourne")
Demo = Tables.subset(GMelb_SA3_SHP, 19:23)

SA_Codes = MixingMatrices.getCodes(GMelb_SA3_SHP)
SA4_Color_dict = Dict(zip(unique(SA_Codes[:,2]), palette(:tab10, length(unique(SA_Codes[:,2])))))
SA4_Colors = get.(Ref(SA4_Color_dict), SA_Codes[:,2], "NA")
plot(GMelb_SA3_SHP.geometry, color = SA4_Colors')


AUS_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp"))
GMelb_SA2_SHP = Tables.subset(AUS_SA2_SHP, AUS_SA2_SHP.GCC_NAME21 .== "Greater Melbourne")
SA_Codes = MixingMatrices.getCodes(GMelb_SA2_SHP)
SA3_Color_dict = Dict(zip(unique(SA_Codes[:,2]), palette(:tab10, length(unique(SA_Codes[:,2])))))
SA3_Colors = get.(Ref(SA3_Color_dict), SA_Codes[:,2], "NA")
plot(GMelb_SA2_SHP.geometry, color = SA3_Colors')
using ColorSchemes

SA3_names =  GMelb_SA2_SHP.SA3_NAME21
SA2_names =  GMelb_SA2_SHP.SA2_NAME21

Demo = Tables.subset(GMelb_SA2_SHP, findall(x -> (x == "Melbourne City" || x == "Essendon" || x == "Maribyrnong" ), SA2_names)) 
Demo = Tables.subset(GMelb_SA2_SHP, findall(
    x -> (
        x == "Flemington" ||
        x == "Ascot Vale" ||
        x == "West Melbourne - Industrial" ||
        x == "Flemington Racecourse" ||
        x == "Kensington (Vic.)" ||
        x == "Footscray"||
        x == "Maribyrnong"), SA2_names))

SA3_cols = Dict(zip(Demo.SA3_NAME21|>unique, palette(:tab10, Demo.SA3_CODE21|>unique) ))

plot(Demo.geometry, color = collect(Demo.SA3_NAME21.|> x -> get(colordict, x, "null"))')

test = MixingMatrices.Spatial_MixingMatrix2(GMelb_SA2_SHP, 0.7)

test2 = (test[:,:,1] + test[:,:,2] + test[:,:,3])/3

test[:,:,1][1,:]|>sum
test2[1,:]|>sum
using StatsPlots
heatmap(log.(test2))