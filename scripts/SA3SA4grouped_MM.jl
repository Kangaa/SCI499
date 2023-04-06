using DrWatson
@quickactivate 

using Shapefile
using Tables
using Plots

include(srcdir("SpatialSEIR_module.jl"))
using .SpatialSEIRutils


AUS_SA3_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA3_2021_AUST_SHP_GDA2020/SA3_2021_AUST_GDA2020.shp"))
GMelb_SA3_SHP = Tables.subset(AUS_SA3_SHP, AUS_SA3_SHP.GCC_NAME21 .== "Greater Melbourne")


MM = MixingMatrices.Spatial_MixingMatrix(GMelb_SA3_SHP)

mixmat = SpatialSEIR.MixingMatrix(MM, GMelb_SA3_SHP.SA3_NAME21, fill(1000, 40))

parameto = SpatialSEIR.Parameters(0.2, 0.1, 0.1, mixmat)

tot_log, inf_log = SpatialSEIR.run_sim(parameto, 1)

using StatsPlots
using DataFrames

@df tot_log Plots.plot(:t, [:S :E :I])

@df inf_log plot(:t, [cols(2:41)], legend = false)

frac_infs = (select(inf_log, Not(:t))./1000)|>Matrix
using ColorSchemes

anim = @animate for i in 1:1000:nrow(inf_log)
    colvec = Vector()
    for j in frac_infs[i, :]
        push!(colvec, ColorSchemes.viridis[j])
    end
    plot(GMelb_SA3_SHP.geometry, color = colvec')
end

gif(anim)