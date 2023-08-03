using DrWatson
@quickactivate 

using Shapefile
using Tables
using Plots

include(srcdir("SpatialSEIR_module.jl"))
using .SpatialSEIRutils
using DataFramesMeta
using StatsPlots
using DataFrames

AUS_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp"))
GMelb_SA2_SHP = Tables.subset(AUS_SA2_SHP, AUS_SA2_SHP.GCC_NAME21 .== "Greater Melbourne")

size(MM,1)

MM = MixingMatrices.Spatial_MixingMatrix2(GMelb_SA2_SHP, 0.5)

MM = MixingMatrices.Spatial_MixingMatrix(GMelb_SA2_SHP, 0.5)

mixmat = SpatialSEIR.MixingMatrix(MM, GMelb_SA2_SHP.SA2_NAME21, fill(10, size(MM,1)))

parameto = SpatialSEIR.Parameters(1, 0.1, 0.1, mixmat)

tot_log, inf_log = SpatialSEIR.run_sim(parameto, 1)



@df tot_log Plots.plot(:t, [:S :E :I])

@df inf_log plot(:t, [cols(2:41)], legend = false)

frac_infs = (select(inf_log, Not(:t))./10)|>Matrix
using ColorSchemes

anim = @animate for i in 1:200:nrow(inf_log)
    colvec = Vector()
    for j in frac_infs[i, :]
        push!(colvec, ColorSchemes.viridis[j])
    end
    plot(GMelb_SA2_SHP.geometry, color = colvec')
end

gif(anim)