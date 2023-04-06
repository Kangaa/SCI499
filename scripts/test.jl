using DrWatson
@quickactivate

include(srcdir("SpatialSEIR_module.jl"))
using .SpatialSEIRutils.MixingMatrices: getCodes
using Shapefile
using Tables
using DataFrames

AUS_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp"))
Shp = Tables.subset(AUS_SA2_SHP, AUS_SA2_SHP.GCC_NAME21 .== "Greater Melbourne")   


Codes = getCodes(Shp)
npatch = size(Codes, 1)
intraregioncoefficient = 0.8
interregioncoefficient = 1.0 - intraregioncoefficient
MM = [prod(
        if Codes[i, l] == Codes[j, l]
            nLLinUL = sum(Codes[i, l] .== Codes[:, l])
            intraregioncoefficient / nLLinUL
        elseif l == ncol(Codes) - 1
            nLLoutUL = sum(Codes[i, l] .!= Codes[:, l])
            interregioncoefficient / nLLoutUL
        else
            interregioncoefficient
        end
        for l in 1:ncol(Codes) - 1
    ) for i in 1:npatch, j in 1:npatch]
return MM
