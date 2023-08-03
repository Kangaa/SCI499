using DrWatson
@quickactivate

include(srcdir("SpatialSEIR_module.jl"))
using .SpatialSEIRutils.MixingMatrices: getCodes
using Shapefile
using Tables
using DataFrames

AUS_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp"))
Shp = Tables.subset(AUS_SA2_SHP, AUS_SA2_SHP.GCC_NAME21 .== "Greater Melbourne")   
Codes = Shp |> getCodes
npatch = nrow(Codes)

MM = fill(1.0, npatch, npatch)

intraregioncoefficient = 0.1
interregioncoefficient = 1.0-intraregioncoefficient

for (i, SAo) in enumerate(Codes[:,1])
    for (j, SAd) in enumerate(Codes[:, 1])
        for l in 1:(ncol(Codes)-1)
            if Codes[i,l] == Codes[j,l]
                nLLinUL = (Codes[i,l] .== Codes[:, l])|> sum
                coef = (intraregioncoefficient/nLLinUL)
                MM[i,j] *= coef
                break
            else
                if l == (ncol(Codes)-1) 
                    nLLoutUL = (Codes[i,l] .!= Codes[:, l])|> sum
                    MM[i,j] *= interregioncoefficient/nLLoutUL
                else
                    coef = interregioncoefficient
                    MM[i,j] *= coef
                end
            end
        end
    end
end






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

MM[9,:]|>sum

log.(MM)
test = Spatial_MixingMatrix2(Shp, 0.9)
using StatsPlots
heatmap(log.(MM))

      α₁ =0.3
      β₁ = 0.7
      
      α₂ =0.3
      β₂ = 0.7

      α₃ =0.3
      β₃ = 0.7


      
      (α₁) + (β₁/4) + (β₁/4) + (β₁/4) + ( β₁/4)

      (α₁) + (β₁*α₂/2) + (β₁*α₂/2) + (β₁*β₂/2) + ( β₁*β₂)/2

      (α₁) + (β₁*α₂/2) + (β₁*α₂/2) + (β₁*β₂*α₃/5)+ (β₁*β₂*α₃/5) + (β₁*β₂*α₃/5) + (β₁*β₂*α₃/5) +(β₁*β₂*α₃/5) + (β₁*β₂*β₃)

