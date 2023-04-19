#non-hierarchical mixing matrix
using DrWatson 

@quickactivate :SCI499
using DataFramesMeta

Gmelb_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp")) |> x ->
    Tables.subset(x, x.GCC_NAME21 .== "Greater Melbourne")|>
    DataFrame|> x-> 
    rename(x,
        "SA2_NAME21" => "SA2_NAME", 
        "SA2_CODE21" => "SA2_CODE" )

Vic_pop = CSV.read(datadir("VicPop21.csv"), DataFrame) 

Gmelb_pop = Vic_pop|> x -> 
    transform(x, "SA2 code" => (y -> string.(y)) => "SA2_CODE") |> x->
    rename(x, "2021" => "popn")

Gmelb =leftjoin(Gmelb_SA2_SHP, Gmelb_pop[:,["SA2_CODE","popn" ]], on = "SA2_CODE")


codes = Gmelb |>
    names |>
    (x-> contains.(x, "CODE")) |>
    (x-> select(Gmelb, x))|>
    (x -> select(x, 1:3))



SMM = SpatialMixingMatrix(codes, 0.5)

heatmap(SMM, title = "Spatial Mixing Matrix", c = :viridis)


## flat mixing matrix(mixing proportional to population)
PPMM = fill(1.0, nrow(Gmelb), nrow(Gmelb))
for i in eachindex(Gmelb.popn)
    for j in eachindex(Gmelb.popn)
        PPMM[i,j] = Gmelb.popn[j]/sum(Gmelb.popn)
    end
end

using StatsPlots
using Animations

anim = @animate for μ in 0:0.01:1
    MMM = ((μ*SMM) + ((1-μ)*PPMM))
    heatmap(log10.(MMM), title = "Mixing Matrix", c = :viridis, clims = (-6,1))
end
gif(anim, "plots\\MixedMixingMatrix.gif", fps = 7)

μ = 0.9

MMM = ((μ*SMM) + ((1-μ)*PPMM))
heatmap(log10.(MMM), title = "Mixing Matrix", c = :viridis, clims = (-5,1))
