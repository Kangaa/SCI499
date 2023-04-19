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



MixingMatrix(SpatialMixingMatrix(codes, 0.5), Gmelb.SA2_NAME, Gmelb.popn).mm




colorder = [1:ncol(codes) fill(0,ncol(codes))]
for i in colorder[:, 1]
    colorder[i,2] = unique(codes[:, i]) |> length
end
sort(colorder; dims =2, rev =true)
