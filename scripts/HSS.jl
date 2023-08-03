using DrWatson
@quickactivate :SCI499
using CSV
using DataFrames
using Shapefile

shps = Shapefile.Table("data\\ASGS_GDA2020\\SA2_2021_AUST_SHP_GDA2020\\SA2_2021_AUST_GDA2020.shp")

codes = shps.SA2_CODE21

VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
codes = VicPop[:,1]|> x -> string.(x)
popns = VicPop[:,3] .+ 1  
names = string.(VicPop[:,1])
mm_parameter = 0.5

CodePops = DataFrame(Codes = codes, Pop = popns)
mixmat = SCI499.MixingMatrices.SpatialMixingMatrix(CodePops, [1/4,1/4,1/4,1/4], 0.5)

using StatsPlots
heatmap(mixmat)

using HssMatrices
mixmat_hss = hss(mixmat)
plotranks(mixmat_hss)


HssMatrices._hssleaf(mixmat_hss)
using LinearAlgebra

norm(mixmat - mixmat_hss)/norm(mixmat)

using BenchmarkTools

x = randn(size(mixmat, 2), 1)

@btime mixmat*x
@btime mixmat_hss*x

generators(mixmat_hss, (1,2))

@hierarchical
