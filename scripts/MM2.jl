using Pkg
Pkg.activate(".")

using CSV
using DataFrames
include("../src/MixingMatrices.jl")


SA2_Pop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame) |> x -> 
    transform(x, :("SA2 code") => ByRow(string) => :("SA2 code"))

Melb_codes = expandCodes(SA2_Pop[:,1])|> x -> 
    leftjoin(x, SA2_Pop, on = (:SA2 => Symbol("SA2 code")))

SA2_pops = Pair.(Melb_codes.SA2, Melb_codes.Pop)|> x -> Dict(x)

SA3_pops =Melb_codes|> x -> 
    groupby(x, :SA3)|> x ->
    combine(x, :Pop => sum => :Pop)
SA3_pops = Pair.(SA3_pops.SA3, SA3_pops.Pop)|> x -> Dict(x)

SA4_pops = Melb_codes|> x -> 
    groupby(x, :SA4)|> x -> 
    combine(x, :Pop => sum => :Pop)
SA4_pops = Pair.(SA4_pops.SA4, SA4_pops.Pop)|> x -> Dict(x)

ST_pop = Melb_codes|> x -> 
    groupby(x, :"S/T")|> x -> 
    combine(x, :Pop => sum => :Pop)
ST_pop = Pair.(ST_pop.:"S/T", ST_pop.Pop)|> x -> Dict(x)
Melb_pops = [SA2_pops, SA3_pops, SA4_pops, ST_pop]


ξ = [1/4, 1/4, 1/4, 1/4]

function SAMMij(i, j, ξ, codes, popns)
    for l in eachindex(ξ)
        if codes[i, l] == codes[j, l]
            l == 1 && return ξ[l] 
            ll_pop = popns[1]|> x -> get(x, codes[j, 1], false) 
            ul_pop = get(popns[l], codes[i, l], false) - get(popns[l-1], codes[i, l-1], false)
            return ξ[l] * (ll_pop/ul_pop)
        end
    end
end

MM = [Mij(i, j, [1/4,1/4,1/4,1/4], Melb_codes, Melb_pops) for i in eachindex(SA2_Codes), j in eachindex(SA2_Codes)]


MM
sum(MM[1, :])

using StatsPlots
heatmap(log.(MM))