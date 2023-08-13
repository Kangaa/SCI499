## Activate project and get src
using Pkg
Pkg.activate("SCI499")
Pkg.instantiate()
include("src/SCI499.jl")
using .SCI499

## load packages
using Distributed
using CSV
using DataFrames
using ProfileView
using JET

function Distrosim(mm_parameter, beta, gamma, SA_scale, nsims)
    VicPop = CSV.read("data/Gmelb$(SA_scale)Pop21.csv", DataFrame)
    codes::Vector{String} = VicPop[:, 1]::Vector{Int64} |> x -> string.(x)
    names::Vector{String} = VicPop[:, 2] |> x -> string.(x)
    popns::Vector{Int64}  = VicPop[:, 3] .+ 1
    CodePops = DataFrame(
        Codes = codes,
        Pop = popns)
    SA_scale == "SA2" ?  ξ = [1/4,1/4,1/4,1/4] : 
    SA_scale == "SA3" ?  ξ = [1/2,1/4,1/4] :
    SA_scale == "SA4" ?  ξ = [3/4,1/4] : error("SA_Scale must be SA2, SA3 or SA4")
    mixmat = SCI499.MixingMatrices.SpatialMixingMatrix(CodePops, ξ, mm_parameter)
      sim_params = (
            patch_names = names::Vector{String},
            population_per_patch = popns::Vector{Int64},
            mixing_matrix = mixmat::Matrix{Float64},
            β = beta::Float64,
            γ = gamma::Float64,
        )::Tuple{Vector{String}, Vector{Int64}, Matrix{Float64}, Float64, Float64}
    #run sims
 map(function batch_sim(sim) 
        tot_data, summary_stats, patch_sus, patch_inf = CompartmentalModels.simulate(sim_params, 10, 1,false, sim)

            CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
            CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchinf_$sim.csv" , patch_inf)
            CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchsus_$sim.csv" , patch_sus)
            CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_summary_$sim.csv" , summary_stats)
        end,
    1:nsims)
end
using Cthulhu
using JET

@report_opt Distrosim(0.5, 1.5, 1.0, "SA2", 1)