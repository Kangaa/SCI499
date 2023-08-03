## Activate project and get src
using Pkg
Pkg.activate("../SCI499")
Pkg.instantiate()
include("src/SCI499.jl")
using .SCI499

using Distributed
using CSV
using DataFrames
## Batch sims function (move to src)

VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
codes = VicPop[:,1]|> x -> string.(x)
popns = VicPop[:,3] .+ 1  
names = string.(VicPop[:,1])
mm_parameter = parse(Float64, ARGS[1])

addprocs(40,   exeflags="--project")

@everywhere begin
    include("src/SCI499.jl")
    using .SCI499
    using CSV
    using DataFrames
    popns = $popns
    codes = $codes 
    names = $names   
    mixmat = SCI499.MixingMatrices.SpatialMixingMatrix(codes, popns, 0.5, $mm_parameter)

    R0 = [(β = 0.1, γ = 0.5),
        (β = 0.5, γ = 0.5),
        (β = 1.0, γ = 0.5),
        (β = 0.5, γ = 1.0),
        (β = 0.1, γ = 1.0),
        (β = 1.0, γ = 0.1),
        (β = 0.5, γ = 0.1)]


    for R0 in R0
        beta = R0[:β]
        gamma = R0[:γ]

        sim_params = (
            patch_names = names  ,
            population_per_patch = popns,
            mixing_matrix = mixmat,
            β = beta,
            γ = gamma
        )

        #run sims
        nsims = 100
        pmap(function batch_sim(sim)
            tot_data, summary_stats, patch_sus, patch_inf = CompartmentalModels.simulate(sim_params, 10, 1,false, sim)

            CSV.write("data/SA2_julia/SA2_$(beta)_$(gamma)_$(mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
            CSV.write("data/SA2_julia/SA2_$(beta)_$(gamma)_$(mm_parameter)_patchinf_$sim.csv" , patch_inf)
            CSV.write("data/SA2_julia/SA2_$(beta)_$(gamma)_$(mm_parameter)_patchsus_$sim.csv" , patch_sus)
            CSV.write("data/SA2_julia/SA2_$(beta)_$(gamma)_$(mm_parameter)_summary_$sim.csv" , summary_stats)
        end,
        1:nsims)
    end
end
rmprocs(workers())

