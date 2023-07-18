## Activate project and get src
using Pkg
Pkg.activate("../SCI499")
Pkg.instantiate()
include("src/SCI499.jl")
using .SCI499


#Add (40) worker  processes an send them the params
using Distributed
addprocs(40, exeflags=["--project"])

## Batch sims function (move to src)
@everywhere begin
    include("src/SCI499.jl")
    using .SCI499
    using CSV
    using DataFrames
    R0 = [  (β = 0.1, γ = 0.5),
    (β = 0.5, γ = 0.5),
    (β = 1.0, γ = 0.5),
    (β = 0.5, γ = 1.0),
    (β = 0.1, γ = 1.0),
    (β = 1.0, γ = 0.1),
    (β = 0.5, γ = 0.1)]

mm = [0.3, 0.5, 0.7, 0.9, 1]

for mm_parameter in mm

    VicPop = CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
    codes = VicPop[:,1]|> x -> string.(x)
    popns = VicPop[:,3] .+ 1    
    mixmat = MixingMatrices.SpatialMixingMatrix(codes, popns, 0.5, mm_parameter)

    for R0 in R0
        beta = R0[:β]
        gamma = R0[:γ]

        sim_params = (
            patch_names =  string.(VicPop[:,1]),
            population_per_patch = popns,
            mixing_matrix = mixmat,
            β = beta,
            γ = gamma
        )
    
        #run sims
        nsims = 100
        pmap(function batch_sim(sim)
            tot_data, summary_stats, patch_sus, patch_inf = CompartmentalModels.simulate(sim_params, 10, 1,false, sim)

            CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
            CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_patchinf_$sim.csv" , patch_inf)
            CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_patchsus_$sim.csv" , patch_sus)
            CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_summary_$sim.csv" , summary_stats)
        end,
        1:nsims)
    end
end

end




rmprocs(workers())

#save results
