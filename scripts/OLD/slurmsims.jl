## Activate project and get src
using Pkg
Pkg.activate("../SCI499")
Pkg.instantiate()
include("src/SCI499.jl")
using .SCI499


## Get cmd line params 
args  = SCI499.parse_commandline()
mm_parameter, pop_file_dir, gamma, beta = args
#create params object
using CSV, DataFrames
mm_parameter = 0.3
VicPop = CSV.read(pop_file_dir[2], DataFrame)
codes = VicPop[:,1]|> x -> string.(x)
popns = VicPop[:,3] .+ 1
mixmat = MixingMatrices.SpatialMixingMatrix(codes, popns, 0.5, 1)

params = (
    patch_names =  string.(VicPop[:,1]),
    population_per_patch = VicPop[:,3] .+ 1,
    mixing_matrix = mixmat,
    β = beta[2],
    γ = gamma[2]
)

#Add (40) worker  processes an send them the params
using Distributed

addprocs(40, exeflags=["--project"])

@everywhere global params = $params 
@everywhere global mm_parameter = 0.3
@everywhere global beta = 0.5
@everywhere global gamma = 0.1
## Batch sims function (move to src)
@everywhere begin
    include("src/SCI499.jl")
    using .SCI499
    using CSV
    function batch_sim(sim)
        tot_data, summary_stats, patch_sus, patch_inf = CompartmentalModels.simulate(params, 10, 1,false, sim)
        
        CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
        CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_patchinf$sim.csv" , patch_inf)
        CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_patchsus$sim.csv" , patch_sus)
        CSV.write("data/SA2/SA2_$(beta)_$(gamma)_$(mm_parameter)_summary$sim.csv" , summary_stats)
    end
end


#run sims
nsims = 1000

pmap(batch_sim, 1:nsims)

rmprocs(workers())

#save results
