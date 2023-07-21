## Activate project and get src
using Pkg
Pkg.activate("../SCI499")
Pkg.instantiate()
include("src/SCI499.jl")
using .SCI499


#Add (40) worker  processes an send them the params
using Distributed
using ClusterManagers
using CSV
using DataFrames
## Batch sims function (move to src)


for mm_parameter in [0.1, 0.3, 0.5, 0.7, 0.9, 1]
    addprocs(SlurmManager(2), t = "00:5:00")

    @everywhere begin
    mm_parameter = $mm_parameter
    my_print = function(x)
        print(mm_parameter)
    end
end
        #run sims
        nsims = 100
        pmap(my_print,
        1:nsims)

    rmprocs(workers())
end

