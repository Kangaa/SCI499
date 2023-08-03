#function main()
    # Parse the command line arguments
    #args = parse_commandline()
    args = Dict("mixing_matrix_parameter" => 1.0,
     "beta" => 0.5,
     "gamma" => 0.1,
     "pop_info_file" => "data/vic_pop2021.csv")


     #mixing matrix is loaded once per parameter value
     mixmat = SpatialMixingMatrix(args["mixing_matrix_parameter"], args["pop_info_file"])
     #create a run sims function which uses Distributed.jl and ClusterManagers.jl to run the simulations in parallel
        function run_sims(args, mixmat)
            #create a cluster manager
            addprocs(4)
            #create a vector of processes
            procs = workers()
            #create a vector of arguments to pass to each process
            args_list = [args for i in 1:length(procs)]
            #run the simulations in parallel
            @sync @distributed for i in 1:length(procs)
                #run the simulation
                run_sim(args_list[i], mixmat)
            end
        end


#end
#main()