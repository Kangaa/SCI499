## Slurm cluster job test
using Distributed, ClusterManagers

#request some tasks
addprocs(SlurmManager(16), exeflags=["--project", "threads=2"])

@everywhere function get_info()
    results = Dict{Symbol, Any}()
    results[:num_threads] = Threads.nthreads()
    results[:id] = myid()
    results[:hostname] = gethostname()
    return results
end

info = []
info_tasks = []
for i in workers()
    worker_info = fetch(@spawnat i get_info())
    push!(info, worker_info)
    println(worker_info)
end

rmprocs(workers())

