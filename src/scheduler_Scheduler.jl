MM_types = ["OD", "HMM", "HPMM"]
SA_scales = ["SA2", "SA3", "SA4"]
Intervention_types = ["None", "local", "travel", "total"]
nsims = 50

#combine permutations into string for running script
params = Base.product(MM_types, SA_scales, Intervention_types)

for param in params
    MM_type= param[1]
    SA_scale = param[2]
    Intervention_type = param[3]
    #run juliaScheduler.jl script
    run(`julia --project=/home/jbender/Documents/SCI499_new src/juliaScheduler.jl $MM_type $SA_scale $Intervention_type $nsims`)
end
