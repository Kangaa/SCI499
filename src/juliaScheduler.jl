#Parse Args
MM_type = ARGS[1]
SA_scale = ARGS[2]
nsims = parse(Int64, ARGS[3])

# Define the range of parameters to use

if MM_type == "HMM"
    β = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    γ = [1.0]
    params = Base.product(β, γ)
    for param in params
        beta = param[1]
        gamma = param[2]
        # Define the Slurm batch job parameters
        jobname = "SIR_$(MM_type)_$(SA_scale)__$(beta)_$(gamma)"
        nodes = 1
        tasks_per_node = 1
        cpus_per_task = 40
        mem = 1000
        time = "20:00:00"
        output = "/home/jbender/Documents/SCI499_new/slurmout/%j.out"
    
        # Define the command to run
        command = "julia --project=/home/jbender/Documents/SCI499_new src/DistributedSims.jl $MM_type $beta $gamma $SA_scale $nsims"
    
        # Generate the Slurm batch job script
        script = "#!/bin/bash
    #SBATCH --job-name=$jobname
    #SBATCH --nodes=$nodes
    #SBATCH --ntasks-per-node=$tasks_per_node
    #SBATCH --cpus-per-task=$cpus_per_task
    #SBATCH --mem=$mem
    #SBATCH --time=$time
    #SBATCH --output=$output
    
    $command
    "
    
        # Write the script to a file
        filename = "slurmjobs/slurmjob_$(SA_scale)_$(gamma).sh"
        open(filename, "w") do f
            write(f, script)
        end
    
        # Submit the job
        run(`sbatch $filename`)
    end

elseif MM_type == "HPMM"
    mmcoef = [0.0, 0.5, 1.0]
    β = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    γ = [1.0]
    params = Base.product(mmcoef, β, γ)
    for param in params
        mm = param[1]
        beta = param[2]
        gamma = param[3]
        # Define the Slurm batch job parameters
        jobname = "SIR_$(MM_type)_$(SA_scale)_$(mm)_$(beta)_$(gamma)"
        nodes = 1
        tasks_per_node = 1
        cpus_per_task = 40
        mem = 1000
        time = "20:00:00"
        output = "/home/jbender/Documents/SCI499_new/slurmout/%j.out"
    
        # Define the command to run
        command = "julia --project=/home/jbender/Documents/SCI499_new src/DistributedSims.jl $MM_type $mm $beta $gamma $SA_scale $nsims"
    
        # Generate the Slurm batch job script
        script = "#!/bin/bash
    #SBATCH --job-name=$jobname
    #SBATCH --nodes=$nodes
    #SBATCH --ntasks-per-node=$tasks_per_node
    #SBATCH --cpus-per-task=$cpus_per_task
    #SBATCH --mem=$mem
    #SBATCH --time=$time
    #SBATCH --output=$output
    
    $command
    "
    
        # Write the script to a file
        filename = "slurmjobs/slurmjob_$(SA_scale)_$(mm)_$(beta)_$(gamma).sh"
        open(filename, "w") do f
            write(f, script)
        end
    
        # Submit the job
        run(`sbatch $filename`)
    end

    
elseif MM_type == "OD"
    deltaH = [0.1, 0.2, 0.4, 0.6, 0.8]
    β = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    γ = [1.0]
    params = Base.product(deltaH, β, γ)
    for param in params
        mm = param[1]
        beta = param[2]
        gamma = param[3]
        # Define the Slurm batch job parameters
        jobname = "SIR_$(MM_type)_$(SA_scale)_$(mm)_$(beta)_$(gamma)"
        nodes = 1
        tasks_per_node = 1
        cpus_per_task = 40
        mem = 1000
        time = "20:00:00"
        output = "/home/jbender/Documents/SCI499_new/slurmout/%j.out"
    
        # Define the command to run
        command = "julia --project=/home/jbender/Documents/SCI499_new src/DistributedSims.jl $MM_type $mm $beta $gamma $SA_scale $nsims"
    
        # Generate the Slurm batch job script
        script = "#!/bin/bash
    #SBATCH --job-name=$jobname
    #SBATCH --nodes=$nodes
    #SBATCH --ntasks-per-node=$tasks_per_node
    #SBATCH --cpus-per-task=$cpus_per_task
    #SBATCH --mem=$mem
    #SBATCH --time=$time
    #SBATCH --output=$output
    
    $command
    "
    
        # Write the script to a file
        filename = "slurmjobs/slurmjob_$(SA_scale)_$(mm)_$(beta)_$(gamma).sh"
        open(filename, "w") do f
            write(f, script)
        end
    
        # Submit the job
        run(`sbatch $filename`)
    end
end


# Loop over the parameters and generate a Slurm batch job for each one
