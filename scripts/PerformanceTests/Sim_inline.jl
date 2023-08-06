## Activate project and get src
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../../src/SCI499.jl")
using .SCI499

## load packages
using Distributed
using CSV
using DataFrames
using ProfileView
using BenchmarkTools

mm_parameter = 0.5
beta = 1.5
gamma = 1.0
SA_scale =  "SA2"
nsims = 1

VicPop = CSV.read("data/Gmelb$(SA_scale)Pop21.csv", DataFrame)
patch_names = convert(Vector{String}, VicPop[:,2]) 
codes = VicPop[:, 1]|> x -> string.(x)
popns  = VicPop[:, 3] .+ 1

CodePops = DataFrame(
    Codes = codes,
    Pop = popns)

SA_scale == "SA2" ?  ξ = [1/4,1/4,1/4,1/4] : 
SA_scale == "SA3" ?  ξ = [1/2,1/4,1/4] :
SA_scale == "SA4" ?  ξ = [3/4,1/4] : error("SA_Scale must be SA2, SA3 or SA4")

mixmat = SCI499.MixingMatrices.SpatialMixingMatrix(CodePops, ξ, mm_parameter)

population_per_patch = popns
mixing_matrix = mixmat
β = beta
γ = gamma 

const  params = (
    patcn

        population_per_patch = popns,
        mixing_matrix = mixmat,
        β = beta,
        γ = gamma
    )

struct CompartmentalModel
    state::NamedTuple{(:S, :I), Tuple{Vector{Int64}, Vector{Int64}}}
    transition_parameters::NamedTuple{(:β, :γ), Tuple{Vector{Float64}, Vector{Float64}}}
    total_population::Int64
    num_patches::Int64
    patch_names::Vector{String}
    population_per_patch::Vector{Int64}
    mixing_matrix::Matrix{Float64}
end



population_per_patch = popns
    num_patches = length(patch_names)
    state = (
        S = Vector{Int64}(popns),
        I = Vector{Int64}(zeros(Int, num_patches)))
    transition_parameters = (
        β = Vector{Float64}(fill(beta, num_patches)),
        γ = Vector{Float64}(fill(gamma, num_patches)))
    total_population = sum(popns)
    num_patches = length(patch_names)
    patch_names = Vector{String}(patch_names)
    population_per_patch = Vector{Int64}(population_per_patch)
    mixing_matrix = mixing_matrix
model =  CompartmentalModel(state, transition_parameters, total_population, num_patches, patch_names, population_per_patch, mixing_matrix)
    

function rand_infection!(model::CompartmentalModel, num_infected::Int64)::CompartmentalModel
    # Randomly select a patch to seed
    patch::Int64 = rand(1:model.num_patches)
    # Set the number of infected individuals in the patch
    model.state.I[patch]::Int64 += num_infected
    # Set the number of susceptible individuals in the patch
    model.state.S[patch]::Int64 -= num_infected
    return model
end

mutable struct EventRates
    rates::Vector{Float64} ## is this the best way to store rates (or as perEvent type)
    net::Float64

    function EventRates(model::CompartmentalModel)
        L = length(model.state)*model.num_patches
        rates = Vector{Float64}(zeros(L))
        net = 0.0
        return new(rates)
    end
end

function gillespie_pick(rates)::Int
    crates = cumsum(rates)
    r = rand()
    ix = findfirst(x -> x > r, crates)
    ix -1
end

## Update the event rates based on state

function update_rates!(model::CompartmentalModel, rates::EventRates)
    ## fraction of susceptible in each patch
    @inbounds s_frac::Vector{Float64} =  model.state.S./model.population_per_patch
    ##force of infection
    @inbounds f::Vector{Float64} = model.transition_parameters.β[1]*(transpose(model.mixing_matrix)*model.state.I)
    ## Infection rate
    @inbounds rates.rates[1:model.num_patches]::Vector{Float64} = s_frac.*f
    ## Recovery rate
    @inbounds rates.rates[model.num_patches+1:2model.num_patches]::Vector{Float64} = model.state.I.*(model.transition_parameters.γ[1])
    @inbounds rates.net = sum(rates.rates)

    return rates
end


## update the model state based on the event rates
function update_state!(model::CompartmentalModel, event_ix, return_event = true)
    ###pick an event to occur using gillespie algorithm
    event_type = (event_ix/model.num_patches)|> floor
    event_location = (event_ix%model.num_patches) +1
    if event_type == 0
        @inbounds    model.state.S[event_location] -= 1
        @inbounds    model.state.I[event_location] += 1
    elseif event_type == 1
        @inbounds   model.state.I[event_location] -= 1
    end
    if return_event
        return  event_type, event_location
    end
end

function exp_sample(net_rate)::Float64
    -(rand() |> log)/net_rate
end

rates = EventRates(model)
t = 0.0
wk::Int64 = 0
update_rates!(model, rates)

## Create log and data frames
if include_log
    sim_log = DataFrame(time = Float64[], event_type = Int[], event_location = Int[])
end

tot_data = DataFrame(
    time = Int64,
    TotalSusceptible = Int64,
    TotalInfected = Int64)

patch_inf =  DataFrame(Dict(name => Vector{Int}() for name in params[:patch_names]))
patch_sus =  DataFrame(Dict(name => Vector{Int}() for name in params[:patch_names]))
    
summary_stats = DataFrame(
    total_infections = 0,
    peak_infections = 0,
    duration = 0)
#run simulation loop

while rates.net != 0.0
    update_rates!(model, rates)
    if rates.net != 0.0
        t += exp_sample(rates.net)
        event_ix = gillespie_pick(rates.rates./rates.net)
        event_type, event_location = update_state!(model, event_ix)

        if event_type == 0
            summary_stats[1,1] += 1
        end 
        if sum(model.state.I) > summary_stats[1, 2]
           summary_stats[1, 2] = sum(model.state.I) 
        end

        if include_log
            push!(sim_log, (time = t, event_type = event_type, event_location = event_location))
        end

        if t/aggregation_time > wk
            wk += 1
            push!(tot_data, (time=wk, TotalSusceptible=sum(model.state.S), TotalInfected=sum(model.state.I)), promote = true)
            push!(patch_sus, model.state.S)
            push!(patch_inf, model.state.I)
        end
    end

end

summary_stats[1, 3] = wk

include_log && CSV.write(datadir("$sim log.csv") , sim_log, append=false, header=[:time, :event_type, :event_location])


CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_$sim.csv" , tot_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchinf_$sim.csv" , patch_inf)
CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_patchsus_$sim.csv" , patch_sus)
CSV.write("data/sims/test/$(SA_scale)/$(SA_scale)_$(beta)_$(gamma)_$(mm_parameter)_summary_$sim.csv" , summary_stats)
