using DrWatson
using DataFrames
using CSV
using LinearAlgebra
# Create a CompartmentalModel type to hold the state and parameters of the model
struct CompartmentalModel
    state::NamedTuple{(:S, :I), Tuple{Vector{Int64}, Vector{Int64}}}
    transition_parameters::NamedTuple{(:β, :γ), Tuple{Vector{Float64}, Vector{Float64}}}
    total_population::Int64
    num_patches::Int64
    patch_names::Vector{String}
    population_per_patch::Vector{Int64}
    mixing_matrix::AbstractArray{Float64}
    foi::Vector{Float64}
    s_frac::Vector{Float64}
    intervention::Vector{Bool}
end

function SIR(patch_names::Vector{String}, population_per_patch::Vector{Int64}, mixing_matrix::Transpose{Float64, Matrix{Float64}}, β::Float64, γ::Float64)
    num_patches = length(patch_names)
    state = (
    S = Vector{Int64}(population_per_patch),
    I = Vector{Int64}(zeros(Int, num_patches)))
    transition_parameters = (
    β = Vector{Float64}(fill(β, num_patches)),
    γ = Vector{Float64}(fill(γ, num_patches)))
    total_population = sum(population_per_patch)
    num_patches = length(patch_names)
    patch_names = Vector{String}(patch_names)
    population_per_patch = Vector{Int64}(population_per_patch)
    mixing_matrix = mixing_matrix
    foi = Vector{Float64}(zeros(num_patches))
    s_frac = Vector{Float64}(ones(num_patches))
    intervention = Vector{Bool}(zeros(num_patches))
    return CompartmentalModel(state, transition_parameters, total_population, num_patches, patch_names, population_per_patch, mixing_matrix, foi, s_frac, intervention)
end

## Function to seed the model with a number of infected individuals at a random patch
function rand_infection!(model::CompartmentalModel, num_infected::Int64)::CompartmentalModel
    # Randomly select a patch to seed
    patch::Int64 = rand(1:model.num_patches)
    # Set the number of infected individuals in the patch
    model.state.I[patch]::Int64 += num_infected
    # Set the number of susceptible individuals in the patch
    model.state.S[patch]::Int64 -= num_infected
    return model
end

## Add a type to hold the event rates
mutable struct EventRates
    rates::Vector{Float64} ## is this the best way to store rates (or as perEvent type)
    net::Float64
    crates::Vector{Float64}
    
    function EventRates(model::CompartmentalModel)
        L = length(model.state)*model.num_patches
        return new(
        Vector{Float64}(zeros(L)),
        0.0,
        Vector{Float64}(zeros(L)))
    end
end

function gillespie_pick(rates)
    rates.crates .= rates.rates./rates.net
    cumsum!(rates.crates, rates.crates)
    r = rand()
    findfirst(x -> x > r, rates.crates) -  1
end

## Update the event rates based on state

function update_rates!(model::CompartmentalModel, rates::EventRates)
    ##force of infection
    mul!(model.foi::Vector{Float64}, model.mixing_matrix::Transpose{Float64, Matrix{Float64}}, model.state.I::Vector{Int64}, model.transition_parameters.β[1]::Float64, 0)
    @inbounds rates.rates[1:model.num_patches] .= (model.s_frac .* model.foi)
    ## Recovery rate
    @inbounds rates.rates[(model.num_patches+1):(2*model.num_patches)] .= (model.state.I .* model.transition_parameters.γ[1])
    @inbounds rates.net = sum(rates.rates)
    
    return rates
end


mutable struct event
    index::Int64
    type::Int64
    location::Int64
end

## update the model state based on the event rates
function update_state!(model::CompartmentalModel, event::event, rates::EventRates,  return_event = true)
    ###pick an event to occur using gillespie algorithm
    event.index = gillespie_pick(rates)::Int64
    event.type = event.index ÷ model.num_patches
    event.location = (event.index%model.num_patches) +1
    if event.type == 0
        @inbounds    model.state.S[event.location] -= 1
        @inbounds    model.state.I[event.location] += 1
    elseif event.type == 1
        @inbounds   model.state.I[event.location] -= 1
    end
    if return_event
        return  event.type, event.location
    end
end

##  update mixmat if 

function exp_sample(net_rate)::Float64
    -(rand() |> log)/net_rate
end

#update mixing matrix to simulate intervention
function intervention!(mixmat, intervention_type, patch)
    if intervention_type == "none"
        return mixmat

    elseif intervention_type == "local"
        mixmat[patch,patch] *= 0.5

    elseif intervention_type == "travel"
        for i in 1:size(mixmat, 2)
            if i != patch
                mixmat[patch,i] *= 0.5
                mixmat[i,patch] *= 0.5
            end
        end

    elseif intervention_type == "total"
        for i in 1:size(mixmat, 2)
            if i != patch
                mixmat[patch,i] *= 0.5
                mixmat[i,patch] *= 0.5
            end
            mixmat[patch,patch] *= 0.5
        end
    else
        error("Intervention type must be none, local, travel or total")
    end
end 


## Define a function to simulate a new model

function sim_loop(model::CompartmentalModel, rates::EventRates, t::Float64, wk::Int64, aggregation_time::Int64, include_log::Bool, Intervention_type,  tot_data, patch_inf, patch_sus)
    event_log = event(0,0,0)
    while rates.net != 0.0
        update_rates!(model, rates)
        if rates.net != 0.0
            t += exp_sample(rates.net)
            event_type, event_location = update_state!(model, event_log, rates)
            #Calculate new s_frac vector
            if event_type == 0 
                model.s_frac[event_location] =  (model.state.S[event_location] / model.population_per_patch[event_location])
            end 
            ## if s_frac has changed to >x% then update mixmat 
            
            if model.intervention[event_location] == 0 && model.s_frac[event_location] > 0.2
                intervention!(model.mixing_matrix, 
                Intervention_type, event_location)
                model.intervention[event_location] = 1
            end
            
            
            if t/aggregation_time > wk
                wk += 1
                push!(tot_data, [wk, sum(model.state.S), sum(model.state.I)])
                push!(patch_sus, model.state.S)
                push!(patch_inf, model.state.I)
            end
        end
    end
end

function simulate(params, i0::Int64, aggregation_time::Int64, intervention_type::String,  include_log::Bool, sim::Int64)
    #setup 
    model = SIR(
    params.patch_names,
    params.population_per_patch,
    transpose(params.mixing_matrix), 
    params.β,
    params.γ)
    intervention_type = intervention_type
    model |> x -> rand_infection!(x, i0)
    
    rates = EventRates(model)
    
    t = 0.0
    
    wk = 0
    
    update_rates!(model, rates)
    
    ## Create log and data frames
    if include_log
        sim_log = DataFrame(time = Float64[], event_type = Int[], event_location = Int[])
    end
    
    tot_data = DataFrame(
    wk = Vector{Int64}(),
    TotalSusceptible = Vector{Int64}(),
    TotalInfected = Vector{Int64}())
    
    patch_inf =  DataFrame([Vector{Int64}() for i in 1:model.num_patches], model.patch_names)
    patch_sus =  DataFrame([Vector{Int64}() for i in 1:model.num_patches], model.patch_names)
    
    #run simulation loop
    sim_loop(model, rates, t, wk, aggregation_time, include_log, intervention_type,  tot_data, patch_inf, patch_sus)
    
    include_log && CSV.write(datadir("$sim log.csv") , sim_log, append=false, header=[:time, :event_type, :event_location])
    return (tot_data, patch_sus, patch_inf)
end