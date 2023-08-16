using DrWatson
using DataFrames
using CSV
# Create a CompartmentalModel type to hold the state and parameters of the model
struct CompartmentalModel
    state::NamedTuple{(:S, :I), Tuple{Vector{Int64}, Vector{Int64}}}
    transition_parameters::NamedTuple{(:β, :γ), Tuple{Vector{Float64}, Vector{Float64}}}
    total_population::Int64
    num_patches::Int64
    patch_names::Vector{String}
    population_per_patch::Vector{Int64}
    mixing_matrix::Matrix{Float64}
    foi::Vector{Float64}
end

function SIR(patch_names::Vector{String}, population_per_patch::Vector{Int64}, mixing_matrix::Matrix{Float64}, β::Float64, γ::Float64)
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
    return CompartmentalModel(state, transition_parameters, total_population, num_patches, patch_names, population_per_patch, mixing_matrix)
end

## Function to seed the model with a number of infected individuals at a random patch
function rand_infection!(model::CompartmentalModel, num_infected::Int64)::CompartmentalModel
    # Randomly select a patch to seed
    patch::Int64 = rand(1:model.num_patches)
    # Set the number of infected individuals in the patch
    model.state.I[patch]::Int64 .+= num_infected
    # Set the number of susceptible individuals in the patch
    model.state.S[patch]::Int64 .-= num_infected
    return model
end

## Add a type to hold the event rates
mutable struct EventRates
    rates::Vector{Float64} ## is this the best way to store rates (or as perEvent type)
    net::Float64
    crates::Vector{Float64}

    function EventRates(model::CompartmentalModel)
        L = length(model.state)*model.num_patches
        rates = Vector{Float64}(zeros(L))
        net = 0.0
        return new(rates)
    end
end

function gillespie_pick(rates)::Int
    rates.crates .= cumsum(rates)
    r = rand()
    ix = findfirst(x -> x > r, crates)
    ix - 1
end

## Update the event rates based on state

function update_rates!(model::CompartmentalModel, rates::EventRates)
    ## fraction of susceptible in each patch
    @inbounds s_frac::Vector{Float64} =  model.state.S./model.population_per_patch
    ##force of infection
    @inbounds model.foi::Vector{Float64} .= model.transition_parameters.β[1]*(transpose(model.mixing_matrix)*model.state.I)
    ## Infection rate
    @inbounds rates.rates[1:model.num_patches]::Vector{Float64} .= s_frac .* f
    ## Recovery rate
    @inbounds rates.rates[model.num_patches+1:2model.num_patches]::Vector{Float64} .= model.state.I .* model.transition_parameters.γ[1]
    @inbounds rates.net = sum(rates.rates)

    return rates
end


struct event
    index::Int64
    type::Int64
    location::Int64
end

## update the model state based on the event rates
function update_state!(model::CompartmentalModel, event::event, rates::EventRates,  return_event = true)
    ###pick an event to occur using gillespie algorithm
    event.index = gillespie_pick(rates.rates./rates.net)
    event.type = (event.index/model.num_patches)|> floor
    event.location = (event.index%model.num_patches) +1
    if event_type == 0
        @inbounds    model.state.S[event.location] -= 1
        @inbounds    model.state.I[event.location] += 1
    elseif event_type == 1
        @inbounds   model.state.I[event.location] -= 1
    end
    if return_event
        return  event.type, event.location
    end
end

function exp_sample(net_rate)::Float64
    -(rand() |> log)/net_rate
end
## Define a function to simulate a new model

function sim_loop(model::CompartmentalModel, rates::EventRates, t::Float64, wk::Int64, aggregation_time::Int64, include_log::Bool, tot_data, patch_inf, patch_sus)
    event = event(0,0,0)
    while rates.net != 0.0
        update_rates!(model, rates)
        if rates.net != 0.0
            t .+= exp_sample(rates.net)
            event_type, event_location = update_state!(model, event)

            if event_type == 0
                summary_stats[1,1] .+= 1
            end 
            if sum(model.state.I) > summary_stats[1, 2]
               summary_stats[1, 2] .= sum(model.state.I) 
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
end

function simulate(params, i0, aggregation_time = 7, include_log = false, sim = 0)
    #setup 
    patch_names, population_per_patch, mixing_matrix, β, γ = params
    model = SIR(
        patch_names,
        population_per_patch,
        mixing_matrix, 
        β,
        γ)
    model |> x -> rand_infection!(x, i0)
    rates::EventRates = EventRates(model)
    t::Float64 = 0.0
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
        
    #run simulation loop
    sim_loop(model, rates, t, wk, aggregation_time, include_log, tot_data, patch_inf, patch_sus)

    include_log && CSV.write(datadir("$sim log.csv") , sim_log, append=false, header=[:time, :event_type, :event_location])
    return (tot_data,
     patch_sus,
     patch_inf)
end


