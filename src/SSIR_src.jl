# New source functions 

using DrWatson
@quickactivate :SCI499
using CSV
using Tables
using StaticArrays
using LinearAlgebra
using ProfileView



# Create a CompartmentalModel type to hold the state and parameters of the model
struct CompartmentalModel{C, P, S, T} ## C is the number of compartments, P is the number of patches
    state::NamedTuple{S,  NTuple{C, MVector{P, Int64}}}
    transition_parameters::NamedTuple{T, NTuple{C, MVector{P, Float64}}}
    total_population::Int
    num_patches::Int
    patch_names::SVector{P, String}
    population_per_patch::SVector{P, Int}
    mixing_matrix::SMatrix{P,P}
end


function SIR(patch_names::SVector{P, String}, population_per_patch::SVector{P, Int}, mixing_matrix::SMatrix{P, P}, beta::Float64, gamma::Float64) where {P}
    state = (
        S = MVector{P, Int}(population_per_patch),
        I = MVector{P, Int}(zeros(Int, P)))
    transition_parameters = (
        β = SVector{P, Float64}(fill(beta, P)),
        γ = SVector{P, Float64}(fill(gamma, P)))
    total_population = sum(population_per_patch)
    num_patches = length(patch_names)
    patch_names = SVector{P, String}(patch_names)
    population_per_patch = SVector{P, Int}(population_per_patch)
    mixing_matrix = SMatrix{P,P}(mixing_matrix)
    return CompartmentalModel{2, P, (:S, :I), (:β, :γ)}(state, transition_parameters, total_population, num_patches, patch_names, population_per_patch, mixing_matrix)
end

## Function to seed the model with a number of infected individuals at a random patch
function rand_infection!(model::CompartmentalModel{C, P, S, T}, num_infected::Int) where {C, P, S, T}
    # Randomly select a patch to seed
    patch = rand(1:P)
    # Set the number of infected individuals in the patch
    model.state.I[patch] += num_infected
    # Set the number of susceptible individuals in the patch
    model.state.S[patch] -= num_infected
    return model
end

## Add a type to hold the event rates
mutable struct EventRates{L}
    rates::MVector{L, Float64} ## is this the best way to store rates (or as perEvent type)
    net::Float64

    function EventRates(Model::CompartmentalModel{C, P, S, T}) where {C, P, S, T}
        L = C*P
        rates = MVector{L, Float64}(zeros(L))
        net = 0.0
        return new{L}(rates)
    end
end

function gillespie_pick(rates)::Int
    crates = cumsum(rates)
    r = rand()
    ix = findfirst(x -> x > r, crates)
    ix - 1
end

## Update the event rates based on state

function update_rates!(model::CompartmentalModel{C, P, S, T}, rates::EventRates) where {C, P, S, T}
    ## fraction of susceptible in each patch
    @inbounds s_frac::MVector{P, Float64} =  model.state.S./model.population_per_patch
    ##force of infection
    @inbounds f::MVector{P, Float64} = model.transition_parameters.β[1]*(transpose(model.mixing_matrix)*model.state.I)
    ## Infection rate
    @inbounds rates.rates[1:P] = s_frac.*f
    ## Recovery rate
    @inbounds rates.rates[P+1:2P] = model.state.I.*model.transition_parameters.γ[1]
    @inbounds rates.net = sum(rates.rates)

    return rates
end


## update the model state based on the event rates
function update_state!(model::CompartmentalModel{C, P, S, T}, event_ix, return_event = true) where {C, P, S, T}
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

## Define a function to simulate a new model
using DataFrames
function simulate(params, i0, savedir = "data/test.csv")
    model = SIR(params.patch_names, params.population_per_patch, params.mixing_matrix, params.beta, params.gamma)
    model|> x -> rand_infection!(x, i0)
    rates = EventRates(model)
    t = 0.0
    update_rates!(model, rates)
    sim_data = DataFrame(time = Float64[], event_type = Int[], event_location = Int[])
    while true
        update_rates!(model, rates)
        if rates.net == 0.0
            break
        end
        t += exp_sample(rates.net)
        event_ix = gillespie_pick(rates.rates./rates.net)
        (event_type, event_location) = update_state!(model, event_ix)
        sim_data = push!(sim_data, (time=t, event_type=event_type, event_location=event_location))
    end
    CSV.write(datadir("test.csv")  , sim_data, append=true, header=[:time, :event_type, :event_location])
    return 
end