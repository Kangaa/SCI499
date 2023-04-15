using DrWatson
@quickactivate :SCI499
using CSV
using Tables
using StaticArrays
using LinearAlgebra

# Create an abstract type for the model compartments
abstract type Compartment{P} end
# Define compartments for SIR, SEIR, and SIS models
struct Susceptible{P} <: Compartment{P}
    population::MVector{P, Int}
end

struct Infected{P} <: Compartment{P}
    population::MVector{P, Int}
end

struct Recovered{P} <: Compartment{P}
    population::MVector{P, Int}
end

struct Exposed{P} <: Compartment{P}
    population::MVector{P, Int}
end

# common compartmental types


# Define transition parameters

abstract type TransitionParameter{P} end

struct β{P} <: TransitionParameter{P}
    value::SVector{P, Float64}
end

struct γ{P} <: TransitionParameter{P}
    value::SVector{P, Float64}
end

struct σ{P} <: TransitionParameter{P}
    value::SVector{P, Float64}
end

# Create a CompartmentalModel type to hold the state and parameters of the model
struct CompartmentalModel{C, P, S, T} ## C is the number of compartments, P is the number of patches
    state::NamedTuple{S,  NTuple{C, Compartment{P}}}
    transition_parameters::NamedTuple{T, NTuple{C, TransitionParameter{P}}}
    total_population::Int
    num_patches::Int
    patch_names::SVector{P, String}
    population_per_patch::SVector{P, Int}
    mixing_matrix::SMatrix{P,P}
end


function SIR(patch_names::SVector{P, String}, population_per_patch::SVector{P, Int}, mixing_matrix::SMatrix{P, P}, beta::Float64, gamma::Float64) where {P}
    state = (
        S = Susceptible{P}(MVector{P, Int}(population_per_patch)),
        I = Infected{P}(MVector{P, Int}(zeros(Int, P))))
    transition_parameters = (
        β = β{P}(SVector{P, Float64}(fill(beta, P))),
        γ = γ{P}(SVector{P, Float64}(fill(gamma, P))))
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
    model.state.I.population[patch] += num_infected
    # Set the number of susceptible individuals in the patch
    model.state.S.population[patch] -= num_infected
    return model
end

## Add a type to hold the event rates
struct EventRates{L}
    rates::MVector{L, Float64} ## is this the best way to store rates (or as perEvent type)

    function EventRates(Model::CompartmentalModel{C, P, S, T}) where {C, P, S, T}
        L = C*P
        rates = MVector{L, Float64}(zeros(L))
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
    s_frac::Float64 = model.state.S.population./model.population_per_patch
    ##force of infection
    f = model.transition_parameters.β.value[1]*(transpose(model.mixing_matrix)*model.state.I.population)
    ## Infection rate
    rates.rates[1:P] = s_frac.*f
    ## Recovery rate
    rates.rates[P+1:2P] = model.state.I.population.*model.transition_parameters.γ.value[1]
    return rates
end


## update the model state based on the event rates
function update_state!(model::CompartmentalModel{C, P, S, T}, event_ix) where {C, P, S, T}
    ###pick an event to occur using gillespie algorithm
    event_type = (event_ix/model.num_patches)|> floor
    event_location = (event_ix%model.num_patches) +1
    if event_type == 0
        model.state.S.population[event_location] -= 1
        model.state.I.population[event_location] += 1
    elseif event_type == 1
        model.state.I.population[event_location] -= 1
    end
    return model   
end
function exp_sample(net_rate)::Float64
    -(rand() |> log)/net_rate
end
## Define a function to simulate a new model
function simulate(params, i0)
    model = SIR(params.patch_names, params.population_per_patch, params.mixing_matrix, params.beta, params.gamma)
    model|> x -> rand_infection!(x, i0)
    rates = EventRates(model)
    t = 0.0
    update_rates!(model, rates)
    w = CSV.open(datadir(savename(params, "csv")), "w") 
    ## write header
    CSV.write(w,[], header=[:time, :event_type, :event_location])
    while true
        update_rates!(model, rates)
        net_rate = sum(rates.rates)
        if net_rate == 0.0
            break
        end
        dt = exp_sample(net_rate)
        t += dt
        event_ix = gillespie_pick(rates.rates./net_rate)
        update_state!(model, event_ix)
        write = Tables.table([t; model.state.I.population]')
        CSV.write(w, write, append=true)
    end
    return 
end



Gmelb_SA3_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA3_2021_AUST_SHP_GDA2020/SA3_2021_AUST_GDA2020.shp")) |> x ->
    Tables.subset(x, x.GCC_NAME21 .== "Greater Melbourne")|>
    DataFrame|> x-> 
    rename(x,
        "SA3_NAME21" => "SA3_NAME", 
        "SA3_CODE21" => "SA3_CODE" )




codes = Gmelb_SA3_SHP |>
    names |>
    (x-> contains.(x, "CODE")) |>
    (x-> select(Gmelb_SA3_SHP, x))|>
    (x -> select(x, 1:2))

##convert missing entries to 0

params =(
    patch_names = SVector{40}(Gmelb_SA3_SHP.SA3_NAME),
    population_per_patch = SVector{40}(fill(1000, 40)),
    mixing_matrix = SMatrix{40,40}(SpatialMixingMatrix(codes, 0.5)),
    beta = 0.5,
    gamma = 0.1)


    using Cthulhu
@descend simulate(params, 3)
