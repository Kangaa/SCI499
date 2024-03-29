using DrWatson
@quickactivate :SCI499
using CSV
using Tables
using StaticArrays
using LinearAlgebra
using ProfileView


# Create a CompartmentalModel type to hold the state and parameters of the model
struct CompartmentalModel## C is the symbols for the compartments, P is the symbols for the parameters
    state::NamedTuple
    transition_parameters::NamedTuple
    total_population::Int64
    num_patches::Int64
    patch_names::Vector{String}
    population_per_patch::Vector{Int64}
    mixing_matrix::Matrix{Float64}
end

CompartmentalModel((S = 10, I = 10), (β = [0.5], γ = [0.5]), 20, 2, ["Patch1", "Patch2"], [10, 10], [0.5 0.5; 0.5 0.5])

function SIR(patch_names::Vector{String}, population_per_patch::Vector{Int64}, mixing_matrix::Matrix{Float64}, β::Float64, γ::Float64)
    num_patches = length(patch_names)
    state = (
        S = Vector{Int}(population_per_patch),
        I = Vector{Int}(zeros(Int, num_patches)))
    transition_parameters = (
        β = Vector{Float64}(fill(β, num_patches)),
        γ = Vector{Float64}(fill(γ, num_patches)))
    total_population = sum(population_per_patch)
    num_patches = length(patch_names)
    patch_names = Vector{String}(patch_names)
    population_per_patch = Vector{Int}(population_per_patch)
    mixing_matrix = mixing_matrix
    return CompartmentalModel(state, transition_parameters, total_population, num_patches, patch_names, population_per_patch, mixing_matrix)
end
model = SIR(["Patch1", "Patch2"], [1000, 1000], [0.5 0.5; 0.5 0.5], 0.5, 0.5)

## Function to seed the model with a number of infected individuals at a random patch
function rand_infection!(model::CompartmentalModel, num_infected::Int)
    # Randomly select a patch to seed
    patch = rand(1:model.num_patches)
    # Set the number of infected individuals in the patch
    model.state.I[patch] += num_infected
    # Set the number of susceptible individuals in the patch
    model.state.S[patch] -= num_infected
    return model
end
model|> x -> rand_infection!(x, 10)

## Add a type to hold the event rates
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
    ix - 1
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
## Define a function to simulate a new model
using DataFrames

function sim_loop(model::CompartmentalModel, rates::EventRates, t::Float64, wk::Int64, aggregation_time::Int64, include_log::Bool, data...)
    sim_data, sim_log = data
    while rates.net != 0.0
        update_rates!(model, rates)
        if rates.net != 0.0
            t += exp_sample(rates.net)
            event_ix = gillespie_pick(rates.rates./rates.net)
            event_type, event_location = update_state!(model, event_ix)
            if include_log
                push!(sim_log, (time = t, event_type = event_type, event_location = event_location))
            end
            if t/aggregation_time > wk
                wk += 1
                push!(sim_data, (time=wk, TotalSusceptible=sum(model.state.S), TotalInfected=sum(model.state.I)), promote = true)
            end
        end
    end
end

function simulate(params, i0, aggregation_time = 7, include_log = false)
    #setup 
    model = SIR(params.patch_names, params.population_per_patch, params.mixing_matrix, params.β, params.γ)
    model |> x -> rand_infection!(x, i0)
    rates = EventRates(model)
    t = 0.0
    wk::Int64 = 0
    update_rates!(model, rates)

    ## Create log and data frames
    if include_log
          sim_log = DataFrame(time = Float64[], event_type = Int[], event_location = Int[])
    end
    sim_data = DataFrame(time = Int64, TotalSusceptible = Int64, TotalInfected = Int64)

    #run simulation loop
    sim_loop(model, rates, t, wk, aggregation_time, include_log, sim_data, sim_log)

    # Write to CSV
    CSV.write(datadir("test.csv") , sim_data, append=false, header=[:time, :TotalSusceptible, :TotalInfected])
    include_log && CSV.write(datadir("test_log.csv") , sim_log, append=false, header=[:time, :event_type, :event_location])
    return 
end

## test params
params = (
    patch_names = ["Patch 1", "Patch 2", "Patch 3"],
    population_per_patch = [1000000, 1000000, 1000000],
    mixing_matrix = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0],
    β = 0.1,
    γ = 0.1)

@time simulate(params, 10, 7, true)

Gmelb_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp")) |> x ->
    Tables.subset(x, x.GCC_NAME21 .== "Greater Melbourne")|>
    DataFrame|> x-> 
    rename(x,
        "SA2_NAME21" => "SA2_NAME", 
        "SA2_CODE21" => "SA2_CODE" )

Vic_pop = CSV.read(datadir("VicPop21.csv"), DataFrame) 

Gmelb_pop = Vic_pop|> x -> 
    transform(x, "SA2 code" => (y -> string.(y)) => "SA2_CODE") |> x->
    rename(x, "2021" => "popn")

Gmelb =leftjoin(Gmelb_SA2_SHP, Gmelb_pop[:,["SA2_CODE","popn" ]], on = "SA2_CODE")

Gmelb.popn = ((Gmelb.popn) .|> floor) |> x -> convert(Vector{Int64}, x)

##! Convert patches with population < 100 to 100 (to avoid zero population patches)
Gmelb.popn = Gmelb.popn .|> x -> x < 100 ? 100 : x

codes = Gmelb |>
    names |>
    (x-> contains.(x, "CODE")) |>
    (x-> select(Gmelb, x))|>
    (x -> select(x, 1:3))

mixmat = SpatialMixingMatrix(codes, 0.5)

params = (
    patch_names = Gmelb.SA2_NAME,
    population_per_patch = Gmelb.popn,
    mixing_matrix = mixmat,
    β = 0.5,
    γ = 0.5)

simulate(params, 10,7,false)