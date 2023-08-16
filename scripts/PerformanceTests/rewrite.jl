## SIR metapopulation model
using CSV 
using DataFrames
using LinearAlgebra
### Model types
abstract type CompartmentalModel end

struct SIR_State
    S::Vector{Int64}
    I::Vector{Int64}
end

struct SIR_Params
    β::Float64
    γ::Float64
end

struct Compartment_Patches
    Name::Vector{String}
    popn::Vector{Int64}
end

mutable struct SIRModel <: CompartmentalModel
    state::SIR_State
    params::SIR_Params 
    patches::Compartment_Patches
    mixing_matrix::Matrix{Float64}
    foi_vector::Vector{Float64}
end



## model constructor

init_SIR(codepops, mixmat, β, γ) = SIRModel(
    SIR_State(
        codepops[:, 2],
        zeros(Int64, length(codepops[:, 2]))
    ),
    SIR_Params(
        β,
        γ
    ),
    Compartment_Patches(
        codepops[:, 1],
        codepops[:, 2]
    ),
    mixmat,
    zeros(Float64, length(codepops[:, 2]))
)
BLAS.set_num_threads(1)


## Function to seed the model with a number of infected individuals at a random patch
function rand_infection!(model::SIRModel, num_infected::Int64)::CompartmentalModel
    # Randomly select a patch to seed
    patch::Int64 = rand(1:length(model.patches.Name))
    # Set the number of infected individuals in the patch
    model.state.I[patch]::Int64 += num_infected
    # Set the number of susceptible individuals in the patch
    model.state.S[patch]::Int64 -= num_infected
    return model
end

# Seed infection

function update_rates!(model::SIRModel, rates::Vector{Float64})::Vector{Float64}
    ## fraction of susceptible in each patch
    s_frac::Vector{Float64} =  model.state.S./model.patches.popn
    ##force of infection
    model.foi_vector::Vector{Float64} .= model.params.β .* (model.mixing_matrix' * model.state.I)
    ## Infection rate
    rates[1:numpatch]::Vector{Float64} = s_frac.*model.foi_vector
    ## Recovery rate
    rates[numpatch+1:(2*numpatch)]::Vector{Float64} = model.state.I.*model.params.γ

    return rates
end


using Random
rng = MersenneTwister(1234)

function exp_sample(net_rate)::Float64
    -(rand(rng) |> log)/net_rate
end

function gillespie_pick(rates)::Int
    crates = cumsum(rates)
    r = rand(rng, Float64)
    ix = findfirst(x -> (x > r), crates)
    ix - 1
end

## update the model state based on the event rates
function update_state!(model::SIRModel, return_event = true)
    ###pick an event to occur using gillespie algorithm
    event_ix::Int64 = gillespie_pick(rates./sum(rates))

    event_type::Int64 = (event_ix/numpatch)|> floor

    event_location = (event_ix%numpatch) +1
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



function sim_loop(model::SIRModel, rates::Vector{Float64}, t::Float64, wk::Int64, aggregation_time::Int64, tot_data, patch_inf, patch_sus, summary_stats)
    while sum(rates) != 0.0
        update_rates!(model, rates)
        if sum(rates) != 0.0
            t += exp_sample(sum(rates))
            update_state!(model, true)


            if t/aggregation_time > wk
                wk += 1
                push!(tot_data, (time=wk, TotalSusceptible=sum(model.state.S), TotalInfected=sum(model.state.I)), promote = true)
                push!(patch_sus, model.state.S)
                push!(patch_inf, model.state.I)
            end
        end

    end
end

## define params

codepops = CSV.read("data/GmelbSA2Pop21.csv", DataFrame, header = 1,  types=[String, String, Int64]) |> 
    x -> select(x, 1, 3)

codepops[:, 2] = codepops[:, 2] .+ 1

mixmat = CSV.read("data/SA2_Mixmat_0.5_0.5.csv", DataFrame)|> 
    x -> Matrix(x)
    
β = 1.5
γ = 1.0
model = init_SIR(codepops, mixmat, β, γ)
global numpatch = length(model.patches.Name)

rand_infection!(model, 10)
rates = Vector{Float64}(zeros(nfields(model.state)*length(model.patches.Name)))

t = 0.0
wk::Int64 = 0
update_rates!(model, rates)
sum(rates)

tot_data = DataFrame(
    time = Int64,
    TotalSusceptible = Int64,
    TotalInfected = Int64)

patch_inf =  DataFrame(Dict(name => Vector{Int}() for name in model.patches.Name))
patch_sus =  DataFrame(Dict(name => Vector{Int}() for name in model.patches.Name))
    
summary_stats = DataFrame(
    total_infections = 0,
    peak_infections = 0,
    duration = 0)

#run simulation loop

using Profile
using PProf 

@pprof sim_loop(model, rates, t, wk, 7, tot_data, patch_inf, patch_sus, summary_stats)



