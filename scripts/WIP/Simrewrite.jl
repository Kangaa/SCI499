# methods for (gillespie) CTMC simulation of Compartmental disease models
using DrWatson

@quickactivate :SCI499


using CSV
using Tables
using DataFrames
##! SHould npatch be a field of the model?
abstract type CompartmentalModel end

struct SIR <: CompartmentalModel
    S::Vector{Int}
    I::Vector{Int}
end

struct SEIR <: CompartmentalModel
    S::Vector{Int}
    E::Vector{Int}
    I::Vector{Int}
end

struct MixingMatrix
    matrix::Array{Real, 2}
    names::Vector{String}
    popns::Array{Int, 1}
end

abstract type CompartmentalModelparameters end

struct SIRparameters <: CompartmentalModelparameters
    β::Float64
    γ::Float64
    mixmat::MixingMatrix
end

struct SEIRparameters <: CompartmentalModelparameters
    β::Float64
    σ::Float64
    γ::Float64
    mixmat::MixingMatrix
end



function initialise_model(params::SIRparameters, i0::Int64)::SIR
    #first block selects a random patch from those in the model
    n_patch = params.mixmat.names |> length
    rand_choice = rand()
    ix = (rand_choice*n_patch) |> trunc |> Integer
    @assert ix < n_patch "Selected patch $ix of $n_patch, from $rand_choice"
    S0 = copy(params.mixmat.popns)
    I0 = zeros(n_patch) 
    ## and sets the number of exposures to I0
    S0[ix] -= i0
    I0[ix] += i0
    ## then creates a 'Model' object with the provided parameters
    m = SIR(S0, I0)
    # and returns the model
    return m
end

function initialise_model(params::SEIRparameters, i0::Int64)::SEIR
    #first block selects a random patch from those in the model
    n_patch = params.mixmat.names |> length
    rand_choice = rand()
    ix = (rand_choice*n_patch) |> trunc |> Integer
    @assert ix < n_patch "Selected patch $ix of $n_patch, from $rand_choice"
    S0 = copy(params.mixmat.popns)
    E0 = zeros(n_patch) 
    I0 = zeros(n_patch)
    ## and sets the number of exposures to I0
    S0[ix] -= i0
    I0[ix] += i0
    ## then creates a 'Model' object with the provided parameters
    m = SEIR(S0, E0, I0)
    # and returns the model
    return m
end

function event_rates(m::SIR, params::SIRparameters)::Array{Float64, 1}
    rates = zeros(2*(m.S|> length))
    s_frac = m.S ./ params.mixmat.popns
    inf_force = params.β*(transpose(params.mixmat.matrix)*m.I)
    n1 = m.S |> length
    inf_rate = s_frac.*inf_force      
    rec_rate = params.γ*m.I
    rates[firstindex(rates):n1] = inf_rate
    rates[(n1+1):lastindex(rates)] = rec_rate
    rates
end



function event_update(m::SIR, event_ix::Int)
    npatch = m.S |> length
    ev_type = (event_ix  / npatch)|> floor
    ev_location = (event_ix % npatch)+1
    if ev_type == 0
        m.S[ev_location] -= 1
        m.I[ev_location] += 1
    elseif ev_type == 1
        m.I[ev_location] -= 1
    else
        error("Event type $ev_type not recognised")
    end
    Event_type = ev_type |> x -> 
        x == 0 ? "I" :
        x == 1 ? "R" : "error"

    return (m, Event_type, ev_location)
end

function exp_sample(net_rate)
    -(rand() |> log)/net_rate
end

function pick(rates::AbstractArray{Float64,1})::Int
    crates = cumsum(rates)
    r = rand()
    ix = findfirst(x -> x > r, crates)
    ix-1
end
function event_update(m::SEIR, params::SEIRparameters, event_ix::Int)
    npatch = m.S |> length
    ev_type = (event_ix  / npatch)|> floor
    ev_location = (event_ix % npatch)+1
    if ev_type == 0
        m.S[ev_location] -= 1
        m.E[ev_location] += 1
    elseif ev_type == 1
        m.E[ev_location] -= 1
        m.I[ev_location] += 1
    elseif ev_type == 2
        m.I[ev_location] -= 1
    else
        error("Event type $ev_type not recognised")
    end
    Event_type = ev_type |> x -> 
        x == 0 ? "E" :
        x == 1 ? "I" :
        x == 2 ? "R" : "error"

    return (m, Event_type, ev_location)
end



function simulate(params::SIRparameters, i0::Int64)

    m = initialise_model(params, i0)
    N = sum(m.S) + sum(m.I)
    t = 0.0
    nameparams = (N, i0, params.β, params.γ)
    w = CSV.open(datadir(savename(nameparams, "csv")), "w") 
    ## write header
    CSV.write(w,[], header=[:time, :event_type, :event_location])
    while true
        rates = event_rates(m, params)
        net_rate = sum(rates)
        net_rate == 0.0 && break
        dt = exp_sample(net_rate)
        t += dt
        event_ix = pick(rates./net_rate)
        m, event_type, event_location = event_update(m, event_ix)
        write = Tables.table([t  event_type  event_location])
        CSV.write(w, write, append=true)
    end
    close(w)
end


params = SIRparameters(0.5, 0.1, MixingMatrix(fill(1/40, 40, 40), fill("test", 40), fill(100000, 40)))
simulate(params, 3)
