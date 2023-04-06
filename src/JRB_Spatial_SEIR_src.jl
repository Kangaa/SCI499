"""
This is the set of types and methods necessary to run a spatial SEIR Model equivalent to that in Moss et al. 2018
"""
module SpatialSEIR
    using DataFrames
    export MixingMatrix
    export Parameters
    export new_model
    export run_sim
    
struct MixingMatrix
    mm::Array{Real, 2}
    names::Vector{String}
    popns::Array{Int, 1}
end

struct Parameters
    β::Float64
    σ::Float64
    γ::Float64
    mixmat::MixingMatrix
end

mutable struct Model
    params::Parameters
    npatch::Int64
    t::Float64
    S::Vector
    E::Vector
    I::Vector
end


function new_model(params::Parameters, e0::Int64)
    #first block selects a random patch from those in the model
    n_patch = params.mixmat.popns |> length
    rand_choice = rand()
    ix = (rand_choice*n_patch) |> trunc |> Integer
    @assert ix < n_patch "Selected patch $ix of $n_patch, from $rand_choice"
    s0 = copy(params.mixmat.popns)
    ##then creates a 'Model' object with the provided parameters
    m = Model(params,n_patch, 0.0, copy(s0), zeros(n_patch),zeros(n_patch))
        ##and sets the number of exposures to e0
    m.S[ix] -= e0 
    m.E[ix] += e0 
    # and returns the model
    m
end


function event_rates(m::Model)::Array{Float64, 1}
        rates = zeros(3*(m.params.mixmat.popns|> length))
        s_frac = m.S ./ m.params.mixmat.popns
        inf_force = m.params.β*(transpose(m.params.mixmat.mm)*m.I)
        n1 = m.npatch
        n2 = 2*m.npatch
        exp_rate = s_frac.*inf_force      
        inf_rate = m.params.σ*m.E
        rec_rate = m.params.γ*m.I
        rates[firstindex(rates):n1] = exp_rate
        rates[(n1+1):n2] = inf_rate
        rates[(n2+1):lastindex(rates)] = rec_rate
        rates
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


function event_occurred(m::Model, event_ix::Int)
    ev_type = (event_ix  / m.npatch)|> floor
    ev_location = (event_ix % m.npatch)+1
    if ev_type == 0 ## exposure
        m.S[ev_location] -= 1
        m.E[ev_location] += 1
    elseif ev_type == 1 ## Infection
        m.E[ev_location] -= 1
        m.I[ev_location] += 1
    elseif ev_type == 2 ## Recovery
        m.I[ev_location] -= 1
    end
    Event_type = (x ->  x == 0 ? "Exposure" :
                        x == 1 ? "Infection" :
                        x == 2 ? "Recovery" : "error")(ev_type)
    
    return (m, Event_type, ev_location)
end

function run_sim(params::Parameters, e0::Integer)
    m = new_model(params, e0)
    t = 0.0
    week = 0
    Event_log = DataFrame(
            t = Vector{Float64}(),
            Event = Vector{String}(),
            location = Vector{Int64}()
        )

    tot_log =  DataFrame(
            t = Vector{Float64}(),
            S = Vector{Int}(),
            E = Vector{Int}(),
            I = Vector{Int}(),
        )

    inf_log = DataFrame(fill(Vector{Number}(), length(params.mixmat.names)+1),    ["t"; params.mixmat.names])

    while true
        rates = m |> event_rates
        net_rate = sum(rates)
        if net_rate == 0.0
            break
        end
        dt = exp_sample(net_rate) #Sample exponential time
        t += dt
        ix = pick(rates./net_rate)
        (m, Event, Location) = event_occurred(m, ix)
        push!(Event_log, [t, Event, Location])
    end
    return Event_log
 #  save(datadir(savename(params, "csv")), tot_log)
end

end