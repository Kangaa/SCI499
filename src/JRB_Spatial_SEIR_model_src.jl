"""
This is the set of types and methods necessary to run a spatial SEIR Model equivalent to that in Moss et al. 2018
"""

struct OdMat
    m::Array{Float64, 2}
    names::Vector{String}
    popns::Array{Int, 1}
end



function read_popns(path::AbstractString, names::Vector{String})::Vector{Int}
    reader = CSV.read(path,DataFrame)
    popn_tbl = Dict{String, Int}()
    for row in eachrow(reader)
        name = row[2]
        if haskey(popn_tbl, "$name")
            error("Multiple populations for $name")
        else
            value = row[3]
            popn_tbl["$name"] = value
        end
    end
    popns = Vector{Int}()
    for name in names
        if haskey(popn_tbl, "$name")
            push!(popns, popn_tbl["$name"])
        else
            error("No population defined for $name")
        end
    end
    return popns
end

function read_OdMat(mm_path::AbstractString, popn_path::AbstractString)::OdMat
    reader = CSV.read(ODdir, DataFrame)
    names = Vector{String}()
    values = Vector{Float64}()
    for row in eachrow(reader)
        push!(names, string(row[1]))
        for val in row[2:end]
            push!(values, val)
        end
    end
    
    n = length(names)
    mat = reshape(values, (n, n))
    popns = read_popns(popn_path, names)
    
     OdMat(mat, names, popns)
end

abstract type FracSelf end
struct Uniform{Float64} <: FracSelf 
    a::Float64
end
struct Variable{Array} <: FracSelf end

struct MixMat
    mm::Array{Real, 2}
    names::Vector{String}
    popns::Array{Int, 1}
    frac_self::FracSelf
    frac_cbd::Real
end


using LinearAlgebra
function scale_diag_var!(mat::Matrix{T}, frac_self::Vector{T}) where T<:AbstractFloat
    @inbounds for i in 1:size(mat, 1)
        mat[i, :] = mat[i, :] .* ((1.0 - frac_self[i]) ./ sum(mat[i, :]))
    end
    diag(mat) .= frac_self
end

function scale_diag_unif!(mat::Matrix{T}, frac_self) where T<:AbstractFloat
    frac_mix = 1.0 - frac_self
    @inbounds for i in 1:size(mat, 1)
        mat[i, :] = mat[i, :] .* (frac_mix ./ sum(mat[i, :]))
    end
    diag(mat) .= frac_self
end

function scale_diag!(mat::Matrix{T}, frac_self) where T<:AbstractFloat
    if isa(frac_self, Uniform)
        scale_diag_unif!(mat, frac_self.a)
    elseif isa(frac_self, Variable)
        scale_diag_var!(mat, frac_self.x)
    end
end

function test_scale_diag(mat::Matrix{T}, frac_self::T) where T<:AbstractFloat
    m1 = copy(mat)
    scale_diag_unif!(m1, frac_self)

    m2 = copy(mat)
    fs2 = fill(frac_self, size(mat, 1))
    scale_diag_var!(m2, fs2)

    m3 = copy(mat)
    fs3 = fill(frac_self, size(mat, 1))
    fs3[2] *= 0.9
    scale_diag_var!(m3, fs3)

    @assert all(abs.(1.0 .- sum(m1, dims=2)) .< 1e-8) "m1 rowsum"
    @assert all(abs.(1.0 .- sum(m2, dims=2)) .< 1e-8) "m2 rowsum"
    @assert all(abs.(1.0 .- sum(m3, dims=2)) .< 1e-8) "m3 rowsum"
    @assert m1 == m2 "mixing matrices differ"
    @assert m1 != m3 "mixing matrices do not differ"
end

function redistribute!(mat::Matrix{Float64}, popns::Vector{Int},
    frac_cbd::Float64, cbd_ix::Int)
    # Determine how much all other regions mix with the CBD.
    mix_w_cbd = copy(mat[:, cbd_ix])
    mix_w_cbd[cbd_ix] = 0.0

    # Weight the mixing by each region's resident population.
    mix_w_cbd .= mix_w_cbd .* Float64.(popns)

    # Normalise these values to conserve the force of infection.
    mix_w_cbd /= sum(mix_w_cbd)

        for rix in 1:size(mat, 1)
            if rix == cbd_ix
                continue
            end

        row = mat[rix, :]
        cbd_mix = row[cbd_ix]
        add_mix = [cbd_mix * (1.0 - frac_cbd) * v for v in mix_w_cbd]
        row += add_mix
        row[cbd_ix] = cbd_mix * frac_cbd
        mat[rix, :] = row
    end
end


function new_MixMat(od::OdMat, frac_self::FracSelf, frac_cbd::Float64)::MixMat
    mat = copy(od.m)
    names = copy(od.names)
    popns = copy(od.popns)
    scale_diag!(mat, frac_self)
    cbd_ix = findfirst(isequal("20604"), names)
    redistribute!(mat, popns, frac_cbd, cbd_ix)
    MixMat(mat, names, popns, frac_self, frac_cbd)
end


struct Params
    β::Float64
    σ::Float64
    γ::Float64
    mixmat::MixMat
end


mutable struct Model
    params::Params
    numpatch::Int
    t::Float64
    s::Vector
    e::Vector
    i::Vector
end


function new_model(params::Params, e0::Int)
    #first block selects a random patch from those in the model
    n_patch = params.mixmat.popns |> length
    rand_choice = rand() |> nextfloat
    ix = (rand_choice*n_patch) |> trunc |> Integer
    @assert ix < n_patch "Selected patch $ix of $n_patch, from $rand_choice"
    s0 = copy(params.mixmat.popns)
    ##then creates a 'Model' object with the provided parameters
    m = Model(params, n_patch, 0.0, copy(s0), zeros(n_patch),zeros(n_patch))
        ##and sets the number of exposures to e0
    m.s[ix] -= e0 
    m.e[ix] += e0 
    # and returns the model
    m
end


function event_rates(m::Model)::Array{Float64, 1}
        rates = zeros(3*(m.params.mixmat.popns|> length))
        s_frac = m.s ./ m.params.mixmat.popns
        inf_force = m.params.β*(transpose(m.params.mixmat.mm)*m.i)
        n1 = m.numpatch
        n2 = 2*m.numpatch
        exp_rate = s_frac.*inf_force      
        inf_rate = m.params.σ*m.e
        rec_rate = m.params.γ*m.i
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


function event_occurred(m::Model, event_ix::Int)::Model
    ev_type = (event_ix  / m.numpatch)|> floor
    ev_location = (event_ix % m.numpatch)+1
    if ev_type == 0 ## exposure
        m.s[ev_location] -= 1
        m.e[ev_location] += 1
    elseif ev_type == 1 ## Infection
        m.e[ev_location] -= 1
        m.i[ev_location] += 1
    elseif ev_type == 2 ## Recovery
        m.i[ev_location] -= 1
    end
    m
end

function run_sim(params::Params, e0)
    m = new_model(params, e0)
    t = 0.0
    week = 0
    weekly_log = DataFrame(
            week = Vector{Float64}(),
            S = Vector{Int}(),
            E = Vector{Int}(),
            I = Vector{Int}(),
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
        m = event_occurred(m, ix)
        push!(tot_log, [t, sum(m.s), sum(m.e), sum(m.i)])
        push!(inf_log, [t; m.i])
    end
    return tot_log, inf_log
 #  save(datadir(savename(params, "csv")), tot_log)
end