using Match
using DataFrames

function log_infs(m::Model, x::Vector{Float64}, mixmat::MixMat, week::Int)
    global Cumulative_infections = DataFrame(
        Fraction_self = Vector{Real}(),
        Fraction_CBD = Vector{Real}(),
        sim = Vector{Int}(),
        week = Vector{Real}(),
        name = Vector{String}(),
        Cumulative_infections = Vector{Int}(),
        proportion_n = Vector{Int}(),
        population_n = Vector{Int}()
    )
    let n1 = m.numpatch
        n2 = 2*m.numpatch
        s = x[firstindex(x):(n1)]
        e = x[(n1+1):n2]
        cum_infs = m.params.popn - (s + e) 
        @assert length(cum_infs) == length(mixmat |> names) "Have $(length(cum_infs)) patches and $(length(mixmat.names)) patch names"
        @assert length(cum_infs) == length(m.params.popn) "Have $(length(cum_infs)) patches and $(length(m.params.popn)) patch names"
        frac_self =  @match mixmat.frac_self begin
            u::Uniform(val) => FracVal
            v::Variable => 0.0
        end
        for ix in cum_infs
            cum_infs[ix] = cinf

                propn = cinf
                name = mixmat.names[ix]
                popn = m.params.popn[ix]
                sim = 0
            push!(Cumulative_infections, [frac_self, mixmat.frac_cbd, sim, week, name, cinf, propn, popn])
        end
    end
end