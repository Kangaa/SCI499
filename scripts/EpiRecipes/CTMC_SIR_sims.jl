using DataFrames
using Distributions
using StatsPlots


function sir(β, γ, N, S0, I0, R0, tf)
    t = 0
    S = S0
    I = I0
    R = R0

    ta = []
    Sa = []
    Ia = []
    Ra = []

    while t < tf

        push!(ta, t)
        push!(Sa, S)
        push!(Ia, I)
        push!(Ra, R)
        
        pf1 = β*S*I
        pf2 = γ*I
        pf = pf1 + pf2
        dt = rand(Exponential(1/pf))
        t = t+dt
        if t > tf
            break
        end
        ru = rand()
        if ru < (pf1/pf)
            S = S-1
            I = I +1
        else    
            I = I-1
            R = R+1
        end

    end

    results = DataFrame()
    results[!, :time] = ta
    results[!, :S] = Sa
    results[!, :I] = Ia
    results[!, :R] = Ra
    return(results)

end

dat = sir(.1, .5, 100, 99, 1, 0, 200)
     
@df dat plot(:time, [:S, :I, :R])
