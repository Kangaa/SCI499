function f(m::Model, x::Array{Float64})::Array{Float64}
    let 
        exp_len = 3*(m.params.popn|>length)
        n1 = m.numpatch
        n2 = 2*m.numpatch
        rates = zeros(exp_len)
        rates_length = rates|>length
        x_length = x|>length
        @assert x|>length == exp_len "Expected length $rates_length, and $x_length"

        s = x[firstindex(x):n1]
        e = x[(n+1):n2]
        i = x[(n2+1):lastindex(x)]
        s_frac = m.s ./ m.params.popn
        inf_force = m.params.β*(transpose(m.params.mixmat)*m.i)
        n1 = m.numpatch
        n2 = 2*m.numpatch
        exp_rate = s_frac.*inf_force      
        inf_rate = m.params.σ*m.e
        rec_rate = m.params.γ*m.i
        ds = -exp_rate
        de = exp_rate - inf_rate
        di = inf_rate - rec_rate
        rates[(n1+1):n2] = ds
        rates[firstindex(rates):n1] = de
        rates[(n2+1):lastindex(rates)] = di
        rates
    end
end