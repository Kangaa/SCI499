function rates(m::Model)::Array{Float64, 1}
    let 
        rates = zeros(3*(m.params.popn|> length))
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