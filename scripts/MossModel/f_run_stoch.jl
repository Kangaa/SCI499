function run_stoch(params::Params, mixmat::MixMat, sim::Int)
    let infs = 10;
        rng = rand()
        m = new_model(params, infs)
        t = 0.0
        week = 0

        while true
            rates = m|>event_rates
            net_rate = sum(rates)
            if net_rate == 0.0 
                break
            end
            
            dt = -(rng |> nextfloat |> log)/net_rate #Sample exponential time
            t += dt

            rand_event = net_rate*(rng |> nextfloat)
            ix = pick(rates, net_rate, length(rates))
            event_occurred(ix)
        end

        if !accept_simulation(m)
            return false
        end

        week = (t/7)|>floor
        log_infs(mixmat, sim, week)
        return true
    end  
end