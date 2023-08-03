function new_model(params::Params, e0::Int)
    #first block selects a random patch from those in the model
    let 
        n_patch = params.popn |> length
        rand_choice = rand() |> nextfloat
        ix = (rand_choice*n_patch) |> trunc
        @assert ix < n_patch "Selected patch $ix of $n_patch, from $rand_choice"
        s0 = params.popn
        ##then creates a 'Model' object with the provided parameters
        m = Model(
            params = params,
            n_patch = n_patch,
            t = 0.0,
            s = s0,
            e = zeros(n_patch),
            i = zeros(n_patch)
            )
            ##and sets the number of exposures to e0
        m.s[ix] -= e0 
        m.e[ix] += e0 
        # and returns the model
        m
    end
end