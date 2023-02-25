using  Random

function run_sims(od::OdMat, f_self::Vector{FracSelf}, f_cbd::Vector{Float64}, n::Int)
    rng = Random.MersenneTwister(1234567)
    r0 = 1.4
    σ = 2.0
    γ = 0.5

    for frac_self in f_self
        for frac_cbd in f_cbd
            mat_info::MixMat = new(od, frac_self, frac_cbd)
            popns = mat_info.popns
            mixmat = mat_info.m
            for sim in 0:n
                accepted = false
                while ! accepted
                    params = Params(r0, σ, γ, popns, mixmat)
                    accepted = run_stoch(params, mat_info, rng, sim)
            end
end