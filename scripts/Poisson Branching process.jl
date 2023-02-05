using Distributions
using StatsPlots

n₀ = 1
t_max = 50
reps = 30
n = repeat(vcat(1, fill(0, t_max-1)), 1, reps)
λ = 1
n_runs = 10

runs = repeat(zeros(t_max),1,n_runs)
for i in 1:n_runs
    for r in 1:reps
        for t in 2:t_max
            n[t, r] = sum(rand(Poisson(λ), n[t-1, r]))
        end
    end
    runs[:,i] = mean(n, dims = 2)
end

plot(runs)