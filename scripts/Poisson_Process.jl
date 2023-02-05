using Distributions


t = 50.0
t₀ = 0
n₀ = 0
λ = 0.75

n = rand(Poisson(λ*t))

A = rand(Exponential(λ), n)
iA = cumsum(A)

plot(iA)
