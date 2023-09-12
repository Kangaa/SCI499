##SIR Metapopulation Model ODE
using ModelingToolkit
using OrdinaryDiffEq
using Plots

include("../src/SCI499.jl")
npatch = 

@parameters t
D = Differential(t)

@parameters t β c[1:npatch] γ;

@variables (S(t))[1:npatch] (I(t))[1:npatch] (R(t))[1:npatch] (λ(t))[1:npatch]

state_eqs = [
[D(S[i]) ~ -λ[i]*S[i],
D(I[i]) ~ λ[i]*S[i]-γ*I[i],
D(R[i]) ~ γ*I[i]]
for i in 1:npatch]
# Turn into a 1D vector
state_eqs = vcat(state_eqs...);

@variables (N(t))[1:npatch] (p(t))[1:npatch,1:npatch]

variable_eqs = [
    [N[i] ~ S[i]+I[i]+R[i] for i in 1:npatch]...,
    [λ[i] ~ sum([β*c[i]*p[i,j]*I[j]/N[j] for j in 1:npatch]) for i in 1:npatch]...,
    [p[i,j] ~ c[j]*N[j]/sum([c[k]*N[k] for k in 1:npatch]) for j in 1:npatch for i in 1:npatch]...
];

@named  sys = ODESystem([state_eqs; variable_eqs])

simpsys = structural_simplify(sys)

u₀ = [[S[i] => 990.0/npatch for i in 1:npatch]...,
      [I[i] => 10.0/npatch for i in 1:npatch]...,
      [R[i] => 0.0 for i in 1:npatch]...];

p = [β => 0.05, γ => 0.25, [c[i] => 10.0 for i in 1:npatch]...];

δt = 0.1
tmax = 40
tspan = (0.0, tmax);

prob = ODEProblem(simpsys, u₀, tspan, p)
sol = solve(prob, Tsit5(), saveat=δt);

all_states = states(simpsys)

indexof(sym,syms) = findfirst(isequal(sym),syms)
S_indexes = [indexof(S[k],all_states) for k in 1:npatch]
I_indexes = [indexof(I[k],all_states) for k in 1:npatch]
R_indexes = [indexof(R[k],all_states) for k in 1:npatch];

Smat = sol[S_indexes,:]
Imat = sol[I_indexes,:]
Rmat = sol[R_indexes,:];

Stotal = sum(Smat,dims=1)'
Itotal = sum(Imat,dims=1)'
Rtotal = sum(Rmat,dims=1)';

times = sol.t
plot(times, Stotal, label="S", xlabel="Time", ylabel="Number");
plot!(times, Itotal, label="I");
plot!(times, Rtotal, label="R")


p2 = [β=>0.05, c[1] => 20, c[2] => 5, γ=>0.25]
prob2 = remake(prob, p=p2)
sol2 = solve(prob2, Tsit5(), saveat=δt);

plot(times, sol2(times, idxs=S_indexes)', labels=["S₁" "S₂"], linecolor=:blue, linestyle=[:solid :dash]);
plot!(times, sol2(times, idxs=I_indexes)', labels=["I₁" "I₂"], linecolor=:red, linestyle=[:solid :dash]);
plot!(times, sol2(times, idxs=R_indexes)', labels=["R₁" "R₂"], linecolor=:green, linestyle=[:solid :dash])
xlabel!("Time")
ylabel!("Number")