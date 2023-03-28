using SparseArrays
using ExponentialUtilities

N = 100
# Define the transition rates
β = 2
βi = β/N # infection rate per susceptible individual
γ = 0.1 # recovery rate per infected individual

# Define the state vector
states = [(S,I,R) for S in 0:N, I in 0:N, R in 0:N if (S+I+R)==N]
n_states = length(states)

# Initialize the Q matrix
Q = zeros(n_states, n_states)

# Fill in the Q matrix
for (i,state_i) in enumerate(states)
    S_i, I_i, R_i = state_i
    for (j,state_j) in enumerate(states)
        S_j, I_j, R_j = state_j
        if (S_i - S_j) == 1
            # Infection event
            Q[i,j] = βi*S_i*I_i/N
        elseif (I_i - I_j) == 1
            # Recovery event
            Q[i,j] = γ*I_i
        end
    end
    # Diagonal entries
    Q[i,i] = -sum(Q[i,:])
end


t = 10
x0 = (99, 1, 0)
x0i = findfirst(x -> x == x0, states)

p0 = zeros(lastindex(states))
p0[x0i] = 1

pt = expv(t, Q, p0)
sum(pt)
xt = (0, 0, 100)
xti = findfirst(x -> x == xt, states)
pt[xti]
using StatsPlots

SImat = zeros(N+1, N+1)
pS = zeros(N+1)
pI = zeros(N+1)
pR= zeros(N+1)

for i in eachindex(pt)
    S, I, R = states[i]
    pS[S+1] += pt[i]
    pI[I+1] += pt[i] 
    pR[R+1] += pt[i] 
    SImat[S+1, I+1] = pt[i]
end
plot([pS, pI, pR])

