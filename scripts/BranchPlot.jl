# Adapted from Korosteleva (2022)

using DataFrames
using Distributions
using Plots
#parameters
Gen_max = 5
p = 0.4
n = 3

# TRAJECTORY SIM
function Level_seg(gen, y, branch_num)

    branch = DataFrame(x = [], y = [], xend = [], yend = [])
    gen_remaining = Gen_max - gen - 1

    if gen_remaining < 0
         return branch
    end
    
    if branch_num > 0 
        branch = vcat(
            branch,
            DataFrame(x= gen, y = y, xend = gen+1, yend = y),
            Level_seg(gen+1, y, rand(Binomial(n, p), 1)[1]))
    end

    if branch_num > 1
        branch = vcat(
            branch,
            DataFrame(x= gen, y = y, xend = gen+1, yend = y + n^gen_remaining),
            Level_seg(gen+1, y + n^gen_remaining, rand(Binomial(n, p), 1)[1])
            )
    end
    if branch_num > 1
        branch = vcat(
            branch,
            DataFrame(x= gen, y = y, xend = gen+1, yend = y - n^gen_remaining),
            Level_seg(gen+1, y - n^gen_remaining, rand(Binomial(n, p), 1)[1])
            )
    end 
    branch
end

db = Level_seg(0,0, 1)

Xs = fill([0,0], nrow(db))
for i in 1:nrow(db)
    Xs[i] = Vector(db[i, [1,3]])
end

Ys = fill([0,0], nrow(db))
for i in 1:nrow(db)
    Ys[i] = Vector(db[i, [2,4]])
end

plot(Xs, Ys)