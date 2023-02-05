using DataFrames
using Distributions
using Plots
#parameters
    ## of generateion the sim will run for
Gen_max = 10

## Prob of reproducing
p = .7

## number of offspring on reproduction
n = 2

# TRAJECTORY SIM
function Level_seg(gen, y, branch_num)
        ##gen keeps track of the x[1] position of line segments (i.e. time)
        ## y is the index of which branch (individual) is being simulated (i.e. its yaxis position)


    ## Branch is the output, containing coordinates of individual line segments (branches)
    branch = DataFrame(x = [], y = [], xend = [], yend = [])

    

    gen_remaining = Gen_max - gen - 1

    if gen_remaining == 0
         return branch
    end
    

    if branch_num > 0
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

db = Level_seg(0,0, 2)

Xs = fill([0,0], nrow(db))
for i in 1:nrow(db)
    Xs[i] = Vector(db[i, [1,3]])
end

Ys = fill([0,0], nrow(db))
for i in 1:nrow(db)
    Ys[i] = Vector(db[i, [2,4]])
end

plot(Xs, Ys)