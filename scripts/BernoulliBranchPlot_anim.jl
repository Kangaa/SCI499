using DataFrames
using Distributions
using Plots
#parameters
    ## of generateion the sim will run for
Gen_max = 20

## Prob of reproducing
p = .6

# TRAJECTORY SIM
function Level_seg(gen, y, branch_num)
        ##gen keeps track of the x[1] position of line segments (i.e. time)
        ## y is the index of which branch (individual) is being simulated (i.e. its yaxis position)


    ## Branch is the output, containing coordinates of individual line segments (branches) in each row (x  y begin and end)
    branch = DataFrame(x = [], y = [], xend = [], yend = [])

    

    gen_remaining = Gen_max - gen - 1

    if gen_remaining == 0
         return branch
    end
    

    if branch_num > 0
        branch = vcat(
            branch,
            DataFrame(x= gen, y = y, xend = gen+1, yend = y + 2^gen_remaining),
            Level_seg(gen+1, y + 2^gen_remaining, rand(Binomial(2, p), 1)[1])
            )
    end
    if branch_num > 1
        branch = vcat(
            branch,
            DataFrame(x= gen, y = y, xend = gen+1, yend = y - 2^gen_remaining),
            Level_seg(gen+1, y - 2^gen_remaining, rand(Binomial(2, p), 1)[1])
            )
    end 
    branch
end

db = Level_seg(0,0, 2)


for j in 1:Gen_max
    bp = subset(db, :x => ByRow(x -> x <= j))

    Xs = fill([0,0], nrow(bp))
    for i in 1:nrow(bp)
        Xs[i] = Vector(bp[i, [1,3]])
    end

    Ys = fill([0,0], nrow(bp))

    for i in 1:nrow(bp)
        Ys[i] = Vector(bp[i, [2,4]])
    end
 display(plot(Xs, Ys, xlim = [0, Gen_max], legend = false))
end

