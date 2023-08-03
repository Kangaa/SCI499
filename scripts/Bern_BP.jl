using Distributions
using Plots

function Brancho(p = 0.8, gen_max = 10)

    branches = [[[0,0], [0.0,0.0]]]

    generation = 0
    gen_pop = 1
    old_gen = branches

    for k in 1: gen_max
            gen_pop = length(old_gen)
            new_gen = [[[0,0], [0.0,0.0]]]
            for i in 1:gen_pop
                reproduce = rand(Binomial(1, p))
                if reproduce == 1
                    childeren = [
                        [
                            [old_gen[i][1][2],
                            old_gen[i][1][2] + 1],
                            [old_gen[i][2][2],
                            old_gen[i][2][2] + (rand(Uniform(0,1),1))[1]]],
                        [
                            [old_gen[i][1][2],
                            old_gen[i][1][2] + 1],
                            [old_gen[i][2][2],
                            old_gen[i][2][2] - (rand(Uniform(0,1),1))[1]]]]
                else 
                    childeren = []
                end

                if i == 1
                    new_gen = childeren
                elseif i > 1 
                    new_gen = [new_gen; childeren]
                end
            end
                branches = [branches; new_gen]
                generation = generation + 1
                old_gen = new_gen
        
    Xs = fill([0,0], length(branches))
    Ys =  fill([0.0,0.0], length(branches))

    for i in 1:length(branches)
        Xs[i] = branches[i][1] 
        Ys[i] = branches[i][2] 
    end
    display(plot(Xs, Ys, legend = false, xlims = [0, gen_max]))

    
    end

    branches


end

db = Brancho(0.9)


