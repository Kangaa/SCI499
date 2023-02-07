using Distributions
using Plots
using Animations


function poisson_BP(λ, gen_max = gen_max)

    branches = [[[0], [0,0], [0.0,0.0]]]

    generation = 0
    gen_pop = 1
    old_gen = branches

    for k in 1: gen_max
            gen_pop = length(old_gen)
            new_gen = [[[k], [0,0], [0.0,0.0]]]
            for i in 1:gen_pop
                nchild = rand(Poisson(λ), 1)[1]
                children = fill([[], [], []], nchild)
                if nchild == 0
                    children = []
                elseif nchild > 0
                    for j in 1:nchild
                        children[j] = [
                            [k],
                                [old_gen[i][2][2],
                                old_gen[i][2][2] + 1],
                                [old_gen[i][3][2],
                                old_gen[i][3][2] + (rand(Uniform(-1,1),1))[1]]]
                    end
                end

                if i <= 1
                    new_gen = children
                end
                if i > 1 
                    new_gen = [new_gen; children]
                end
            end
                branches = [branches; new_gen]
                generation = generation + 1
                old_gen = new_gen

    
    end
    gens = fill([0.0], length(branches))
    Xs = fill([0.0,0.0], length(branches))
    Ys =  fill([0.0,0.0], length(branches))

    for i in 1:length(branches)
        gens[i] = branches[i][1] 
        Xs[i] = branches[i][2] 
        Ys[i] = branches[i][3] 
    end

    [gens, Xs, Ys]

end

gen_max = 150
db = poisson_BP(1.1, gen_max)

plot(db[2], db[3], legend = false)
