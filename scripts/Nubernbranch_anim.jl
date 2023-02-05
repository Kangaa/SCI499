using Pkg
Pkg.activate("Juliaenv")

using Distributions
using Plots
using Animations



gen_max = 50

function Brancho(p = 0.8, gen_max = gen_max)

    branches = [[[0], [0,0], [0.0,0.0]]]

    generation = 0
    gen_pop = 1
    old_gen = branches

    for k in 1: gen_max
            gen_pop = length(old_gen)
            new_gen = [[[k], [0,0], [0.0,0.0]]]
            for i in 1:gen_pop
                reproduce = rand(Binomial(1, p))
                if reproduce == 1
                    childeren = [
                        [ [k],
                            [old_gen[i][2][2],
                            old_gen[i][2][2] + 1],
                            [old_gen[i][3][2],
                            old_gen[i][3][2] + (rand(Uniform(0,1),1))[1]]],
                        [ [k],
                            [old_gen[i][2][2],
                            old_gen[i][2][2] + 1],
                            [old_gen[i][3][2],
                            old_gen[i][3][2] - (rand(Uniform(0,1),1))[1]]]]
                else 
                    childeren = []
                end

                if i == 1
                    new_gen = childeren
                end
                if i > 1 
                    new_gen = [new_gen; childeren]
                end
            end
                branches = [branches; new_gen]
                generation = generation + 1
                old_gen = new_gen

    
    end
    gens = fill([0], length(branches))
    Xs = fill([0,0], length(branches))
    Ys =  fill([0.0,0.0], length(branches))

    for i in 1:length(branches)
        gens[i] = branches[i][1] 
        Xs[i] = branches[i][2] 
        Ys[i] = branches[i][3] 
    end

    [gens, Xs, Ys]

end


db = Brancho(0.5, gen_max)

anim = @animate for i in 1:gen_max
    index = minimum.(db[2]) .<  i
    plot(db[2][index], db[3][index], legend = false)
end

gif(anim)