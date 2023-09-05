

levels = names(data)
    xi = fill(1/length(levels), length(levels))
    npatch = nrow(data)
    level_vec = Vector{String}(undef, npatch)
    mixmat = zeros(npatch, npatch)
    for i in 1:npatch
        patch_i = data[i,:]
        norm = zeros(length(levels))
        for j in 1:npatch
            patch_j = data[j,:]
            for l in 1:length(levels)
                code_l = levels[l]
                if patch_i[Symbol(code_l)] == patch_j[Symbol(code_l)]
                    level_vec[j] = levels[l]
                    norm[l:end] .+= 1
                    break
                end
            end
        end
        for l in 1:length(levels)
            for j in 1:npatch
                if level_vec[j] in levels[1:l]
                    mixmat[i,j] += xi[l]/norm[l]
                end
            end
        end
    end
    return mixmat
end