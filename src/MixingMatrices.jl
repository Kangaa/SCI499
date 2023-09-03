#Functions for constructing mixing matrices
using CSV
using DataFrames

struct MixingMatrix
    mm::Array{Real, 2}
    names::Vector{String}
    popns::Array{Int, 1}
end

function HMixingMatrix(codes::Vector{String}, ξ::Vector{Float64})
    Codes = expandCodes(codes)
    npatch = size(Codes, 1)
    HMM = fill(0.0, npatch, npatch)

    levels = names(Codes)
    nlevels = length(levels)
    norm_vec = fill(0, nlevels)
    level_vec = fill("", npatch)

    for i in 1:npatch
        norm_vec .= 0
        for j in 1:npatch
            for l in 1:nlevels
                if Codes[i, l] == Codes[j, l]
                    @views norm_vec[l:nlevels] .+= 1
                    level_vec[j] = levels[l]
                    break
                end
            end
      end
        for l in 1:nlevels
            for j in 1:npatch
               if level_vec[j] ∈ levels[1:l]
                   HMM[i, j] += (ξ[l]/norm_vec[l])
                end
            end
        end
    end
    HMM
end


function ODMixingMatrix(SA, δᴴ)
    OD = CSV.read("data/GMelb_$(SA)_URxPOW_2016.csv", DataFrame)
    OD_mat = Matrix(OD[:,2:end])

    ## remove diag
    for i in 1:size(OD_mat,1)
        OD_mat[i,i] = 0
    end
    OD_mat = convert(Array{Float64}, OD_mat)

    #normalise rows
    for i in eachslice(OD_mat, dims = 1)
        sumi = sum(i)
        if sumi == 0
            i .= 0.0
        else
            i ./= sumi
        end
    end

    δᴬ = 1-δᴴ

    for i in 1:size(OD_mat,1)
        for j in 1:size(OD_mat,2)
            if i == j
                OD_mat[i,j] = δᴴ
            else
                OD_mat[i,j] = δᴬ*OD_mat[i,j]
            end
        end
    end
    OD_mat
end

function HPMixingMatrix(llCodepop::DataFrame, ξ::Vector{Float64}, μ::Float64)
    Codes = expandCodes(llCodepop[:,1])
    
    Popns = Vector{Dict{String, Int}}()

    for i in 1:ncol(Codes)
        Codes |> x -> 
            leftjoin(x, llCodepop, on = (Symbol(names(Codes)[1]) => Symbol(names(llCodepop)[1]))) |> x ->
            groupby(x, i) |> x ->
            combine(x, Symbol(names(llCodepop)[2]) => sum => :Pop) |> x ->
            Pair.(x[:,1], x.Pop) |> x ->
            Dict(x) |> x ->
            push!(Popns, x)
    end
    
    SAMM = [SAMMij(i, j, ξ, Codes, Popns) for i in eachindex(Codes[:,1]), j in eachindex(Codes[:,1])]

    PPMM = [PPMMij(i, j, llCodepop[:,2]) for i in eachindex(Codes[:,1]), j in eachindex(Codes[:,1])]

    MMM = ((μ*SAMM) + ((1-μ)*PPMM))

    return MMM
end

function PPMMij(i, j, popns)
    popns[j]/sum(popns)
end

function SAMMij(i, j, ξ, Codes, Popns)
    for l in eachindex(ξ)
        if Codes[i, l] == Codes[j, l]

            l == 1 && return ξ[1]

            ll_pop = Popns[1]|> x -> get(x, Codes[j, 1], false) 
            ul_pop = get(Popns[l], Codes[i, l], false) - get(Popns[l-1], Codes[i, l-1], false)
            return ξ[l] * (ll_pop/ul_pop)
        end
    end
end


function  expandCodes(codes::Vector{String})
    lowerCodes = string.(codes)

    Code_DF = DataFrame(getSARegion(lowerCodes[1])  => lowerCodes)
        
    while getSARegion(lowerCodes[1]) != "S/T"
        upperCodes = DataFrame(new = getNextLevelUpCode.(lowerCodes))
        rename!(upperCodes, 1 => getSARegion(upperCodes[1,1]))
        Code_DF = hcat(Code_DF, upperCodes)
        lowerCodes = upperCodes[:,1]
    end
    Code_DF
end


function getSARegion(SAcode)
    x = SAcode|>length
    x == 1 ? "S/T" :
    x == 3 ? "SA4" :
    x == 5 ? "SA3" :
    x == 9 ? "SA2" :
    x == 11 ? "SA1" :
    "Invalid Code"
end

function getfirstn(string::String, n::Int64)
    string[1:n]
end

function getNextLevelUpCode(Code)
    SA_level = getSARegion(Code)
    if SA_level == "S/T"
        throw(DomainError(SA_level, "State/Territory level (no higher level)"))
    elseif SA_level== "SA2"
        upcode = Code[1:(length(Code)-4)]
    else
        upcode = Code[1:(length(Code)-2)]
    end

    return upcode
end