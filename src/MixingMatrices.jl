#Functions for constructing mixing matrices

struct MixingMatrix
    mm::Array{Real, 2}
    names::Vector{String}
    popns::Array{Int, 1}
end

function SpatialMixingMatrix(llCodepop, ξ, μ = 0.5)
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

PPMMij(i, j, popns) = popns[j]/sum(popns)


function  expandCodes(codes::Vector)
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