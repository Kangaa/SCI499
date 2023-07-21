#Functions for constructing mixing matrices

struct MixingMatrix
    mm::Array{Real, 2}
    names::Vector{String}
    popns::Array{Int, 1}
end

function SpatialMixingMatrix(llCodepop, ξ, μ = 0.5)
    Codes = expandCodes(llCodepop[:,1])
    
    Popns = Vector{Dict{String, Int}}()

    for i in 1:ncol(Codes)-1
        Popns[i] = Codes |> x -> 
            leftjoin(x, llCodepop, on = (:SA2 => Symbol("SA2 code"))) |> x ->
            groupby(x, i) |> x ->
            combine(x, :Pop => sum => :Pop) |> x ->
            Pair.(x[:,i], x.Pop) 
    end
    
    SAMM = [SAMMij(i, j, ξ, Codes, Popns) for i in eachindex(Codes[:,1]), j in eachindex(Codes[:,1])]

    PPMM = [PPMMij(i, j, Popns) for i in eachindex(Codes[:,1]), j in eachindex(Codes[:,1])]

    MMM = ((μ*MM) + ((1-μ)*PPMM))

    return MMM
end

function SAMMij(i, j, ξ, codes, popns)
    for l in eachindex(ξ)
        if codes[i, l] == codes[j, l]
            l == 1 && return ξ[l] 
            ll_pop = popns[1]|> x -> get(x, codes[j, 1], false) 
            ul_pop = get(popns[l], codes[i, l], false) - get(popns[l-1], codes[i, l-1], false)
            return ξ[l] * (ll_pop/ul_pop)
        end
    end
end

PPMMij(i, j, popns) = popns[j]/sum(popns)

function  expandCodes(shapetable::Shapefile.Table)
    shapefields = shapetable |>
    Tables.columnnames .|>
    String 
    lowerCodes = Tables.getcolumn(shapetable, findfirst(endswith.(shapefields, "CODE21")))

    Code_DF = DataFrame(ll = lowerCodes)

    rename!(Code_DF, 1 => getSARegion(Code_DF[1,1]))
    
    while getSARegion(lowerCodes[1]) != "S/T"
        upperCodes = DataFrame(new = getNextLevelUpCode.(lowerCodes))
        rename!(upperCodes, 1 => getSARegion(upperCodes[1,1]))
        Code_DF = hcat(Code_DF, upperCodes)
        lowerCodes = upperCodes[:,1]
    end
    Code_DF
end

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