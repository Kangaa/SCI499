

module MixingMatrices

using Shapefile
using Tables
using DataFrames

export SpatialMixingMatrix

function SpatialMixingMatrix(Shp, intra = 0.9) 
    
    Codes = Shp |> getCodes
    npatch = nrow(Codes)
    
    MM = fill(1.0, npatch, npatch)
    
    intraregioncoefficient = intra
    interregioncoefficient = 1.0-intraregioncoefficient
    
    for i in eachindex(Codes[:,1])
        for j in eachindex(Codes[:, 1])
            for l in 1:(ncol(Codes)-1)
                if Codes[i,l] == Codes[j,l]
                    nLLinUL = (Codes[i,l] .== Codes[:, l])|> sum
                    coef = (intraregioncoefficient/nLLinUL)
                    MM[i,j] *= coef
                    break
                else
                    if l == (ncol(Codes)-1) 
                        nLLoutUL = (Codes[i,l] .!= Codes[:, l])|> sum
                        MM[i,j] *= interregioncoefficient/nLLoutUL
                    else
                        coef = interregioncoefficient
                        MM[i,j] *= coef
                    end
                end
            end
        end
    end
    
    return MM
end

function  getCodes(shapetable)
    shapefields = shapetable |>
    Tables.columnnames .|>
    String 
    lowerCodes = Tables.getcolumn(shapetable, findfirst(endswith.(shapefields, "CODE21")))

    Code_DF = DataFrame(ll = lowerCodes)

    rename!(Code_DF, 1 => MixingMatrices.getSARegion(Code_DF[1,1]))
    
    while MixingMatrices.getSARegion(lowerCodes[1]) != "S/T"
        upperCodes = DataFrame(new = MixingMatrices.getNextLevelUpCode.(lowerCodes))
        rename!(upperCodes, 1 => MixingMatrices.getSARegion(upperCodes[1,1]))
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
end