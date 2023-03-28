

module MixingMatrices

using Shapefile
using Tables
using Meshes

export Spatial_MixingMatrix

function shape2mesh(geom::Shapefile.Polygon)

    if length(geom.parts) > 1
        ringstart = last(geom.parts) + 1
    else
        ringstart = 1
    end

    point_vec = Vector{Tuple}()

    for i in ringstart:length(geom.points)
        x = geom.points[i].x
        y = geom.points[i].y
        push!(point_vec,(x,y))
    end
    poly = point_vec |> PolyArea 
end

function Spatial_MixingMatrix(Shp, Mixing  = "none")
    Codes = Shp |> getCodes
    LLSA_length = Codes[1]|>length 

    UL_Codes = Codes .|> getNextLevelUpCode
    ULSA_length = UL_Codes[1]|>length

    npatch = length(Codes)
    
    MM = Matrix(undef, npatch, npatch)
    
    intraregioncoefficient = 0.999
    interregioncoefficient = 1.0-intraregioncoefficient

    for (i, origin) in enumerate(Codes)
        originUL = getfirst(origin, ULSA_length)
        nLLinUL = (UL_Codes .== originUL )|> sum
        nLLnotinUL= npatch - nLLinUL
        for (j, destination) in enumerate(Codes)
            if getfirst(origin, ULSA_length) == getfirst(destination, ULSA_length)
                MM[i,j] = intraregioncoefficient/nLLinUL
            else
                MM[i,j] = interregioncoefficient/nLLnotinUL
            end
        end
    end
    return MM
end

function  getCodes(shapetable)
    shapefields = shapetable |>
    Tables.columnnames .|>
    String 
    lowerCodeFieldIndex = findfirst(endswith.(shapefields, "CODE21"))
    lowerCodeFieldName = shapefields[lowerCodeFieldIndex]
    lowerCodes = Tables.getcolumn(shapetable, Symbol(lowerCodeFieldName))

    lowerCodes
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

function getfirst(string::String, n::Int64)
    string[1:n]
end

function getNextLevelUpCode(Code)
    r = getSARegion(Code)
    if r == "S/T"
        throw(DomainError(r, "top level"))
    elseif r == "SA2"
        upcode = Code[1:(length(Code)-4)]
    else
        upcode = Code[1:(length(Code)-2)]
    end
end


function nregionsinSA4(SACode)
    (SA4_Codes .== SACode[1:3] )|> sum
end

end