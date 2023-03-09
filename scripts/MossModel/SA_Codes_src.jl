using DrWatson
@quickactivate


module SAcodes
using Tables

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