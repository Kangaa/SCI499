using CSV

struct OdMat
    m::Array{Float64, 2}
    names::Vector{String}
    popns::Array{Int, 1}
end

function read_popns(path::AbstractString, names::Vector{String})::Vector{Int}
    reader = CSV.read(path, delim=' ', header=false)
    popn_tbl = Dict{String,Int}()
    for row in eachrow(reader)
        name = row[2]
        if haskey(popn_tbl, name)
            error("Multiple populations for $name")
        else
            value = parse(Int, row[3])
            popn_tbl[name] = value
        end
    end
    popns = Vector{Int}()
    for name in names
        if haskey(popn_tbl, name)
            push!(popns, popn_tbl[name])
        else
            error("No population defined for $name")
        end
    end
    return popns
end

function read_OdMat(mm_path::AbstractString, popn_path::AbstractString)::OdMat
    reader = CSV.read(mm_path, DataFrame)
    names = Vector{String}()
    values = Vector{Float64}()
    for row in eachrow(reader)
        push!(names, "$row[1]")
        for val in row[2:end]
            push!(values, val)
        end
    end
    
    n = length(names)
    mat = reshape(values, (n, n))
    popns = read_popns(popn_path, names)
    
     OdMat(mat, names, popns)
end