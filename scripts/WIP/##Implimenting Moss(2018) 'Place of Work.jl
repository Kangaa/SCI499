##Implimenting Moss(2018) 'Place of Work' (POW) Mixing Matrix

## import the OD data
using CSV
using DataFrames
using LinearAlgebra

    OD = CSV.read("data/GMelb_SA2_URxPOW_2016.csv", DataFrame)
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

OD[:,1]

Hierarchy = CSV.read("data/MB_2021_AUST.csv", DataFrame)
OD = CSV.read("data/GMelb_SA2_URxPOW_2016.csv", DataFrame)

Grouped_DF = subset(Hierarchy, :GCCSA_NAME_2021 => y -> y .== "Greater Melbourne") |>
 	x -> groupby(x, :SA2_NAME_2021) |>
    x -> combine(x,
        :SA3_NAME_2021 => unique => :SA3_NAME,
        :SA4_NAME_2021 => unique => :SA4_NAME)

stacked_OD = OD |>
    x -> stack(x, Not(1)) |>
    x -> rename(x, 1 => :UR, 2 => :POW, 3 => :N)

leftjoin(stacked_OD, Grouped_DF, on = :UR => :SA2_NAME_2021, renamecols = "" => "_UR") |>
    x -> leftjoin(x, Grouped_DF, on = :POW => :SA2_NAME_2021, renamecols = "" => "_POW") |>
    x -> groupby(x, [:SA3_NAME_UR, :SA3_NAME_POW])     |>
    x -> combine(x, :N => sum => :N) |>
    x -> unstack(x, :SA3_NAME_UR, :SA3_NAME_POW, :N, allowmissing = true) |>
    x -> x[1:40, 1:41]|>
    x -> CSV.write("data/GMelb_SA3_URxPOW_2016.csv", x)

leftjoin(stacked_OD, Grouped_DF, on = :UR => :SA2_NAME_2021, renamecols = "" => "_UR") |>
    x -> leftjoin(x, Grouped_DF, on = :POW => :SA2_NAME_2021, renamecols = "" => "_POW") |>
    x -> groupby(x, [:SA4_NAME_UR, :SA4_NAME_POW])     |>
    x -> combine(x, :N => sum => :N) |>
    x -> unstack(x, :SA4_NAME_UR, :SA4_NAME_POW, :N, allowmissing = true) |>
    x -> x[1:9, 1:10]|>
    x -> CSV.write("data/GMelb_SA4_URxPOW_2016.csv", x)