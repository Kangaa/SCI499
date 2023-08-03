struct MixMat
    m::Array{Real, 1}
    names::Vector{String}
    popns::Array{Int, 1}
    frac_self::FracSelf
    frac_cbd::Real
end

using LinearAlgebra
function scale_diag_var!(mat::Matrix{T}, frac_self::Vector{T}) where T<:AbstractFloat
    @inbounds for i in 1:size(mat, 1)
        mat[i, :] = mat[i, :] .* ((1.0 - frac_self[i]) ./ sum(mat[i, :]))
    end
    diag(mat) .= frac_self
end

function scale_diag_unif!(mat::Matrix{T}, frac_self) where T<:AbstractFloat
    frac_mix = 1.0 - frac_self
    @inbounds for i in 1:size(mat, 1)
        mat[i, :] = mat[i, :] .* (frac_mix ./ sum(mat[i, :]))
    end
    diag(mat) .= frac_self
end

function scale_diag!(mat::Matrix{T}, frac_self) where T<:AbstractFloat
    if isa(frac_self, Uniform)
        scale_diag_unif!(mat, frac_self.a)
    elseif isa(frac_self, Variable)
        scale_diag_var!(mat, frac_self.x)
    end
end

function test_scale_diag(mat::Matrix{T}, frac_self::T) where T<:AbstractFloat
    m1 = copy(mat)
    scale_diag_unif!(m1, frac_self)

    m2 = copy(mat)
    fs2 = fill(frac_self, size(mat, 1))
    scale_diag_var!(m2, fs2)

    m3 = copy(mat)
    fs3 = fill(frac_self, size(mat, 1))
    fs3[2] *= 0.9
    scale_diag_var!(m3, fs3)

    @assert all(abs.(1.0 .- sum(m1, dims=2)) .< 1e-8) "m1 rowsum"
    @assert all(abs.(1.0 .- sum(m2, dims=2)) .< 1e-8) "m2 rowsum"
    @assert all(abs.(1.0 .- sum(m3, dims=2)) .< 1e-8) "m3 rowsum"
    @assert m1 == m2 "mixing matrices differ"
    @assert m1 != m3 "mixing matrices do not differ"
end

function redistribute!(mat::Matrix{Float64}, popns::Vector{Int},
    frac_cbd::Float64, cbd_ix::Int)
    # Determine how much all other regions mix with the CBD.
    mix_w_cbd = copy(mat[:, cbd_ix])
    mix_w_cbd[cbd_ix] = 0.0

    # Weight the mixing by each region's resident population.
    mix_w_cbd .= mix_w_cbd .* Float64.(popns)

    # Normalise these values to conserve the force of infection.
    mix_w_cbd /= sum(mix_w_cbd)

        for rix in 1:size(mat, 1)
            if rix == cbd_ix
                continue
            end

        row = mat[rix, :]
        cbd_mix = row[cbd_ix]
        add_mix = [cbd_mix * (1.0 - frac_cbd) * v for v in mix_w_cbd]
        row += add_mix
        row[cbd_ix] = cbd_mix * frac_cbd
        mat[rix, :] = row
    end
end
function new_MixMat(od::OdMat, frac_self::FracSelf, frac_cbd::Float64)::MixMat
    mat = copy(od.m)
    names = copy(od.names)
    popns = copy(od.popns)
    scale_diag!(mat, frac_self)
    cbd_ix = findfirst(isequal("20604"), names)
    redistribute!(mat, popns, frac_cbd, cbd_ix)
    MixMat(mat, names, popns, frac_self, frac_cbd)
end