### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d2feb6a9-ac62-4658-8493-38cf33b96ac1
using DrWatson

# ╔═╡ 6f465bf0-1277-421e-8d33-24f3794be148
@quickactivate :SCI499

# ╔═╡ b54f10e8-633c-4b9c-9ac3-53bf315ae2e9
using DataFramesMeta

# ╔═╡ 6bbf27b9-e0d1-4a42-9f00-b5e3566d24b2
using Shapefile, CSV

# ╔═╡ aff389f2-c2b6-4246-a107-0989ca88a330
using PlutoUI

# ╔═╡ 39692a93-d97c-4f92-b5aa-5397ef931f44
Gmelb_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp")) |> x ->
    Tables.subset(x, x.GCC_NAME21 .== "Greater Melbourne")|>
    DataFrame|> x-> 
    rename(x,
        "SA2_NAME21" => "SA2_NAME", 
        "SA2_CODE21" => "SA2_CODE" )


# ╔═╡ 502acc84-04ef-47fd-b18a-81dfeee77774
Vic_pop = CSV.read(datadir("VicPop21.csv"), DataFrame) 


# ╔═╡ 4720484c-709c-48e5-8adf-f2225be14f3a

Gmelb_pop = Vic_pop|> x -> 
    transform(x, "SA2 code" => (y -> string.(y)) => "SA2_CODE") |> x->
    rename(x, "2021" => "popn")


# ╔═╡ 7046cd20-e3fb-4d85-91cf-af7cb147298a

Gmelb =leftjoin(Gmelb_SA2_SHP, Gmelb_pop[:,["SA2_CODE","popn" ]], on = "SA2_CODE")



# ╔═╡ b9c8ae4b-e92e-4cf3-8d74-b5e0580720cc
codes = Gmelb |>
    names |>
    (x-> contains.(x, "CODE")) |>
    (x-> select(Gmelb, x))|>
    (x -> select(x, 1:3))

# ╔═╡ e6e474a7-f421-4a0d-9528-f31d38d4eaba
SMM = MixingMatrix(SpatialMixingMatrix(codes, 0.5), Gmelb.SA2_NAME, Gmelb.popn).mm

# ╔═╡ 084b2183-add2-4bb8-b50e-71651f90441a
## flat mixing matrix(mixing proportional to population)
MM = fill(1.0, nrow(Gmelb), nrow(Gmelb))

# ╔═╡ 51a59e35-daa4-4cb1-9e75-b5ba6897f03c
for i in eachindex(Gmelb.popn)
    for j in eachindex(Gmelb.popn)
        MM[i,j] = Gmelb.popn[j]/sum(Gmelb.popn)
    end
end

# ╔═╡ 25bab25b-986d-4ad8-9f73-01659a859a6c
@bind μ html"<input type=range 0:1>"

# ╔═╡ aa920b12-94b4-4c35-a673-a693cae97e1e
MMM = (μ * )

# ╔═╡ f554a50b-0ef7-411b-a123-18176b6a7a7f


# ╔═╡ Cell order:
# ╠═d2feb6a9-ac62-4658-8493-38cf33b96ac1
# ╠═6f465bf0-1277-421e-8d33-24f3794be148
# ╠═b54f10e8-633c-4b9c-9ac3-53bf315ae2e9
# ╠═6bbf27b9-e0d1-4a42-9f00-b5e3566d24b2
# ╠═39692a93-d97c-4f92-b5aa-5397ef931f44
# ╠═502acc84-04ef-47fd-b18a-81dfeee77774
# ╠═4720484c-709c-48e5-8adf-f2225be14f3a
# ╠═7046cd20-e3fb-4d85-91cf-af7cb147298a
# ╠═b9c8ae4b-e92e-4cf3-8d74-b5e0580720cc
# ╠═e6e474a7-f421-4a0d-9528-f31d38d4eaba
# ╠═084b2183-add2-4bb8-b50e-71651f90441a
# ╠═51a59e35-daa4-4cb1-9e75-b5ba6897f03c
# ╠═aff389f2-c2b6-4246-a107-0989ca88a330
# ╠═25bab25b-986d-4ad8-9f73-01659a859a6c
# ╠═aa920b12-94b4-4c35-a673-a693cae97e1e
# ╠═f554a50b-0ef7-411b-a123-18176b6a7a7f
