using DrWatson
@quickactivate "SCI499"

using Shapefile
using Tables: subset
using Meshes
using MeshViz
using Makie

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

