import ArchGDAL as AG

AUS_SA3_SHP = AG.read(datadir("ASGS_GDA2020/SA3_2021_AUST_SHP_GDA2020/SA3_2021_AUST_GDA2020.shp"))
layer = AG.getlayer(AUS_SA3_SHP, 0)

featuredefn = AG.layerdefn(layer)
geomdefn = AG.getgeomdefn(featuredefn, 0)

geom = ArchGDAL.getfeature(layer, 45) do feature
    AG.getgeom(feature, 0)
end


AG.simplify(geom, 0.000001)|>plot

using GeoMakie

Makie.plot(GMelb_SA3_SHP.geometry)