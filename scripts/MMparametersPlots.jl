## Hierarchical mixing matrices parameter variation

using DrWatson 
@quickactivate :SCI499

AUS_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp"))
GMelb_SA2_SHP = Tables.subset(AUS_SA2_SHP, AUS_SA2_SHP.GCC_NAME21 .== "Greater Melbourne")

deltas = zip(Iterators.cycle("δ"), [0.8, 0.7, 0.5, 0.2, 0.05])

using StatsPlots

for i in deltas
    SpatialMixingMatrix(GMelb_SA2_SHP, i[2]) .|>
    log10|>
    heatmap|>
    x -> savefig(x, plotsdir("$i MMheatmap.png"))
end


