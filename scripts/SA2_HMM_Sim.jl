#Hierarchical Mixing matrix SA2 Greater Melbourne region sim

using DrWatson 

@quickactivate :SCI499

AUS_SA3_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA3_2021_AUST_SHP_GDA2020/SA3_2021_AUST_GDA2020.shp"))
GMelb_SA3_SHP = Tables.subset(AUS_SA3_SHP, AUS_SA3_SHP.GCC_NAME21 .== "Greater Melbourne")


SA3_Mixmat = MixingMatrix(SpatialMixingMatrix(GMelb_SA3_SHP), GMelb_SA3_SHP.SA3_NAME21, fill(70000, size(GMelb_SA3_SHP.SA3_NAME21)))
SA3_params = Parameters(0.5,0.1,0.1, SA3_Mixmat)


@time ev_log = run_sim(SA3_params, 10)
#1047.119495 seconds (42.83 G allocations: 685.346 GiB, 5.36% gc time, 0.04% compilation time)


using CSV
CSV.write(datadir("sims\\GMelb_SA3_1000_10_0.5_0.1_0.1_0.9.csv"), ev_log)


AUS_SA2_SHP = Shapefile.Table(datadir("ASGS_GDA2020/SA2_2021_AUST_SHP_GDA2020/SA2_2021_AUST_GDA2020.shp"))
GMelb_SA2_SHP = Tables.subset(AUS_SA2_SHP, AUS_SA2_SHP.GCC_NAME21 .== "Greater Melbourne")


SA2_Mixmat = MixingMatrix(SpatialMixingMatrix(GMelb_SA2_SHP), GMelb_SA2_SHP.SA2_NAME21, fill(10000, size(GMelb_SA2_SHP.SA2_NAME21)))
SA2_params = Parameters(0.5,0.1,0.1, SA2_Mixmat)


@time ev_log = run_sim(SA2_params, 10)