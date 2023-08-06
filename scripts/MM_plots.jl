include("src/MixingMatrices.jl")

using CSV
using DataFrames
using CairoMakie


SA2_data= CSV.read("data/GmelbSA2Pop21.csv", DataFrame)
SA2_CodePop = DataFrame(Codes = string.(SA2_data[:, 1]), Pop = SA2_data[:,3])
MM2 = SpatialMixingMatrix(SA2_CodePop, [1/4,1/4,1/4,1/4], 1.0)

mm2_plot = Figure();
mm2_axis, mm2_hm = Makie.heatmap(mm2_plot[1,1],  log10.(transpose(MM2)), colormap = :plasma)
colorbar = Colorbar(mm2_plot[1,2], mm2_hm)
mm2_plot
save("MMSA2.png", mm2_plot)

SA3_codePop = CSV.read("data/GmelbSA3Pop21.csv", DataFrame)
transform!(SA3_codePop, :SA3_Code => ByRow(x -> string(x)) => :SA3_Code)
CodePop = DataFrame(Codes = SA3_codePop[:, 1], Pop = SA3_codePop[:,2])
MM3 = SpatialMixingMatrix(CodePop, [1/2,1/4,1/4],1.0 )

mm3_plot = Figure();
mm3_axis, mm3_hm = Makie.heatmap(mm3_plot[1,1],  log10.(transpose(MM3)), colormap = :plasma)
colorbar = Colorbar(mm3_plot[1,2], mm3_hm)
mm3_plot
Makie.save("MMSA3.png", mm3_plot)


SA4_CodePop = expandCodes(SA3_codePop[:, 1]) |> x ->
leftjoin(x, SA3_codePop, on = (:SA3 => :SA3_Code)) |> x ->
groupby(x, :SA4)|> x -> 
combine(x, :SA3_pop => sum => :Pop)

MM4 = SpatialMixingMatrix(SA4_CodePop, [3/4, 1/4], 1.0 )
mm4_plot = Figure()
ax, hm = Makie.heatmap(mm4_plot[1,1], log10.(MM4'))
cb = Colorbar(mm4_plot[1, 2], hm)
Makie.save("MMSA4.png", mm4_plot)



Moss_mm = CSV.read("data/Moss_new/Moss_MM_5.csv", DataFrame)|> x -> Matrix(x)

for i in 1:size(Moss_mm, 1)
    for j in 1:size(Moss_mm, 2)
        if  iszero(Moss_mm[i,j])
            Moss_mm[i,j] = 1e-4
        end
    end
end


fig = Figure();
ax, hm = Makie.heatmap(fig[1,1], log10.(Moss_mm'), colormap = :plasma)
cb = Colorbar(fig[1, 2],hm)
text!("μ = 0.5");
cb.scale = x ->2.0^x
maximum(Moss_mm)







fig, ax, hm = Makie.heatmap(log10.(MM2'), colormap = :plasma);

colorbar = Colorbar(fig[1, 2], hm)
nframes = 120;
framerate = 60;
μ_iterator = [range(1, 0, length = nframes); range(0, 1, length = nframes)]
record(fig, "test.mp4", μ_iterator; framerate = framerate) do i
    mm = SpatialMixingMatrix(SA2_CodePop, [1/4,1/4,1/4,1/4], i)
    lhm = log10.(mm')
    hm.values = lhm
 end