module post_process
using Plots
using CSV, DataFrames

msize = 3
KRs = 0.1:0.1:2.5
mass_plots = plot(xlabel="k̄ [-]", ylabel="Ā₃₃ [-]")
damping_plots = plot(xlabel="k̄ [-]", ylabel="B̄₃₃ [-]")
data1 = CSV.read("data/exp_raw/results_rect_agfem2_order2.csv",DataFrame)
data2 = CSV.read("data/exp_raw/results_rect_cutfem2_order2.csv",DataFrame)
data3 = CSV.read("data/exp_raw/results_rect_sbm2_order2.csv",DataFrame)

# reference data
# experimental
data_me = CSV.read("data/exp_pro/ref_rect_me.csv",DataFrame)
data_de = CSV.read("data/exp_pro/ref_rect_de.csv",DataFrame)

# HPC
data_mhpc = CSV.read("data/exp_pro/ref_rect_mhpc.csv",DataFrame)
data_dhpc = CSV.read("data/exp_pro/ref_rect_dhpc.csv",DataFrame)

# linear FEM
data_mfem1 = CSV.read("data/exp_pro/ref_rect_mfem1.csv",DataFrame)
data_dfem1 = CSV.read("data/exp_pro/ref_rect_dfem1.csv",DataFrame)

# quadratic FEM
data_mfem2 = CSV.read("data/exp_pro/ref_rect_mfem2.csv",DataFrame)
data_dfem2 = CSV.read("data/exp_pro/ref_rect_dfem2.csv",DataFrame)

# linear XFEM
data_mxfem1 = CSV.read("data/exp_pro/ref_rect_mxfem1.csv",DataFrame)
data_dxfem1 = CSV.read("data/exp_pro/ref_rect_dxfem1.csv",DataFrame)

# quadratic XFEM
data_mxfem2 = CSV.read("data/exp_pro/ref_rect_mxfem2.csv",DataFrame)
data_dxfem2 = CSV.read("data/exp_pro/ref_rect_dxfem2.csv",DataFrame)


plot!(mass_plots, KRs, data2[!,1], label="CutFEM p=2",linecolor="#332288")#,markershape=:circle)
plot!(mass_plots, KRs, data1[!,1], label="AgFEM p=2",linecolor="#CC6677")#,markershape=:cross)
plot!(mass_plots, KRs, data3[!,1], label="SBM p=2",linecolor="#117733")#,markershape=:diamond)
scatter!(mass_plots, data_me[!,1], data_me[!,2], label="experiments",markershape=:cross,markercolor=:black)#,markershape=:diamond)
plot!(mass_plots,data_mfem1[!,1],data_mfem1[!,2],label="FEM p=1",linecolor="#44AA99", linestyle=:dash,markershape=:+,markercolor="#44AA99",markersize=msize)
plot!(mass_plots,data_mfem2[!,1],data_mfem2[!,2],label="FEM p=2",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
plot!(mass_plots,data_mxfem1[!,1],data_mxfem1[!,2],label="XFEM p=1",linecolor="#88CCEE", linestyle=:dash,markershape=:utriangle,markercolor="#88CCEE",markersize=msize)
plot!(mass_plots,data_mxfem2[!,1],data_mxfem2[!,2],label="XFEM p=2",linecolor="#DDCC77", linestyle=:dash,markershape=:dtriangle,markercolor="#DDCC77",markersize=msize)
plot!(mass_plots,data_mhpc[!,1],data_mhpc[!,2],label="HPC",linecolor="#882255", linestyle=:dash,markershape=:circle,markercolor="#882255",markersize=msize)

plot!(damping_plots, KRs, data2[!,2], label="CutFEM p=2",linecolor="#332288")#,markershape=:circle)
plot!(damping_plots, KRs, data1[!,2], label="AgFEM p=2",linecolor="#CC6677")#,markershape=:cross)
plot!(damping_plots, KRs, data3[!,2], label="SBM p=2",linecolor="#117733")#,markershape=:diamond)
scatter!(damping_plots, data_de[!,1], data_de[!,2], label="experiments",markershape=:cross,markercolor=:black)#,markershape=:diamond)
plot!(damping_plots,data_dfem1[!,1],data_dfem1[!,2],label="FEM p=1",linecolor="#44AA99", linestyle=:dash,markershape=:+,markercolor="#44AA99",markersize=msize)
plot!(damping_plots,data_dfem2[!,1],data_dfem2[!,2],label="FEM p=2",linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)
plot!(damping_plots,data_dxfem1[!,1],data_dxfem1[!,2],label="XFEM p=1",linecolor="#88CCEE", linestyle=:dash,markershape=:utriangle,markercolor="#88CCEE",markersize=msize)
plot!(damping_plots,data_dxfem2[!,1],data_dxfem2[!,2],label="XFEM p=2",linecolor="#DDCC77", linestyle=:dash,markershape=:dtriangle,markercolor="#DDCC77",markersize=msize)
plot!(damping_plots,data_dhpc[!,1],data_dhpc[!,2],label="HPC",linecolor="#882255", linestyle=:dash,markershape=:circle,markercolor="#882255",markersize=msize)


display(mass_plots)
display(damping_plots)

savefig(mass_plots, "Added_Mass_vs_KR.png")
savefig(damping_plots, "Added_Damping_vs_KR.png")

end # module