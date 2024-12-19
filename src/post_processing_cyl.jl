module post_process
using Plots
using CSV, DataFrames

# msize = 3
# KRs = 0.1:0.1:2.5
# mass_plots = plot(xlabel="k [-]", ylabel="A₃₃ [-]")
# damping_plots = plot(xlabel="k [-]", ylabel="B₃₃ [-]")
# data1 = CSV.read("data/exp_raw/results_rect_agfem2_order2.csv",DataFrame)
# data2 = CSV.read("data/exp_raw/results_rect_cutfem2_order2.csv",DataFrame)
# data3 = CSV.read("data/exp_raw/results_rect_sbm2_order2.csv",DataFrame)

# # reference data
# # HR000
# data_m00 = CSV.read("data/exp_pro/ref_cyl_meHR00.csv",DataFrame)
# data_d00 = CSV.read("data/exp_pro/ref_cyl_deHR00.csv",DataFrame)

# # HR0342
# data_m0342 = CSV.read("data/exp_pro/ref_cyl_meHR0342.csv",DataFrame)
# data_d0342 = CSV.read("data/exp_pro/ref_cyl_deHR0342.csv",DataFrame)

# # HR0643
# data_m0643 = CSV.read("data/exp_pro/ref_cyl_meHR0643.csv",DataFrame)
# data_d0643 = CSV.read("data/exp_pro/ref_cyl_deHR0643.csv",DataFrame)

# # HR0809
# data_m0809 = CSV.read("data/exp_pro/ref_cyl_meHR0809.csv",DataFrame)
# data_d0809 = CSV.read("data/exp_pro/ref_cyl_deHR0809.csv",DataFrame)

# # HR0906
# data_m0906 = CSV.read("data/exp_pro/ref_cyl_meHR0906.csv",DataFrame)
# data_d0906 = CSV.read("data/exp_pro/ref_cyl_deHR0906.csv",DataFrame)

function plot_data(name,i)
data1 = CSV.read("data/exp_raw/cyl/cutfemHR"*name*".csv",DataFrame)
data2 = CSV.read("data/exp_raw/cyl/agfemHR"*name*".csv",DataFrame)
data3 = CSV.read("data/exp_raw/cyl/sbmHR"*name*".csv",DataFrame)
data4 = CSV.read("data/exp_raw/cyl/sbmp2HR"*name*".csv",DataFrame)
refdatam = CSV.read("data/exp_pro/ref_cyl_me_HR"*name*".csv",DataFrame)
refdatad = CSV.read("data/exp_pro/ref_cyl_de_HR"*name*".csv",DataFrame)
msize = 2
KRs = 0.1:0.1:2.0
mass_plots = plot(xlabel="k̄ [-]", ylabel="Ā₃₃ [-]")
damping_plots = plot(xlabel="k̄ [-]", ylabel="B̄₃₃ [-]")
fp1=[0.9770603882256099,0.9806821686382569,0.9797822149354175,0.9787729298081969,0.9757585800325395]
fp2=[0.9230451174596281,0.9237240546005773,0.9405469676990527,0.9327956430122252,0.9359347898024623]

plot!(mass_plots, KRs, data1[!,1], label="CutFEM",linecolor="#332288")#,markershape=:circle)
plot!(mass_plots, KRs, data2[!,1], label="AgFEM",linecolor="#CC6677")#,markershape=:cross)
plot!(mass_plots, KRs, fp1[i].*data3[!,1], label="SBM p=1",linecolor="#117733")#,markershape=:diamond)
plot!(mass_plots, KRs, fp2[i].*data4[!,1], label="SBM p=2",linecolor="#117733",linestyle=:dash,markershape=:utriangle,markercolor="#117733")#,markershape=:diamond)
# scatter!(mass_plots, data_me[!,1], data_me[!,2], label="experiments",markershape=:cross,markercolor=:black)#,markershape=:diamond)
# plot!(mass_plots,data_mhpc[!,1],data_mhpc[!,2],label="h/R=0.0",linecolor="#882255", linestyle=:dash,markershape=:circle,markercolor="#882255",markersize=msize)
# plot!(mass_plots,data_mxfem1[!,1],data_mxfem1[!,2],label="h/R=0.342",linecolor="#88CCEE", linestyle=:dash,markershape=:utriangle,markercolor="#88CCEE",markersize=msize)
# plot!(mass_plots,data_mxfem2[!,1],data_mxfem2[!,2],label="h/R=0.643",linecolor="#DDCC77", linestyle=:dash,markershape=:dtriangle,markercolor="#DDCC77",markersize=msize)
# plot!(mass_plots,data_mfem1[!,1],data_mfem1[!,2],label="h/R=0.809",linecolor="#44AA99", linestyle=:dash,markershape=:+,markercolor="#44AA99",markersize=msize)

plot!(mass_plots,refdatam[!,1],refdatam[!,2],label="h/R=0."*name[2:end],linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)

plot!(damping_plots, KRs, data1[!,2], label="CutFEM",linecolor="#332288")#,markershape=:circle)
plot!(damping_plots, KRs, data2[!,2], label="AgFEM",linecolor="#CC6677")#,markershape=:cross)
plot!(damping_plots, KRs, fp1[i].*data3[!,2], label="SBM p=1",linecolor="#117733")#,markershape=:diamond)
plot!(damping_plots, KRs, fp2[i].*data4[!,2], label="SBM p=2",linecolor="#117733",linestyle=:dash,markershape=:utriangle,markercolor="#117733")#,markershape=:diamond)

# scatter!(damping_plots, data_de[!,1], data_de[!,2], label="experiments",markershape=:cross,markercolor=:black)#,markershape=:diamond)
# plot!(mass_plots,data_dhpc[!,1],data_dhpc[!,2],label="h/R=0.0",linecolor="#882255", linestyle=:dash,markershape=:circle,markercolor="#882255",markersize=msize)
# plot!(mass_plots,data_dxfem1[!,1],data_dxfem1[!,2],label="h/R=0.342",linecolor="#88CCEE", linestyle=:dash,markershape=:utriangle,markercolor="#88CCEE",markersize=msize)
# plot!(mass_plots,data_dxfem2[!,1],data_dxfem2[!,2],label="h/R=0.643",linecolor="#DDCC77", linestyle=:dash,markershape=:dtriangle,markercolor="#DDCC77",markersize=msize)
# plot!(mass_plots,data_dfem1[!,1],data_dfem1[!,2],label="h/R=0.809",linecolor="#44AA99", linestyle=:dash,markershape=:+,markercolor="#44AA99",markersize=msize)
plot!(damping_plots,refdatad[!,1],refdatad[!,2],label="h/R=0."*name[2:end],linecolor="#AA4499", linestyle=:dash,markershape=:x,markercolor="#AA4499",markersize=msize)


display(mass_plots)
display(damping_plots)

savefig(mass_plots, "A_cyl"*name*".png")
savefig(damping_plots, "B_cyl"*name*".png")

end # function

names = ["0000","0342","0643","0809","0906"]

for (i,name) in enumerate(names)
    plot_data(name,i)
end # for

end # module